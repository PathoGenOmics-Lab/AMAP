#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract protein sequences from a reference genome + annotation (GFF3 or GenBank).

Produces:
  1. {prefix}_ref_proteins.faa  — all CDS translations (trusted reference proteome)
  2. {prefix}_ref_genome.fasta  — copy of the reference genome (for reverse alignment)

Each protein header includes the gene coordinates for downstream reconciliation:
  >GENE_ID_START-END
"""

import argparse
import re
import sys
from pathlib import Path
from urllib.parse import unquote
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import seq1


def append_stop_codon(translation, start, end):
    """Match the convention the bundled H37Rv proteome uses.

    There, the coordinates in the header span the stop codon and the protein
    carries it as a trailing '*', so 3 * len(protein) equals the span. Both a
    GenBank /translation qualifier and Bio.Seq.translate(to_stop=True) leave the
    stop out, which makes the protein three bases shorter than the region its
    own header names.

    That matters because miniprot aligns the protein: a query without the stop
    ends three bases before the CDS that Bakta called, so classify_cases can
    never see an exact coordinate match. Every gene lands in case 4 instead of
    case 2, and the case counts stop meaning anything.

    Only appends when the span is exactly one codon longer than the protein,
    which is the stop-codon case and nothing else.
    """
    if not translation or translation.endswith('*'):
        return translation
    if (end - start + 1) == (len(translation) + 1) * 3:
        return translation + '*'
    return translation


def extract_from_genbank(gbk_path, output_faa):
    """Extract CDS translations from a GenBank file."""
    proteins = []
    for record in SeqIO.parse(gbk_path, "genbank"):
        for feat in record.features:
            if feat.type != "CDS":
                continue
            # Get gene ID
            gene_id = (
                feat.qualifiers.get("locus_tag", [None])[0] or
                feat.qualifiers.get("gene", [None])[0] or
                feat.qualifiers.get("protein_id", [None])[0]
            )
            if not gene_id:
                continue

            start = int(feat.location.start) + 1  # 0-based → 1-based
            end = int(feat.location.end)
            # feat.strand was removed in Biopython 1.85; read it off the location.
            strand = "+" if feat.location.strand == 1 else "-"

            # Use pre-computed translation if available
            translation = feat.qualifiers.get("translation", [None])[0]
            if not translation:
                # Extract and translate from sequence
                seq = feat.extract(record.seq)
                try:
                    translation = str(seq.translate(to_stop=True))
                except Exception:
                    continue

            translation = append_stop_codon(translation, start, end)

            if translation:
                header = f"{gene_id}_{start}-{end}"
                proteins.append((header, translation, strand))

    with open(output_faa, 'w') as f:
        for header, seq, strand in proteins:
            f.write(f">{header}\n{seq}\n")

    print(f"Extracted {len(proteins)} proteins from GenBank → {output_faa}")
    return len(proteins)


# Shorter proteins say too little to be placed on sequence identity alone; the
# same floor as the exact search in classify/exact.py.
MIN_RESIDUES = 8


def _gff3_cds(gff_path):
    """The CDS of a GFF3 file, one per ID, with all the lines that share it.

    A CDS can take several lines: a join across an insertion, or a programmed
    frameshift such as prfB's. Lines without an ID fall back on their Parent,
    and on being their own CDS.
    """
    cds, order = {}, []
    with open(gff_path) as f:
        for n, line in enumerate(f):
            if line.startswith("##FASTA"):
                break  # Stop at embedded FASTA
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split('\t')
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = {}
            for attr in parts[8].split(';'):
                if '=' in attr:
                    k, v = attr.split('=', 1)
                    attrs[k.strip()] = v.strip()
            key = (parts[0], attrs.get("ID") or attrs.get("Parent") or f"line {n}")
            if key not in cds:
                cds[key] = {"seq_id": parts[0], "strand": parts[6], "attrs": attrs, "segments": []}
                order.append(key)
            phase = int(parts[7]) if parts[7].isdigit() else 0
            cds[key]["segments"].append((int(parts[3]), int(parts[4]), phase))  # GFF3 is 1-based
    return [cds[k] for k in order]


def _exceptions(attrs, segments, strand, phase):
    """Residues the annotation reads differently from the code, by codon:
    transl_except=(pos:1547986..1547988,aa:Sec) is the selenocysteine of
    E. coli fdoG, read through a UGA that is otherwise a stop."""
    out = {}
    value = unquote(attrs.get("transl_except", ""))
    for m in re.finditer(r"pos:(complement\()?(\d+)\.\.(\d+)\)?,aa:(\w+)", value):
        lo, hi, aa = int(m.group(2)), int(m.group(3)), seq1(m.group(4))
        if not aa or aa == "X":
            continue
        offset = 0
        for s, e, _ in (segments if strand != "-" else sorted(segments, key=lambda x: -x[1])):
            if s <= (lo if strand != "-" else hi) <= e:
                offset += (lo - s) if strand != "-" else (e - hi)
                break
            offset += e - s + 1
        out[(offset - phase) // 3] = aa
    return out


def extract_from_gff3(gff_path, fasta_path, output_faa):
    """Extract CDS translations from a GFF3 file + FASTA genome.

    Each line used to be translated on its own, which cut a CDS of several
    lines into fragments, and everything under 11 residues was dropped to hide
    them, taking small proteins with it (E. coli has some of 8 and 9). The
    lines of a CDS are now joined, as GenBank's join() is, and the protein
    translated from the whole, with the bacterial code, so that a GTG or TTG
    start reads as the methionine the GenBank /translation gives.
    """
    genome = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

    proteins = []
    for c in _gff3_cds(gff_path):
        attrs = c["attrs"]
        gene_id = (
            attrs.get("locus_tag") or
            attrs.get("ID") or
            attrs.get("gene") or
            attrs.get("Name")
        )
        if not gene_id or c["seq_id"] not in genome:
            continue

        segments = sorted(c["segments"])
        seq = Seq("".join(str(genome[c["seq_id"]].seq[s - 1:e]) for s, e, _ in segments))
        first = segments[-1] if c["strand"] == "-" else segments[0]
        if c["strand"] == "-":
            seq = seq.reverse_complement()
        seq = seq[first[2]:]

        exceptions = _exceptions(attrs, segments, c["strand"], first[2])
        try:
            if exceptions:
                raise ValueError("translated codon by codon below")
            translation = str(seq.translate(table=11, cds=True))
        except Exception:
            # Not a clean CDS (a frameshift, a partial gene, a selenocysteine):
            # codon by codon, with the annotation's exceptions, to the first stop.
            aa = list(str(seq[:len(seq) - len(seq) % 3].translate(table=11)))
            for i, residue in exceptions.items():
                if 0 <= i < len(aa):
                    aa[i] = residue
            if exceptions and aa and str(seq[:3]).upper() in ("GTG", "TTG", "CTG", "ATT", "ATC", "ATA"):
                aa[0] = "M"
            translation = "".join(aa).split("*")[0]

        start, end = segments[0][0], max(e for _, e, _ in segments)
        translation = append_stop_codon(translation, start, end)

        if len(translation.rstrip('*')) >= MIN_RESIDUES:
            proteins.append((f"{gene_id}_{start}-{end}", translation))

    with open(output_faa, 'w') as f:
        for header, seq in proteins:
            f.write(f">{header}\n{seq}\n")

    print(f"Extracted {len(proteins)} proteins from GFF3 → {output_faa}")
    return len(proteins)


def main():
    parser = argparse.ArgumentParser(
        description="Extract reference proteins from GFF3/GenBank + FASTA for AMAP.")
    parser.add_argument("--ref-annotation", required=True,
                        help="Reference annotation file (GFF3 or GenBank format).")
    parser.add_argument("--ref-genome", required=True,
                        help="Reference genome FASTA file.")
    parser.add_argument("--prefix", default="custom_ref",
                        help="Output file prefix (default: custom_ref).")
    args = parser.parse_args()

    ann_path = Path(args.ref_annotation)
    genome_path = Path(args.ref_genome)
    prefix = args.prefix

    if not ann_path.exists():
        sys.exit(f"Error: Annotation file not found: {ann_path}")
    if not genome_path.exists():
        sys.exit(f"Error: Genome file not found: {genome_path}")

    output_faa = f"{prefix}_ref_proteins.faa"

    # Detect format
    suffix = ann_path.suffix.lower()
    if suffix in ('.gbk', '.gb', '.genbank', '.gbff'):
        n = extract_from_genbank(ann_path, output_faa)
    elif suffix in ('.gff', '.gff3'):
        n = extract_from_gff3(ann_path, genome_path, output_faa)
    else:
        # Try GenBank first, fall back to GFF3
        try:
            n = extract_from_genbank(ann_path, output_faa)
        except Exception:
            n = extract_from_gff3(ann_path, genome_path, output_faa)

    if n == 0:
        sys.exit("Error: No proteins extracted. Check annotation format.")

    print(f"Done. {n} reference proteins ready for AMAP.")


if __name__ == "__main__":
    main()
