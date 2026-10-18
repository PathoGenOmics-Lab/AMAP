#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract protein sequences from a reference genome + annotation (GFF3 or GenBank).

Produces:
  1. {prefix}_ref_proteins.faa  — all CDS translations (trusted reference proteome)
  2. {prefix}_ref_genome.fasta  — copy of the reference genome (for reverse alignment)

Each protein header includes the gene coordinates for downstream reconciliation,
and names the sequence of the genome FASTA the gene is on:
  >GENE_ID_START-END contig=SEQUENCE
A gene across the origin of a circular sequence ends past its last base, the
way GFF3 writes it (see unroll).
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


def unroll(segments, length):
    """The segments of a feature, in the order it runs along the forward
    strand, on its sequence unrolled past the origin.

    A feature across the origin of a circular sequence runs on from its last
    base to its first: tagA on plasmid pO157 (92721 bp) of E. coli O157:H7
    Sakai is join(92527..92721,1..2502). Its header said 1-92721, the whole
    plasmid, from the lowest and the highest base its pieces reach; read as a
    span, that put the gene over every other gene of the plasmid, and 92721 bp
    is not 3 x 899 residues, so the protein lost its stop. Every segment after the one where
    the feature turns back is moved on by the sequence's length, the way GFF3
    and Bakta write such a feature (end = stop + length): 92527..95223. Spans
    and midpoints read off the header are then right, and a feature that does
    not cross the origin keeps its coordinates. Each segment is (start, end,
    anything else), 1-based.
    """
    out, shift = [], 0
    for seg in segments:
        start, end = seg[0] + shift, seg[1] + shift
        if out and start < out[-1][0]:
            shift += length
            start, end = start + length, end + length
        out.append((start, end) + tuple(seg[2:]))
    return out


def _bases(seq, start, end):
    """Bases start..end (1-based) of a sequence, read on across its origin
    for coordinates past its last base."""
    n, out, pos = len(seq), [], start - 1
    while pos < end:
        i = pos % n
        take = min(end - pos, n - i)
        out.append(str(seq[i:i + take]))
        pos += take
    return "".join(out)


def _read_genome(fasta_path):
    return SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta")) if fasta_path else {}


def _genome_name(record, genome):
    """The name the genome FASTA gives a GenBank record's sequence.

    Bakta annotates the FASTA, and find_bakta_false_positives.py compares
    each of its calls with the genes of the same sequence, so the proteins
    name their sequence the way the FASTA does: by the record's
    accession.version or LOCUS name when the FASTA has it, or else by the one
    sequence of the FASTA identical to it (of the same length, for a record
    without its sequence). None if there is no such sequence.
    """
    for name in (record.id, record.name):
        if name in genome:
            return name
    try:
        seq = str(record.seq).upper()
    except Exception:      # a CON record without its sequence
        seq = None
    same = [name for name, r in genome.items()
            if len(r) == len(record) and (seq is None or str(r.seq).upper() == seq)]
    return same[0] if len(same) == 1 else None


def _write(proteins, output_faa):
    with open(output_faa, 'w') as f:
        for header, seq, contig in proteins:
            f.write(f">{header} contig={contig}\n{seq}\n")


def extract_from_genbank(gbk_path, output_faa, fasta_path=None):
    """Extract CDS translations from a GenBank file."""
    genome = _read_genome(fasta_path)
    proteins = []
    for record in SeqIO.parse(gbk_path, "genbank"):
        contig = _genome_name(record, genome)
        if contig is None:
            if genome:
                print(f"Warning: {record.id} of the annotation is not a sequence of the genome "
                      f"FASTA, by name or by sequence; its proteins carry its own name", file=sys.stderr)
            contig = record.id
        length = len(genome[contig]) if contig in genome else len(record)
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

            # The parts of the location as the feature runs along the forward
            # strand: Biopython lists those of complement(join(...)) from the
            # 5' end, the highest first. (feat.strand was removed in Biopython
            # 1.85; the strand is read off the location.)
            parts = [(int(p.start) + 1, int(p.end)) for p in feat.location.parts]  # 0-based → 1-based
            if feat.location.strand == -1:
                parts.reverse()
            segments = unroll(parts, length)
            start, end = segments[0][0], max(e for _, e in segments)

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
                proteins.append((header, translation, contig))

    _write(proteins, output_faa)

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


def _gff3_order(segments, length):
    """The lines of a GFF3 CDS in the order it runs along the forward strand.

    GFF3 does not order them, and writes a CDS across the origin of a
    circular sequence as one line that ends past it (end = stop + length), as
    NCBI writes tagA of pO157: 92527..95223. A CDS cut in two there instead,
    one line ending on the sequence's last base and another starting on its
    first, starts after the widest gap between its lines, which is where the
    rest of the sequence lies.
    """
    segments = sorted(segments)
    if (len(segments) > 1 and segments[0][0] == 1
            and any(seg[1] == length for seg in segments)):
        gaps = [segments[i + 1][0] - segments[i][1] for i in range(len(segments) - 1)]
        i = gaps.index(max(gaps)) + 1
        segments = segments[i:] + segments[:i]
    return segments


def _exceptions(attrs, segments, strand, phase, length=0):
    """Residues the annotation reads differently from the code, by codon:
    transl_except=(pos:1547986..1547988,aa:Sec) is the selenocysteine of
    E. coli fdoG, read through a UGA that is otherwise a stop."""
    out = {}
    value = unquote(attrs.get("transl_except", ""))
    for m in re.finditer(r"pos:(complement\()?(\d+)\.\.(\d+)\)?,aa:(\w+)", value):
        lo, hi, aa = int(m.group(2)), int(m.group(3)), seq1(m.group(4))
        if not aa or aa == "X":
            continue
        # Past the origin, the unrolled segments sit one sequence length on
        pos = lo if strand != "-" else hi
        if length and not any(s <= pos <= e for s, e, _ in segments):
            lo, hi = lo + length, hi + length
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

    A CDS across the origin of a circular sequence ends past the sequence's
    last base, and was read only up to it: tagA of pO157 came out as its first
    65 residues of 898.
    """
    genome = _read_genome(fasta_path)

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

        sequence = genome[c["seq_id"]].seq
        segments = unroll(_gff3_order(c["segments"], len(sequence)), len(sequence))
        seq = Seq("".join(_bases(sequence, s, e) for s, e, _ in segments))
        first = segments[-1] if c["strand"] == "-" else segments[0]
        if c["strand"] == "-":
            seq = seq.reverse_complement()
        seq = seq[first[2]:]

        exceptions = _exceptions(attrs, segments, c["strand"], first[2], len(sequence))
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
            proteins.append((f"{gene_id}_{start}-{end}", translation, c["seq_id"]))

    _write(proteins, output_faa)

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
        n = extract_from_genbank(ann_path, output_faa, genome_path)
    elif suffix in ('.gff', '.gff3'):
        n = extract_from_gff3(ann_path, genome_path, output_faa)
    else:
        # Try GenBank first, fall back to GFF3
        try:
            n = extract_from_genbank(ann_path, output_faa, genome_path)
        except Exception:
            n = extract_from_gff3(ann_path, genome_path, output_faa)

    if n == 0:
        sys.exit("Error: No proteins extracted. Check annotation format.")

    print(f"Done. {n} reference proteins ready for AMAP.")


if __name__ == "__main__":
    main()
