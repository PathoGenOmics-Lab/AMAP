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
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


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


def extract_from_gff3(gff_path, fasta_path, output_faa):
    """Extract CDS translations from a GFF3 file + FASTA genome."""
    # Load genome
    genome = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

    proteins = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("##FASTA"):
                if line.startswith("##FASTA"):
                    break  # Stop at embedded FASTA
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != "CDS":
                continue

            seq_id = parts[0]
            start = int(parts[3])       # GFF3 is 1-based
            end = int(parts[4])
            strand = parts[6]

            # Parse attributes
            attrs = {}
            for attr in parts[8].split(';'):
                if '=' in attr:
                    k, v = attr.split('=', 1)
                    attrs[k.strip()] = v.strip()

            gene_id = (
                attrs.get("locus_tag") or
                attrs.get("ID") or
                attrs.get("gene") or
                attrs.get("Name")
            )
            if not gene_id or seq_id not in genome:
                continue

            # Extract sequence
            seq = genome[seq_id].seq[start - 1:end]
            if strand == "-":
                seq = seq.reverse_complement()

            try:
                translation = str(seq.translate(to_stop=True))
            except Exception:
                continue

            translation = append_stop_codon(translation, start, end)

            if translation and len(translation) > 10:  # Skip tiny fragments
                header = f"{gene_id}_{start}-{end}"
                proteins.append((header, translation))

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
