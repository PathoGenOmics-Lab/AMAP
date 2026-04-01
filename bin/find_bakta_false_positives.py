#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find Bakta false positives by comparing Bakta annotation against a trusted
curated annotation of the SAME genome.

A Bakta CDS is a "false positive" if NONE of its coordinates overlap with
ANY CDS in the curated annotation.

Input:
  --bakta-tsv      Bakta TSV output (from annotating the reference genome)
  --curated-faa    Curated proteome FAA (from extract_ref_proteins.py)
  --bakta-faa      Bakta protein FASTA (from annotating the reference genome)
  --prefix         Output prefix

Output:
  {prefix}_bad_bakta.faa   — protein sequences of Bakta false positives
"""

import argparse
import sys
from pathlib import Path
from intervaltree import IntervalTree
from Bio import SeqIO
import pandas as pd


def load_curated_coords(curated_faa):
    """
    Extract coordinates from curated protein headers.
    Expected format: >GENE_ID_START-END (e.g., >Rv0001_1-1524)
    """
    import re
    tree = IntervalTree()
    count = 0
    for record in SeqIO.parse(curated_faa, "fasta"):
        m = re.search(r'_(\d+)-(\d+)$', record.id)
        if m:
            start = int(m.group(1))
            end = int(m.group(2))
            if start < end:
                tree[start:end] = record.id
                count += 1
    print(f"Loaded {count} curated gene coordinates")
    return tree


def load_bakta_cds(bakta_tsv):
    """Read Bakta TSV and return CDS entries with coordinates."""
    df = pd.read_csv(bakta_tsv, comment='#', sep='\t', header=None)
    df.columns = [
        "Sequence Id", "Type", "Start", "End", "Strand",
        "Locus Tag", "Gene", "Product", "DbXrefs"
    ]
    cds = df[df['Type'] == 'cds'].copy()
    return cds


def find_false_positives(bakta_cds, curated_tree, min_overlap_frac=0.1):
    """
    Find Bakta CDS that don't overlap with ANY curated gene.
    
    A Bakta CDS is a false positive if there is no curated gene
    overlapping its coordinates by at least min_overlap_frac of the
    Bakta gene length (default 10%). This prevents large spurious
    Bakta genes that barely touch a real gene from being kept.
    """
    false_positives = []
    for _, row in bakta_cds.iterrows():
        start = int(row['Start'])
        end = int(row['End'])
        locus_tag = row['Locus Tag']
        gene_len = end - start

        # Query the interval tree for overlapping curated genes
        overlaps = curated_tree[start:end]

        if len(overlaps) == 0:
            # No overlap at all → definitely a false positive
            false_positives.append(locus_tag)
        else:
            # Check if any overlap is substantial relative to Bakta gene length
            has_substantial = False
            min_bp = max(30, gene_len * min_overlap_frac)
            for iv in overlaps:
                overlap_bp = min(end, iv.end) - max(start, iv.begin)
                if overlap_bp >= min_bp:
                    has_substantial = True
                    break
            if not has_substantial:
                false_positives.append(locus_tag)

    return set(false_positives)


def extract_false_positive_proteins(bakta_faa, false_positive_ids, output_faa):
    """Extract protein sequences of false positives from Bakta FAA."""
    count = 0
    with open(output_faa, 'w') as out:
        for record in SeqIO.parse(bakta_faa, "fasta"):
            # Bakta FAA headers can be just the locus tag or more complex
            record_id = record.id.split()[0]
            if record_id in false_positive_ids:
                SeqIO.write(record, out, "fasta")
                count += 1
    return count


def main():
    parser = argparse.ArgumentParser(
        description="Find Bakta false positives by comparison to curated annotation.")
    parser.add_argument("--bakta-tsv", required=True,
                        help="Bakta TSV from annotating the reference genome.")
    parser.add_argument("--curated-faa", required=True,
                        help="Curated proteome (from extract_ref_proteins.py).")
    parser.add_argument("--bakta-faa", required=True,
                        help="Bakta protein FASTA from annotating the reference genome.")
    parser.add_argument("--prefix", default="custom_ref",
                        help="Output prefix.")
    args = parser.parse_args()

    # Load curated coordinates
    curated_tree = load_curated_coords(args.curated_faa)
    if len(curated_tree) == 0:
        sys.exit("Error: No coordinates found in curated FAA. Check header format (GENE_START-END).")

    # Load Bakta CDS
    bakta_cds = load_bakta_cds(args.bakta_tsv)
    print(f"Bakta annotated {len(bakta_cds)} CDS on the reference genome")

    # Find false positives
    fp_ids = find_false_positives(bakta_cds, curated_tree)
    print(f"Found {len(fp_ids)} Bakta false positives (no match to curated annotation)")

    # Extract their protein sequences
    output_faa = f"{args.prefix}_bad_bakta.faa"
    n = extract_false_positive_proteins(args.bakta_faa, fp_ids, output_faa)
    print(f"Wrote {n} false positive proteins → {output_faa}")

    if n == 0:
        # Write an empty file so downstream doesn't crash
        Path(output_faa).touch()
        print("Warning: No false positives found. Using empty bad_bakta file.")


if __name__ == "__main__":
    main()
