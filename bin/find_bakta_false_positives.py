#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find Bakta false positives by comparing Bakta annotation against a trusted
curated annotation of the SAME genome.

A Bakta CDS is a "false positive" if NONE of its coordinates overlap with
ANY CDS the curated annotation has on the same sequence.

Input:
  --bakta-tsv      Bakta TSV output (from annotating the reference genome,
                   with --keep-contig-headers so that it names its sequences
                   as the FASTA does)
  --curated-faa    Curated proteome FAA (from extract_ref_proteins.py)
  --bakta-faa      Bakta protein FASTA (from annotating the reference genome)
  --ref-genome     The reference genome FASTA Bakta annotated
  --prefix         Output prefix

Output:
  {prefix}_bad_bakta.faa   — protein sequences of Bakta false positives
"""

import argparse
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from intervaltree import IntervalTree
from Bio import SeqIO
import pandas as pd

# Allow importing from bin/classify/ when run as a script, as classify_cases.py does.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classify import circular

# The sequence extract_ref_proteins.py names after each protein's coordinates
_CONTIG = re.compile(r'(?:^|\s)contig=(\S+)')


def load_curated_coords(curated_faa):
    """
    Extract coordinates from curated protein headers, by sequence.
    Expected format: >GENE_ID_START-END contig=SEQUENCE
    (e.g., >Rv0001_1-1524 contig=NC_000962.3)

    Returns {sequence: [(start, end, protein)]}. A gene across the origin of
    a circular sequence ends past its last base (see extract_ref_proteins.py).
    Genes whose header names no sequence go under None, compared with the
    calls on every sequence, as all genes were before.
    """
    curated = defaultdict(list)
    count = 0
    for record in SeqIO.parse(curated_faa, "fasta"):
        m = re.search(r'_(\d+)-(\d+)$', record.id)
        if m:
            start = int(m.group(1))
            end = int(m.group(2))
            if start < end:
                contig = _CONTIG.search(record.description)
                curated[contig.group(1) if contig else None].append((start, end, record.id))
                count += 1
    print(f"Loaded {count} curated gene coordinates")
    return dict(curated)


def load_bakta_cds(bakta_tsv):
    """Read Bakta TSV and return CDS entries with coordinates."""
    # Sequence names as text: a FASTA whose sequences are called 1, 2, 3 would
    # otherwise give numbers that no curated header's name equals
    df = pd.read_csv(bakta_tsv, comment='#', sep='\t', header=None, dtype={0: str})
    df.columns = [
        "Sequence Id", "Type", "Start", "End", "Strand",
        "Locus Tag", "Gene", "Product", "DbXrefs"
    ]
    cds = df[df['Type'] == 'cds'].copy()
    return cds


def _as_bakta_writes(start, end, contig):
    """A curated gene the way Bakta and classify/circular.py write one across
    the origin: start after end."""
    try:
        n = circular.contig_length(contig)
    except ValueError:
        return start, end
    return (start, end - n) if end > n else (start, end)


def _trees(curated):
    """An IntervalTree of the curated genes of each sequence. A gene across
    the origin goes in as both of its copies on the unrolled circle (see
    classify/circular.py), so that a call on either side of the origin finds it."""
    trees = {}
    for contig, genes in curated.items():
        tree = trees[contig] = IntervalTree()
        for start, end, _ in genes:
            gene = _as_bakta_writes(start, end, contig)
            for s, e in circular.copies(*gene, contig):
                tree[s:e + 1] = gene
    return trees


def find_false_positives(bakta_cds, curated, min_overlap_frac=0.1):
    """
    Find Bakta CDS that don't overlap with ANY curated gene of their sequence.

    A Bakta CDS is a false positive if there is no curated gene
    overlapping its coordinates by at least min_overlap_frac of the
    Bakta gene length (default 10%). This prevents large spurious
    Bakta genes that barely touch a real gene from being kept.

    All the curated genes used to be one tree of coordinates, whatever their
    sequence: a call on a plasmid counted as matched by a chromosome gene at
    the same numbers. And a call across the origin of a circular sequence,
    which Bakta writes start after end (stcE on plasmid pO157 of E. coli
    O157:H7 Sakai is 92563..2502), matched nothing at all, so its protein was
    listed as a false positive.

    The calls on a sequence the curated annotation has no gene on are not
    judged: it says nothing about them. Returns the false positives and the
    number of calls left unjudged on each such sequence.
    """
    trees = _trees(curated)
    everywhere = trees.get(None)
    false_positives, unjudged = [], Counter()
    for _, row in bakta_cds.iterrows():
        contig = row['Sequence Id']
        call = (int(row['Start']), int(row['End']))
        locus_tag = row['Locus Tag']
        own = trees.get(contig)
        if own is None and everywhere is None:
            unjudged[contig] += 1
            continue
        gene_len = circular.length(*call, contig)

        # Query the trees for overlapping curated genes
        overlaps = set()
        for tree in (own, everywhere):
            if tree is not None:
                for s, e in circular.copies(*call, contig):
                    overlaps |= {iv.data for iv in tree[s:e + 1]}

        # Check if any overlap is substantial relative to Bakta gene length
        min_bp = max(30, gene_len * min_overlap_frac)
        if not any(circular.shared(call, gene, contig) >= min_bp for gene in overlaps):
            false_positives.append(locus_tag)

    return set(false_positives), unjudged


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
    parser.add_argument("--ref-genome", required=True,
                        help="The reference genome FASTA Bakta annotated, for the length of each "
                             "sequence: genes across the origin of a circular one take it.")
    parser.add_argument("--prefix", default="custom_ref",
                        help="Output prefix.")
    args = parser.parse_args()

    circular.set_lengths(circular.read_lengths(args.ref_genome))

    # Load curated coordinates
    curated = load_curated_coords(args.curated_faa)
    if not curated:
        sys.exit("Error: No coordinates found in curated FAA. Check header format (GENE_START-END).")

    # Load Bakta CDS
    bakta_cds = load_bakta_cds(args.bakta_tsv)
    print(f"Bakta annotated {len(bakta_cds)} CDS on the reference genome")

    # The curated genes and Bakta's calls must name the sequences alike, or
    # nothing could be compared
    named = set(curated) - {None}
    if None not in curated and len(bakta_cds) and not named & set(bakta_cds['Sequence Id']):
        sys.exit(f"Error: Bakta's table names sequences ({', '.join(sorted(set(bakta_cds['Sequence Id']))[:5])}) "
                 f"that the curated proteome does not ({', '.join(sorted(named)[:5])}). "
                 f"Run Bakta with --keep-contig-headers on the genome the proteome was extracted with.")

    # Find false positives
    fp_ids, unjudged = find_false_positives(bakta_cds, curated)
    print(f"Found {len(fp_ids)} Bakta false positives (no match to curated annotation)")
    for contig, n in sorted(Counter(bakta_cds.loc[bakta_cds['Locus Tag'].isin(fp_ids), 'Sequence Id']).items()):
        print(f"    {contig}: {n}")
    for contig, n in sorted(unjudged.items()):
        print(f"Note: the curated annotation has no gene on {contig}; its {n} Bakta CDS are not judged")

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
