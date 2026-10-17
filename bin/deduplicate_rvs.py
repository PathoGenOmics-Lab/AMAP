#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Deduplicate genes with the same Rv assignment.

When multiple genes map to the same single Rv identifier, keeps only the best
candidate based on: (1) non-pseudogene over pseudogene, (2) longer CDS.

Also removes genes without an Rv that overlap >90% with an Rv-assigned gene.
"""

import argparse
import sys
import os
from collections import defaultdict


def parse_amap_tsv(path):
    """Parse AMAP TSV (with header) into list of dicts."""
    genes = []
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            genes.append(dict(zip(header, fields)))
    return header, genes


def overlap_fraction(s1, e1, s2, e2):
    """Return fraction of gene1 overlapping gene2."""
    ov = max(0, min(e1, e2) - max(s1, s2) + 1)
    length = e1 - s1 + 1
    return ov / length if length > 0 else 0


def read_bakta_cds(path):
    """(start, end, strand, contig) of every CDS in a Bakta TSV."""
    cds = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) > 4 and f[1] == 'cds':
                cds.add((int(f[2]), int(f[3]), f[4], f[0]))
    return cds


def deduplicate(genes, bakta_cds=frozenset()):
    """Remove duplicate Rv assignments and overlapping no-Rv genes."""
    # Index single-Rv genes by Rv
    rv_genes = defaultdict(list)
    for i, g in enumerate(genes):
        rvs = g.get('rvs', '').strip()
        if rvs and ',' not in rvs:
            rv_genes[rvs].append(i)

    remove = set()
    # A table from before contigs were carried through has no contig column;
    # then a Bakta call on any contig at those coordinates counts.
    bakta_any_contig = {c[:3] for c in bakta_cds}

    def bakta_called(g):
        key = (int(g['start']), int(g['end']), g['strand'])
        contig = g.get('contig', '')
        return key + (contig,) in bakta_cds if contig else key in bakta_any_contig

    # 1. Deduplicate: same single Rv mapped to multiple genes
    for rv, indices in rv_genes.items():
        if len(indices) <= 1:
            continue
        # Score each: prefer non-pseudogene, then longer CDS. Between copies
        # that tie, prefer one Bakta called at exactly the same coordinates:
        # two methods agreeing on a locus beats coming first in the genome,
        # which is what decided before. In E. coli K-12 that kept a
        # miniprot-only copy of b2635 at 266553 and dropped the real one.
        def score(idx):
            g = genes[idx]
            is_pseudo = 1 if g.get('pseudogene', 'no') == 'yes' else 0
            length = int(g['end']) - int(g['start']) + 1
            return (-is_pseudo, length, bakta_called(g))  # higher is better

        ranked = sorted(indices, key=score, reverse=True)
        # Keep best, remove rest
        for idx in ranked[1:]:
            remove.add(idx)

    # 2. Remove no-Rv genes that overlap >90% with an Rv-assigned gene
    rv_assigned = []
    no_rv = []
    for i, g in enumerate(genes):
        if i in remove:
            continue
        rvs = g.get('rvs', '').strip()
        s, e = int(g['start']), int(g['end'])
        if rvs:
            rv_assigned.append((i, s, e, g.get('contig', '')))
        else:
            no_rv.append((i, s, e, g.get('contig', '')))

    # Only genes on the same contig can overlap: every contig starts at 1.
    for ni, ns, ne, nc in no_rv:
        for _, rs, re_, rc in rv_assigned:
            if nc == rc and overlap_fraction(ns, ne, rs, re_) > 0.9:
                remove.add(ni)
                break

    return [g for i, g in enumerate(genes) if i not in remove], remove


def record_key(gene):
    """Rebuild the FASTA header check_cds.py wrote for this row.

    Headers are 'id;strand;start-end' with a 0-based start, while the table
    carries it 1-based.
    """
    return "%s;%s;%d-%s" % (gene['id'], gene['strand'],
                            int(gene['start']) - 1, gene['end'])


def filter_fasta(path, remove_keys):
    """Rewrite a FASTA file, dropping the records named in remove_keys.

    Matching is on the whole header, not the part before the first ';'. A
    multi-copy gene has one record per placement, all sharing the same gene id
    (Sakai_b4570 appears six times in E. coli O157:H7), so cutting the header
    at the ';' deleted every copy as soon as one of them was removed. The gene
    stayed in the table with no sequence anywhere, and the families this hits
    are exactly the ones the duplication analysis is about.
    """
    if not os.path.exists(path):
        return
    lines = []
    skip = False
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                seq_id = line[1:].split()[0]
                skip = seq_id in remove_keys
            if not skip:
                lines.append(line)
    with open(path, 'w') as fh:
        fh.writelines(lines)


def main():
    parser = argparse.ArgumentParser(description='Deduplicate Rv assignments in AMAP output')
    parser.add_argument('amap_tsv', help='AMAP TSV file (with header)')
    parser.add_argument('--faa', help='AMAP FAA file to filter', default=None)
    parser.add_argument('--fna', help='AMAP FNA file to filter', default=None)
    parser.add_argument('--bakta', help='Bakta TSV, to break ties between copies of a gene '
                                        'in favour of one Bakta also called', default=None)
    args = parser.parse_args()

    header, genes = parse_amap_tsv(args.amap_tsv)
    n_before = len(genes)
    bakta_cds = read_bakta_cds(args.bakta) if args.bakta else frozenset()
    deduped, removed_indices = deduplicate(genes, bakta_cds)
    n_removed = n_before - len(deduped)

    if n_removed > 0:
        # Collect the removed genes for FASTA filtering, keyed by placement so
        # the other copies of a multi-copy gene survive.
        all_genes_orig = parse_amap_tsv(args.amap_tsv)[1]
        removed_keys = {record_key(all_genes_orig[i]) for i in removed_indices}
        removed_ids = {all_genes_orig[i]['id'] for i in removed_indices}

        # Rewrite TSV
        with open(args.amap_tsv, 'w') as fh:
            fh.write('\t'.join(header) + '\n')
            for g in deduped:
                fh.write('\t'.join(g[h] for h in header) + '\n')

        # Filter FASTAs
        if args.faa:
            filter_fasta(args.faa, removed_keys)
        if args.fna:
            filter_fasta(args.fna, removed_keys)

        print(f"Deduplicated: removed {n_removed} gene(s)", file=sys.stderr)
        for rid in sorted(removed_ids):
            print(f"  - {rid}", file=sys.stderr)
    else:
        print("No duplicate Rv assignments found", file=sys.stderr)


if __name__ == '__main__':
    main()
