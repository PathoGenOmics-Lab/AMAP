#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Deduplicate genes with the same Rv assignment.

When several genes map to the same single Rv identifier, they are ranked by:
(1) non-pseudogene over pseudogene, (2) longer CDS, (3) called by Bakta at the
same coordinates, (4) closer to where the protein sits in the reference, judged
from the genes around it (with --ref-proteome). The first keeps its place.

The rest are removed unless they are further copies of the gene: intact, at
least nearly as long as the first, and clear of every copy kept. A genome can
carry a gene more than once, and those copies are what the duplication analysis
counts; E. coli O157:H7 Sakai has several of K-12's prophage genes many times
over, and every copy but one used to be deleted here.

Also removes genes without an Rv that overlap >90% with an Rv-assigned gene.
"""

import argparse
import os
import re
import sys
from collections import defaultdict

# A further copy must be at least this long next to the first, so that pieces
# of one gene do not pass for copies of it.
MIN_COPY_LENGTH = 0.9


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


def read_reference_positions(path):
    """Where each reference protein sits in the reference, from its header
    (NAME_START-END, as classify_cases.py reads them)."""
    positions = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                m = re.match(r'^(.+)_(\d+)-(\d+)$', line[1:].split()[0])
                if m:
                    positions.setdefault(m.group(1), (int(m.group(2)), int(m.group(3))))
    return positions


def synteny_distance(genes, ref_positions):
    """How far each gene sits from where its reference protein belongs.

    Genes whose protein appears once in the table are the anchors that give
    the sample's local offset from the reference, the way classify/synteny.py
    places miniprot hits. Returns {index: distance} for the genes it can judge.
    """
    import pandas as pd
    from classify.synteny import predict_ref_positions

    rows = {}
    for i, g in enumerate(genes):
        rv = g.get('rvs', '').strip()
        if rv in ref_positions:
            rows[i] = {'Prot Id': rv, 'Start': int(g['start']), 'End': int(g['end']),
                       'Contig': g.get('contig', ''),
                       'Ref Start': ref_positions[rv][0], 'Ref End': ref_positions[rv][1]}
    if not rows:
        return {}
    df = pd.DataFrame.from_dict(rows, orient='index')
    distance = (predict_ref_positions(df) - (df['Ref Start'] + df['Ref End']) / 2).abs()
    return distance.dropna().to_dict()


def _length(g):
    return int(g['end']) - int(g['start']) + 1


def _is_pseudo(g):
    return g.get('pseudogene', 'no') == 'yes'


def _overlap(a, b):
    return (a.get('contig', '') == b.get('contig', '')
            and int(a['start']) <= int(b['end']) and int(b['start']) <= int(a['end']))


def _covered(gene, others, share=0.5):
    """Whether any of others, on gene's strand, covers share of it."""
    s, e = int(gene['start']), int(gene['end'])
    for o in others:
        if (o['strand'] == gene['strand'] and _overlap(gene, o)
                and min(e, int(o['end'])) - max(s, int(o['start'])) + 1 >= share * (e - s + 1)):
            return True
    return False


def is_copy(gene, kept, others=()):
    """Whether gene is a further copy of the gene whose copies are kept, the
    first of them its best: intact, nearly as long as it, clear of the copies
    kept, and not on the locus of another reference gene (others). In
    M. bovis an intact copy of Rv2177c lies inside Rv2424c, a pseudogene where
    it belongs, and two names would cover one locus."""
    return (not _is_pseudo(gene)
            and _length(gene) >= MIN_COPY_LENGTH * _length(kept[0])
            and not any(_overlap(gene, k) for k in kept)
            and not _covered(gene, others))


def deduplicate(genes, bakta_cds=frozenset(), distance=None):
    """Remove duplicate Rv assignments and overlapping no-Rv genes."""
    distance = distance or {}
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

    by_contig = defaultdict(list)
    for g in genes:
        by_contig[g.get('contig', '')].append(g)

    # 1. Deduplicate: same single Rv mapped to multiple genes
    for rv, indices in rv_genes.items():
        if len(indices) <= 1:
            continue
        # Rank: non-pseudogene, then longer CDS. Between copies that tie,
        # prefer one Bakta called at exactly the same coordinates, then the
        # one nearer where the protein belongs: two methods agreeing on a
        # locus, or the reference's own order of genes, beat coming first in
        # the genome, which is what decided before. In E. coli K-12 that put
        # a miniprot-only copy of b2635 at 266553 ahead of the real one.
        def score(idx):
            g = genes[idx]
            return (not _is_pseudo(g), _length(g), bakta_called(g),
                    -distance.get(idx, float('inf')))  # higher is better

        ranked = sorted(indices, key=score, reverse=True)
        kept = [genes[ranked[0]]]
        for idx in ranked[1:]:
            g = genes[idx]
            others = [o for o in by_contig[g.get('contig', '')]
                      if o.get('rvs', '').strip() not in ('', rv)]
            if is_copy(g, kept, others):
                kept.append(g)
            else:
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


def unique_ids(genes):
    """New ids for genes that share one, as {index: id}.

    Copies of a gene that only miniprot placed are all named after the
    reference protein, but the GFF3 needs one ID per feature and the master
    table joins on it, keeping one row per id. The first keeps the name and
    the others get _2, _3 and so on, the way Prokka numbers repeated names.
    """
    taken = {g['id'] for g in genes}
    seen = defaultdict(int)
    renamed = {}
    for i, g in enumerate(genes):
        seen[g['id']] += 1
        if seen[g['id']] == 1:
            continue
        n = seen[g['id']]
        while f"{g['id']}_{n}" in taken:
            n += 1
        renamed[i] = f"{g['id']}_{n}"
        taken.add(renamed[i])
    return renamed


def record_key(gene):
    """Rebuild the FASTA header check_cds.py wrote for this row.

    Headers are 'id;strand;start-end' with a 0-based start, while the table
    carries it 1-based.
    """
    return "%s;%s;%d-%s" % (gene['id'], gene['strand'],
                            int(gene['start']) - 1, gene['end'])


def filter_fasta(path, remove_keys, renames=None):
    """Rewrite a FASTA file, dropping the records named in remove_keys and
    renaming those in renames ({old header: new header}).

    Matching is on the whole header, not the part before the first ';'. A
    multi-copy gene has one record per placement, all sharing the same gene id
    (Sakai_b4570 appears six times in E. coli O157:H7), so cutting the header
    at the ';' deleted every copy as soon as one of them was removed. The gene
    stayed in the table with no sequence anywhere, and the families this hits
    are exactly the ones the duplication analysis is about.
    """
    if not os.path.exists(path):
        return
    renames = renames or {}
    lines = []
    skip = False
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                seq_id = line[1:].split()[0]
                skip = seq_id in remove_keys
                if seq_id in renames:
                    line = '>' + renames[seq_id] + line[1 + len(seq_id):]
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
    parser.add_argument('--ref-proteome', help='Reference proteome with NAME_START-END headers, '
                                               'to break the remaining ties in favour of the copy '
                                               'where the protein belongs', default=None)
    args = parser.parse_args()

    header, genes = parse_amap_tsv(args.amap_tsv)
    bakta_cds = read_bakta_cds(args.bakta) if args.bakta else frozenset()
    distance = (synteny_distance(genes, read_reference_positions(args.ref_proteome))
                if args.ref_proteome else {})
    # Keyed by placement so the other copies of a multi-copy gene survive.
    removed_keys = set()
    deduped, removed_indices = deduplicate(genes, bakta_cds, distance)
    removed_keys = {record_key(genes[i]) for i in removed_indices}
    removed_ids = sorted({genes[i]['id'] for i in removed_indices})

    renames = {}
    for i, new_id in unique_ids(deduped).items():
        old_key = record_key(deduped[i])
        deduped[i]['id'] = new_id
        renames[old_key] = record_key(deduped[i])

    if removed_indices or renames:
        with open(args.amap_tsv, 'w') as fh:
            fh.write('\t'.join(header) + '\n')
            for g in deduped:
                fh.write('\t'.join(g[h] for h in header) + '\n')
        for path in (args.faa, args.fna):
            if path:
                filter_fasta(path, removed_keys, renames)

    if removed_indices:
        print(f"Deduplicated: removed {len(removed_indices)} gene(s)", file=sys.stderr)
        for rid in removed_ids:
            print(f"  - {rid}", file=sys.stderr)
    else:
        print("No duplicate Rv assignments found", file=sys.stderr)
    if renames:
        print(f"Renamed {len(renames)} further copies to keep ids unique", file=sys.stderr)


if __name__ == '__main__':
    main()
