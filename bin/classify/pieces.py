"""
The other pieces of the genes a contig end cut.

A draft assembly cuts the genes its contig breaks run through, and miniprot
places only the larger piece of each: it reports the placements of a protein
that score within 5% of its best. bin/contig_end_pieces.py writes out what
each such placement lacks, as a protein of its own named
PROTEIN|FIRST-LAST|LENGTH, and MINIPROT places those too. A placement of
theirs is a piece of the gene when it runs off a contig edge on the side where
the rest of the protein lies: at its 3' end for the protein's first residues,
at its 5' end for its last ones. It goes in as a placement of the protein, so
the piece is annotated under the gene's name, and check_cds.py marks it
partial at the cut end.
"""

import re

import pandas as pd

from .ids import clean_id
from .io import _read_paf

# How close to its contig's edge a piece must end on the side of the cut, as
# bin/contig_end_pieces.py asks of the placements it starts from.
EDGE = 30

# Fewer aligned residues than this place nothing, as in classify/exact.py.
MIN_RESIDUES = 8

# Share of a piece that a placement of its protein from the first pass, on the
# same strand of the same contig, must cover for the piece to be that one.
MIN_COVER = 0.5

_NAME = re.compile(r'^(.+)\|(\d+)-(\d+)\|(\d+)$')

_COLUMNS = ['Prot Id', 'Prot len', 'Strand', 'Start', 'End', 'Match Nt', 'Num nt',
            'Contig', 'Complete', 'Ref Start', 'Ref End']


def runs_off(strand, start, end, contig_len, holds_n_terminus, edge=EDGE):
    """Whether a piece ends within edge bases of its contig's edge on the side
    of the cut: the 3' side when it holds the protein's first residues, the 5'
    side when it holds its last ones. Coordinates 0-based, end excluded."""
    right, left = contig_len - end <= edge, start <= edge
    if holds_n_terminus:
        return right if strand == '+' else left
    return left if strand == '+' else right


def read_pieces(path):
    """The placements in a PAF of the pieces pass that are pieces of their
    gene, with the columns read_miniprot(..., ref_coords=True) gives: the
    protein's own ID and length, genome coordinates 1-based, and never
    complete."""
    rows = []
    for _, r in _read_paf(path).iterrows():
        name = _NAME.match(str(r['Prot Id']))
        if name is None or str(r['fasta Tag']) == '*':
            continue
        protein = name.group(1)
        first, length = int(name.group(2)), int(name.group(4))
        start, end = int(r['Start']), int(r['End'])
        if (int(r['ProtStop']) - int(r['Prot Start']) < MIN_RESIDUES
                or not runs_off(r['Strand'], start, end, int(r['fasta len']), first == 0)):
            continue
        ref = re.search(r'_(\d+)-(\d+)$', protein)
        rows.append({'Prot Id': clean_id(protein), 'Prot len': length, 'Strand': r['Strand'],
                     'Start': start + 1, 'End': end, 'Match Nt': int(r['Match Nt']),
                     'Num nt': int(r['Num nt']), 'Contig': str(r['fasta Tag']), 'Complete': False,
                     'Ref Start': float(ref.group(1)) if ref else float('nan'),
                     'Ref End': float(ref.group(2)) if ref else float('nan')})
    return (pd.DataFrame(rows, columns=_COLUMNS)
            .drop_duplicates(['Prot Id', 'Contig', 'Start', 'End']).reset_index(drop=True))


def add_contig_end_pieces(miniprot_df, path):
    """miniprot_df with the pieces in path added as placements of their
    proteins, except where the first pass placed the protein already: a cut
    that leaves two pieces of much the same size gets both placed there.

    Expects read_miniprot(..., ref_coords=True) and returns its columns plus
    'Piece', true for the placements added here. A piece is a part of its
    protein, not a copy of it, which reconciliation.recover_missing_genes
    needs to know.
    """
    pieces = read_pieces(path).assign(Piece=True)
    miniprot_df = miniprot_df.assign(Piece=False)
    placed = miniprot_df[(miniprot_df['Start'] < miniprot_df['End']) & (miniprot_df['Contig'] != '*')]

    def already_placed(p):
        o = placed[(placed['Prot Id'] == p['Prot Id']) & (placed['Contig'] == p['Contig'])
                   & (placed['Strand'] == p['Strand'])
                   & (placed['Start'] <= p['End']) & (placed['End'] >= p['Start'])]
        cover = o['End'].clip(upper=p['End']) - o['Start'].clip(lower=p['Start']) + 1
        return bool((cover >= MIN_COVER * (p['End'] - p['Start'] + 1)).any())

    keep = pd.Series([not already_placed(p) for _, p in pieces.iterrows()],
                     index=pieces.index, dtype=bool)
    new = pieces[keep]
    if new.empty:
        return miniprot_df
    return pd.concat([miniprot_df, new.reindex(columns=miniprot_df.columns)], ignore_index=True)
