"""
Paralog placements: telling a protein's own locus from a copy of it elsewhere.

Miniprot reports every placement of a reference protein that scores within 5%
of its best, so a protein with a paralog is placed at both loci, and each
placement then competes for its locus with the protein that actually belongs
there. Which copy survives used to come down to insertion order and
tie-breaks. Annotating H37Rv against itself, that lost Rv1150, Rv1945 and
Rv2828c to their paralogs, and labelled every IS6110 copy with all seventeen
IS6110 genes.

Reference protein IDs carry the protein's own coordinates (Rv1150_1278269-
1278820), but those are reference coordinates: any indel upstream shifts every
position after it, by 6.4 kb by the end of the CP110619 assembly of H37Rv. So
each placement's reference position is predicted from the nearest proteins
placed only once, which follow the sample's local offset from the reference,
and the placement closest to where its protein belongs is the syntenic one.

A placement that is not syntenic is dropped only where another protein's
syntenic placement covers it, i.e. where the locus is already accounted for,
and matches there at least as well. A copy at a locus no protein claims, such
as an IS6110 insertion at a site the reference does not have, or a duplication
in the sample, is left alone; so is a protein that fits a locus better than its
owner does, which happens in tandem repeats whose copy number differs: E. coli
Sakai carries two of K-12's three ldr toxin copies, and position alone would
label the b4421-type copy as b4423, two residues off.
"""

import numpy as np
import pandas as pd
from intervaltree import IntervalTree

# How far from its predicted position a placement may sit and still count as
# its protein's own locus. Absorbs the error of an offset estimated next to an
# indel, and is small next to the distance between paralogs in a genome.
SYNTENY_TOLERANCE = 10_000

# Single-placement proteins on each side used to estimate the local offset.
N_ANCHORS = 5

# Share of a placement another protein's syntenic placement must cover to claim
# its locus. Genes that merely overlap at their ends stay well below it.
MIN_COVER = 0.5

_HELPER_COLS = ('Ref Start', 'Ref End')


def _midpoint(start, end):
    return (start + end) / 2


def predict_ref_positions(df):
    """Predict the reference position of every placement's locus.

    Uses proteins with a single placement as anchors, per sequence, so a
    multi-contig assembly never mixes offsets between contigs. A sequence that
    runs against the reference (reverse-complemented contig) maps positions
    backwards, so its orientation is read off its anchors first. Placements on
    a sequence with too few anchors get NaN and are left alone downstream.
    """
    pred = pd.Series(np.nan, index=df.index, dtype=float)
    n_place = df.groupby('Prot Id')['Prot Id'].transform('size')
    anchors = df[(n_place == 1) & df['Ref Start'].notna()]

    for contig, a in anchors.groupby('Contig'):
        if len(a) < 3:
            continue
        a = a.sort_values('Start')
        x = _midpoint(a['Start'], a['End']).to_numpy(dtype=float)
        r = _midpoint(a['Ref Start'], a['Ref End']).to_numpy(dtype=float)
        sign = -1.0 if np.median(np.sign(np.diff(r))) < 0 else 1.0

        rows = df[df['Contig'] == contig]
        mids = _midpoint(rows['Start'], rows['End']).to_numpy(dtype=float)
        for idx, m in zip(rows.index, mids):
            i = np.searchsorted(x, m)
            lo, hi = max(0, i - N_ANCHORS), min(len(x), i + N_ANCHORS)
            pred[idx] = np.median(r[lo:hi] + sign * (m - x[lo:hi]))
    return pred


def mark_syntenic(df):
    """Flag, per protein, the placement that sits where the protein belongs.

    The closest placement to the protein's predicted reference position is its
    syntenic one, provided it is within SYNTENY_TOLERANCE. Every other
    placement of the protein, and every placement that could not be predicted,
    is not syntenic.
    """
    pred = predict_ref_positions(df)
    dist = (pred - _midpoint(df['Ref Start'], df['Ref End'])).abs()

    syntenic = pd.Series(False, index=df.index)
    known = dist.notna()
    if known.any():
        best = dist[known].groupby(df.loc[known, 'Prot Id']).idxmin()
        syntenic[best.values] = (dist[best.values] <= SYNTENY_TOLERANCE).values
    return syntenic, known


def drop_paralog_placements(df):
    """Drop placements made at a locus that another protein syntenically owns.

    Expects read_miniprot(..., ref_coords=True) and returns its usual columns.
    When no protein carries reference coordinates there is nothing to go on,
    and the frame comes back unchanged.
    """
    out_cols = [c for c in df.columns if c not in _HELPER_COLS]
    if 'Ref Start' not in df.columns or df['Ref Start'].isna().all():
        return df[out_cols]

    # miniprot -u writes a row for every protein it could not place, with '*'
    # for the target and zeros for coordinates: Rv3427c in CDC1551, 31 proteins
    # in M. bovis. They are neither anchors nor placements to judge, and go
    # through untouched as they always have.
    placed = df[(df['Start'] < df['End']) & (df['Contig'] != '*') & df['Ref Start'].notna()]
    syntenic, known = mark_syntenic(placed)

    identity = (placed['Match Nt'] / placed['Num nt'].where(placed['Num nt'] > 0)).fillna(0)
    # A locus is a reading frame on a strand, not just a stretch of sequence:
    # a paralog aligns over the same codons as the gene it lands on, while
    # genes that really overlap (Rv1041c and Rv1042c) sit in different frames
    # and must never cost each other a place. Read off the 5' end.
    frame = placed['Start'].where(placed['Strand'] != '-', placed['End']) % 3

    def covering(tree, idx, beats):
        """An entry of tree, for another protein in the same frame, that
        covers MIN_COVER of placement idx and satisfies beats(its identity).

        Two proteins whose reference spans overlap are neighbours, not
        competitors: Rv2879c and Rv2880c overlap in H37Rv, and an indel in
        CDC1551 and M. bovis brings them into one frame. Only a protein that
        lives elsewhere in the reference can be squatting on a locus.
        """
        row = placed.loc[idx]
        length = row['End'] - row['Start'] + 1
        for iv in tree[row['Start']:row['End'] + 1]:
            prot, ident, fr, ref_s, ref_e = iv.data
            if prot == row['Prot Id'] or fr != frame[idx] or not beats(ident):
                continue
            if ref_s <= row['Ref End'] and row['Ref Start'] <= ref_e:
                continue
            if min(row['End'] + 1, iv.end) - max(row['Start'], iv.begin) >= MIN_COVER * length:
                return True
        return False

    # Two proteins can both sit where they belong and still land on one
    # locus: after a deletion, a lost gene's nearest place is often its paralog
    # next door. CDC1551 lacks Rv3427c-Rv3428c, and Rv3425 and Rv3426 then sit
    # on Rv3429 as close to home as anything left. The better match keeps the
    # claim; the others give it up, and are dropped below like any paralog copy.
    owners = {}
    for idx in identity[syntenic].sort_values(ascending=False).index:
        row = placed.loc[idx]
        tree = owners.setdefault((row['Contig'], row['Strand']), IntervalTree())
        if covering(tree, idx, lambda ident: ident > identity[idx]):
            syntenic[idx] = False
            continue
        tree[row['Start']:row['End'] + 1] = (row['Prot Id'], identity[idx], frame[idx],
                                             row['Ref Start'], row['Ref End'])

    drop = []
    for idx in placed.index[known & ~syntenic]:
        tree = owners.get((placed.at[idx, 'Contig'], placed.at[idx, 'Strand']))
        if tree is not None and covering(tree, idx, lambda ident: ident >= identity[idx]):
            drop.append(idx)

    return df.drop(index=drop)[out_cols].reset_index(drop=True)
