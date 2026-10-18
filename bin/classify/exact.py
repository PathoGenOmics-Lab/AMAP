"""
Exact placements for the proteins miniprot leaves unplaced, or places only
with mismatches.

miniprot needs a seed of several residues and enough score to report an
alignment, which a protein of a few dozen residues rarely gives it: annotating
E. coli K-12 with its own proteome, it placed none of 33 small proteins of 8 to
29 residues, a third of the genes the annotation then got wrong. A protein the
genome encodes exactly needs no alignment. It is looked for in the six-frame
translation, preceded by a start codon and followed by a stop, which pins the
whole gene down.
"""

import re

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from .ids import clean_id

START_CODONS = {'ATG', 'GTG', 'TTG'}

# Shorter proteins say too little to be placed on sequence identity alone.
MIN_RESIDUES = 8


def _six_frames(genome_path):
    """(contig, strand, frame offset, translation, length) for every frame."""
    frames = []
    for rec in SeqIO.parse(genome_path, 'fasta'):
        seq = str(rec.seq).upper()
        for strand, s in (('+', seq), ('-', str(Seq(seq).reverse_complement()))):
            for f in range(3):
                end = f + (len(s) - f) // 3 * 3
                frames.append((rec.id, strand, f, s, str(Seq(s[f:end]).translate()), len(seq)))
    return frames


def find_exact(protein, frames):
    """Every place protein is encoded exactly: a start codon, its other
    residues, then a stop. Yields (contig, start, end, strand), 1-based on the
    forward strand, stop codon included."""
    body = protein.rstrip('*')
    rest = body[1:]
    for contig, strand, f, s, aa, length in frames:
        i = aa.find(rest)
        while i != -1:
            first, stop = i - 1, i + len(rest)
            p = f + 3 * first
            if first >= 0 and stop < len(aa) and aa[stop] == '*' and s[p:p + 3] in START_CODONS:
                q = p + 3 * (len(body) + 1) - 1
                yield (contig, p + 1, q + 1, '+') if strand == '+' else (contig, length - q, length - p, '-')
            i = aa.find(rest, i + 1)


def place_exact_matches(miniprot_df, genome_path, proteome_path, min_residues=MIN_RESIDUES):
    """Add a placement for every exact copy of a protein miniprot did not
    place exactly.

    That is a protein miniprot left unplaced, and also one it only placed with
    mismatches: the ldr toxins of E. coli are 19 residues that differ in three
    or four, and miniprot put b4665 on b4664's locus and nowhere else, so once
    that placement went to its owner, b4665 had none, though its own locus
    carries it exactly. A copy found where miniprot already placed the protein
    is left to miniprot's alignment.

    Expects read_miniprot(..., ref_coords=True). The placements go in as
    complete alignments, with the reference coordinates from the protein's
    header, and replace the protein's unplaced row.
    """
    placed = miniprot_df[(miniprot_df['Start'] < miniprot_df['End'])
                         & (miniprot_df['Contig'] != '*')]
    exact = placed['Match Nt'] == placed['Num nt']
    if 'Complete' in placed.columns:
        exact &= placed['Complete'].astype(bool)
    have = set(placed.loc[exact, 'Prot Id'])
    wanted = {}
    for rec in SeqIO.parse(proteome_path, 'fasta'):
        name = clean_id(rec.id)
        seq = str(rec.seq).upper()
        if name not in have and name not in wanted and len(seq.rstrip('*')) >= min_residues:
            wanted[name] = (rec.id, seq)
    if not wanted:
        return miniprot_df

    def already_placed(name, contig, strand, start, end):
        p = placed[(placed['Prot Id'] == name) & (placed['Contig'] == contig)
                   & (placed['Strand'] == strand) & (placed['Start'] <= end) & (placed['End'] >= start)]
        cover = (p['End'].clip(upper=end) - p['Start'].clip(lower=start) + 1)
        return bool((cover >= 0.5 * (end - start + 1)).any())

    frames = _six_frames(genome_path)
    rows = []
    for name, (full_id, seq) in wanted.items():
        ref = re.search(r'_(\d+)-(\d+)$', full_id)
        for contig, start, end, strand in find_exact(seq, frames):
            if already_placed(name, contig, strand, start, end):
                continue
            row = {'Prot Id': name, 'Prot len': len(seq), 'Strand': strand,
                   'Start': start, 'End': end, 'Match Nt': end - start + 1,
                   'Num nt': end - start + 1, 'Contig': contig, 'Complete': True}
            if 'Ref Start' in miniprot_df.columns:
                row['Ref Start'] = float(ref.group(1)) if ref else float('nan')
                row['Ref End'] = float(ref.group(2)) if ref else float('nan')
            rows.append(row)
    if not rows:
        return miniprot_df

    found = {r['Prot Id'] for r in rows}
    unplaced = miniprot_df['Prot Id'].isin(found) & ~miniprot_df.index.isin(placed.index)
    added = pd.DataFrame(rows)[[c for c in miniprot_df.columns]]
    return pd.concat([miniprot_df[~unplaced], added], ignore_index=True)
