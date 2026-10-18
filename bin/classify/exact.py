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

A small protein the genome does not encode exactly is looked for again,
allowing a few substitutions, since miniprot cannot place it either: E. coli
O157:H7 Sakai carries ten of K-12's micro-proteins with one to four residues
changed, and none of them was annotated.
"""

import re

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from .ids import clean_id

START_CODONS = {'ATG', 'GTG', 'TTG'}

# Shorter proteins say too little to be placed on sequence identity alone.
MIN_RESIDUES = 8

# A copy with substitutions is looked for only for proteins this short, which
# miniprot cannot seed, and may differ in at most this share of the residues
# after the first. A stop codon in place of a residue counts as a difference:
# the gene is then a pseudogene, as check_cds.py will say.
NEAR_MAX_RESIDUES = 50
NEAR_MAX_DIFFERENCES = 0.2


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


def find_near(protein, frames, max_differences=NEAR_MAX_DIFFERENCES):
    """Every place protein is encoded with a few substitutions: a start codon,
    one codon for each of its other residues, at most max_differences of them
    different, then a stop. Yields (contig, start, end, strand, differences)
    as find_exact does.

    With at most k differences, one of k + 1 pieces of the protein is encoded
    exactly, so only the places where a piece is found are compared."""
    body = protein.rstrip('*')
    rest = body[1:]
    n = len(rest)
    k = int(max_differences * n)
    cuts = [round(j * n / (k + 1)) for j in range(k + 2)]
    for contig, strand, f, s, aa, length in frames:
        tried = set()
        for a, b in zip(cuts, cuts[1:]):
            i = aa.find(rest[a:b])
            while i != -1:
                at = i - a
                if at not in tried:
                    tried.add(at)
                    first, stop = at - 1, at + n
                    p = f + 3 * first
                    if first >= 0 and stop < len(aa) and aa[stop] == '*' and s[p:p + 3] in START_CODONS:
                        differences = sum(x != y for x, y in zip(aa[at:stop], rest))
                        if differences <= k:
                            q = p + 3 * (len(body) + 1) - 1
                            yield ((contig, p + 1, q + 1, '+', differences) if strand == '+'
                                   else (contig, length - q, length - p, '-', differences))
                i = aa.find(rest[a:b], i + 1)


def place_exact_matches(miniprot_df, genome_path, proteome_path, min_residues=MIN_RESIDUES):
    """Add a placement for every exact copy of a protein miniprot did not
    place exactly.

    That is a protein miniprot left unplaced, and also one it only placed with
    mismatches: the ldr toxins of E. coli are 19 residues that differ in three
    or four, and miniprot put b4665 on b4664's locus and nowhere else, so once
    that placement went to its owner, b4665 had none, though its own locus
    carries it exactly. A copy found where miniprot already placed the protein
    is left to miniprot's alignment.

    A protein of up to NEAR_MAX_RESIDUES with no exact copy is looked for with
    a few substitutions (find_near); each one counts as three bases that do
    not match.

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
        copies = [hit + (0,) for hit in find_exact(seq, frames)]
        if not copies and len(seq.rstrip('*')) <= NEAR_MAX_RESIDUES:
            copies = list(find_near(seq, frames))
        for contig, start, end, strand, differences in copies:
            if already_placed(name, contig, strand, start, end):
                continue
            row = {'Prot Id': name, 'Prot len': len(seq), 'Strand': strand,
                   'Start': start, 'End': end, 'Match Nt': end - start + 1 - 3 * differences,
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
