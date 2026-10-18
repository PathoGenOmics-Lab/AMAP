#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The residues of reference proteins that a contig end cut off, written out for
miniprot to look for at the other contig ends.

A draft assembly cuts the genes its contig breaks run through. miniprot, run
with --outs=0.95, reports only the placements of a protein that score within
5% of its best one, so it places the larger piece of a cut gene and not the
smaller piece at the end of another contig: in CDC1551 cut into 40 contigs,
none of the 25 smaller pieces had a placement.

A placement was cut short by its contig when its alignment leaves out residues
at one end of the protein, ends within EDGE bases of the contig's edge on that
side, and the residues it leaves out would not fit in the bases left there.
Those residues are written as a protein of their own. The piece that holds
them is then that protein's best placement, which miniprot does report, and
classify/pieces.py keeps the placements that run off a contig edge.

Records are named PROTEIN|FIRST-LAST|LENGTH: the ID of the reference protein,
the residues the record holds, 0-based and end excluded as miniprot counts
them, and the length of the whole protein.

Usage: contig_end_pieces.py SAMPLE.paf REF_PROTEOME > SAMPLE-pieces.faa
"""

import argparse
import sys

from Bio import SeqIO

__author__ = 'Paula Ruiz-Rodriguez'
__credits__ = ['Paula Ruiz-Rodriguez', 'PathoGenOmics']
__license__ = ['GPL 3']
__version__ = '1.0.1'
__maintainer__ = '@paururo'
__email__ = 'paula.ruiz.rodriguez@csic.es'
__status__ = 'finished'
__lab__ = 'PathoGenOmics (I2SysBio)'

# How close to a contig's edge an alignment must end to have been cut there:
# ten codons, the distance within which check_cds.py takes a gene's frame to
# run off its contig. classify/pieces.py holds its pieces to the same.
EDGE = 30

# The fewest residues worth looking for. As with the exact search of
# classify/exact.py, fewer say too little to be placed on sequence alone, and
# miniprot cannot seed an alignment on them anyway.
MIN_RESIDUES = 8


def cut_ends(qlen, qstart, qend, strand, tlen, tstart, tend, edge=EDGE):
    """The ends of its protein, 'N' and 'C', that a placement lacks because
    its contig ends there. Coordinates as PAF gives them: 0-based with the end
    excluded, the protein's in residues and the contig's in bases."""
    before = tstart if strand == '+' else tlen - tend   # bases on the gene's 5' side
    after = tlen - tend if strand == '+' else tstart    # and on its 3' side
    ends = []
    if qstart > 0 and before <= edge and 3 * qstart > before:
        ends.append('N')
    if qend < qlen and after <= edge and 3 * (qlen - qend) > after:
        ends.append('C')
    return ends


def missing_parts(paf):
    """(protein, first, last, length) for each stretch of a reference protein
    that one of its placements in paf lacks because a contig ends there, each
    stretch once."""
    parts = {}
    with open(paf) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            # miniprot -u writes a row with '*' for a protein it did not place
            if line.startswith('#') or len(f) < 12 or f[5] == '*':
                continue
            qlen, qstart, qend = int(f[1]), int(f[2]), int(f[3])
            for end in cut_ends(qlen, qstart, qend, f[4], int(f[6]), int(f[7]), int(f[8])):
                part = (f[0], 0, qstart, qlen) if end == 'N' else (f[0], qend, qlen, qlen)
                parts.setdefault(part, None)
    return list(parts)


def main():
    parser = argparse.ArgumentParser(
        description="Write the residues of reference proteins that a contig end cut off, "
                    "as proteins of their own, for miniprot to place.")
    parser.add_argument('paf', help="miniprot PAF of the reference proteome on the genome")
    parser.add_argument('proteome', help="the reference proteome miniprot aligned")
    args = parser.parse_args()

    parts = missing_parts(args.paf)
    wanted = {name for name, _, _, _ in parts}
    sequences = {r.id: str(r.seq) for r in SeqIO.parse(args.proteome, 'fasta') if r.id in wanted}
    written = 0
    for name, first, last, length in parts:
        seq = sequences.get(name, '')
        # miniprot counts a final '*' as a residue; it is not one to look for
        if len(seq) != length or len(seq[first:last].rstrip('*')) < MIN_RESIDUES:
            continue
        sys.stdout.write(f">{name}|{first}-{last}|{length}\n{seq[first:last]}\n")
        written += 1
    print(f"{len(parts)} stretches of reference proteins cut off by a contig end, "
          f"{written} long enough to look for", file=sys.stderr)


if __name__ == '__main__':
    main()
