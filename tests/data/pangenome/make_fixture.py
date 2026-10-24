#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cut the test fixture of PANGENE out of a real AMAP run.

    python3 tests/data/pangenome/make_fixture.py RESULTS PROTEOME OUT

RESULTS is a run's --outdir and PROTEOME the reference proteome the AMAP
container ships, /reference/aa_H37Rv.translated.fas, kept byte for byte, its
Windows line ends included.

The fixture keeps two stretches of the reference: 1.693 to 1.710 Mb, where
M. bovis lacks RD4, and 1.984 to 1.990 Mb, around the RvD2 deletion of H37Rv,
where CDC1551 and M. bovis carry three genes no reference gene matches. Every
row is the run's own; what changes is which rows there are:

- each genome keeps, stretch by stretch, the sequence from the first to the
  last of its genes that name a single reference gene of that stretch, plus
  MARGIN on either side, one FASTA record per stretch;
- its corrected table and proteins keep the genes that lie there, those
  naming no reference gene among them, with their coordinates in the whole
  genome;
- the proteome keeps the reference genes of the stretches.

The run behind tests/data/pangenome/ was AMAP at 6351c0c (main after PR #68)
under -profile docker, with a stand-in for eggNOG-mapper, which none of these
files depends on, on the first 2 Mb of CDC1551 (NC_002755.2), M. bovis
AF2122/97 (NC_002945.4) and H37Rv (AL123456.3).
"""

import csv
import os
import re
import sys

STRETCHES = {'RD4': (1693000, 1710000), 'RvD2': (1984000, 1990000)}
SAMPLES = ['L4_CDC1551', 'A4_Mbovis', 'L4_H37Rv']
MARGIN = 500        # bp kept beyond the outermost genes of a stretch
HEADER = re.compile(r'^(.+)_(\d+)-(\d+)$')


def read_fasta(path):
    """(header line, sequence lines) of each record, line ends as they are."""
    records = []
    with open(path, newline='') as fh:
        for line in fh:
            if line.startswith('>'):
                records.append([line, []])
            elif records:
                records[-1][1].append(line)
    return records


def reference_genes(proteome):
    """Each reference gene and the stretch it lies in, from the proteome's
    headers (Rv1506c_1695943-1696443)."""
    genes = {}
    for header, _ in read_fasta(proteome):
        m = HEADER.match(header[1:].strip())
        if not m:
            continue
        start, end = int(m.group(2)), int(m.group(3))
        for name, (lo, hi) in STRETCHES.items():
            if lo <= start and end <= hi:
                genes[m.group(1)] = name
    return genes


def main(results, proteome, out):
    genes = reference_genes(proteome)
    os.makedirs(os.path.join(out, 'genomes'), exist_ok=True)
    os.makedirs(os.path.join(out, 'corrected'), exist_ok=True)

    with open(os.path.join(out, 'reference.faa'), 'w', newline='') as dst:
        for header, lines in read_fasta(proteome):
            m = HEADER.match(header[1:].strip())
            if m and m.group(1) in genes:
                dst.write(header)
                dst.writelines(lines)

    for sample in SAMPLES:
        amap = os.path.join(results, 'AMAP', sample)
        with open(os.path.join(amap, f'{sample}.AMAP.tsv'), newline='') as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
            columns = list(rows[0].keys())
        # Where the stretches lie on this genome: from the first to the last of
        # its genes that name a single reference gene of each
        spans = {}
        for row in rows:
            rvs = [r.strip() for r in row['rvs'].split(',') if r.strip()]
            if len(rvs) == 1 and rvs[0] in genes:
                lo, hi = sorted((int(row['start']), int(row['end'])))
                name = genes[rvs[0]]
                old = spans.get(name, (lo, hi))
                spans[name] = (min(old[0], lo), max(old[1], hi))
        spans = {name: (max(lo - MARGIN, 1), hi + MARGIN) for name, (lo, hi) in spans.items()}

        def kept(row):
            lo, hi = sorted((int(row['start']), int(row['end'])))
            return any(a <= lo and hi <= b for a, b in spans.values())

        keep = [row for row in rows if kept(row)]
        with open(os.path.join(out, 'corrected', f'{sample}.AMAP.tsv'), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=columns, delimiter='\t', lineterminator='\n')
            w.writeheader()
            w.writerows(keep)
        ids = {row['id'] for row in keep}
        with open(os.path.join(out, 'corrected', f'{sample}.AMAP.faa'), 'w', newline='') as dst:
            for header, lines in read_fasta(os.path.join(amap, f'{sample}.AMAP.faa')):
                if header[1:].split(';')[0] in ids:
                    dst.write(header)
                    dst.writelines(lines)

        ((header, lines),) = read_fasta(os.path.join(results, 'PATHOTYPR', 'fastas', f'{sample}.fasta'))
        seq = ''.join(line.strip() for line in lines)
        with open(os.path.join(out, 'genomes', f'{sample}.fasta'), 'w') as fh:
            for name in sorted(spans, key=lambda n: spans[n]):
                lo, hi = spans[name]
                part = seq[lo - 1:hi]
                fh.write(f'>{sample}_{name} {lo}-{hi}\n')
                fh.writelines(part[i:i + 80] + '\n' for i in range(0, len(part), 80))
        print(sample, {name: f'{lo}-{hi}' for name, (lo, hi) in spans.items()},
              len(keep), 'genes')


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
