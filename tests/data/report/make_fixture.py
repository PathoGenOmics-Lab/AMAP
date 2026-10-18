#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cut the test fixture of amap_report.py out of a real AMAP run.

    python3 tests/data/report/make_fixture.py RESULTS PROTEOME MYCOBROWSER OUT

RESULTS is a run's --outdir; PROTEOME and MYCOBROWSER are the reference files
the AMAP container ships, /reference/aa_H37Rv.translated.fas and
/reference/Mycobacterium_tuberculosis_H37Rv_txt_v4.txt.

The fixture keeps four genomes and two stretches of the reference: 1.690 to
1.712 Mb, where M. bovis lacks RD4, and 3.780 to 3.940 Mb, where these genomes
have a gene in every state the report draws. Every row is the run's own; what
changes is which rows there are:

- each genome keeps the rows of its tables that lie where the genes naming a
  reference gene of those stretches lie, contig by contig;
- the proteome keeps those reference genes only, so that every gene of it has
  all the genes naming it in the fixture, and its state there is its state in
  the run (an IS6110 copy still names the copies outside the stretches);
- the case summaries count the rows kept, so the counts and the tables agree;
- Bakta's tables keep the CDS rows of the genes in the final table, without
  DbXrefs, Mycobrowser's the columns the report reads, and Bakta's GFF3 its
  comment lines, which name every contig and its length.

The run behind tests/data/report/ was AMAP at fa48972 (main after PR #46)
under -profile docker --skip_functional_annotation, on genomes from NCBI:
CP110619 (H37Rv), CDC1551 cut into 40 contigs, M. bovis AF2122/97 and
CCDC5079 (GCA_000270345.1).
"""

import csv
import os
import re
import shutil
import sys
from collections import defaultdict

WINDOWS = [(1690000, 1712000), (3780000, 3940000)]
SAMPLES = ['L4_H37Rv', 'L4_CDC1551draft', 'A4_AF2122', 'L2_CCDC5079']
MARGIN = 1500       # bp kept beyond the outermost genes of a stretch
GAP = 20000         # genes further apart than this start a new stretch
MYCO_COLUMNS = ['Locus', 'Name', 'Feature', 'Product', 'Functional_Category',
                'Is_Pseudogene', 'Drug Resistance Mutations']
HEADER = re.compile(r'^(.+)_(\d+)-(\d+)$')


def read_tsv(path):
    with open(path, newline='') as fh:
        rows = list(csv.reader(fh, delimiter='\t'))
    return rows[0], rows[1:]


def write_tsv(path, header, rows):
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(header)
        w.writerows(rows)


def overlaps(intervals, contig, start, end):
    return any(c == contig and start <= b and end >= a for c, a, b in intervals)


def main(results, proteome, mycobrowser, out):
    os.makedirs(out, exist_ok=True)
    proteins, order = {}, []
    name = None
    with open(proteome) as fh:
        for line in fh:
            if line.startswith('>'):
                name = line[1:].strip()
                order.append(name)
                proteins[name] = []
            else:
                proteins[name].append(line)
    coords = {}
    for header in order:
        m = HEADER.match(header)
        coords[m.group(1)] = (int(m.group(2)), int(m.group(3)), header)
    window = {g for g, (s, e, _) in coords.items() if any(s <= b and e >= a for a, b in WINDOWS)}

    for sample in SAMPLES:
        amap = os.path.join(results, 'AMAP', sample, f'{sample}.AMAP.tsv')
        head, rows = read_tsv(amap)
        col = {c: i for i, c in enumerate(head)}
        anchors = defaultdict(list)
        for r in rows:
            rvs = [x.strip() for x in r[col['rvs']].split(',') if x.strip()]
            if window & set(rvs):
                anchors[r[col['contig']]].append((int(r[col['start']]), int(r[col['end']])))
        intervals = []
        for contig, spans in anchors.items():
            spans.sort()
            a, b = spans[0]
            for s, e in spans[1:]:
                if s - b > GAP:
                    intervals.append((contig, a - MARGIN, b + MARGIN))
                    a = s
                b = max(b, e)
            intervals.append((contig, a - MARGIN, b + MARGIN))
        kept = [r for r in rows if overlaps(intervals, r[col['contig']], int(r[col['start']]), int(r[col['end']]))]
        os.makedirs(os.path.join(out, 'amap'), exist_ok=True)
        write_tsv(os.path.join(out, 'amap', f'{sample}.AMAP.tsv'), head, kept)
        final_ids = {r[col['id']] for r in kept}
        kept_rvs = {x.strip() for r in kept for x in r[col['rvs']].split(',') if x.strip()}

        # Case tables: the calls in the same stretches
        os.makedirs(os.path.join(out, 'cases'), exist_ok=True)
        counts, tags = {}, {}
        for case in '1234':
            h, rs = read_tsv(os.path.join(results, 'CDS_CASES', sample, f'{sample}-case{case}.tsv'))
            c = {x: i for i, x in enumerate(h)}
            rs = [r for r in rs if overlaps(intervals, r[c['Contig']], int(r[c['Start']]), int(r[c['End']]))]
            write_tsv(os.path.join(out, 'cases', f'{sample}-case{case}.tsv'), h, rs)
            counts[case] = len(rs)
            tags[case] = {r[c['Locus Tag']] for r in rs}
        h, rs = read_tsv(os.path.join(results, 'CDS_CASES', sample, f'{sample}-rvmissing.tsv'))
        c = {x: i for i, x in enumerate(h)}
        missing = [r for r in rs if HEADER.sub(r'\1', r[c['Prot Id']]) in kept_rvs | window]
        write_tsv(os.path.join(out, 'cases', f'{sample}-rvmissing.tsv'), h, missing)
        # The summary counts what was kept. Which case 1 calls the filter kept
        # is not written down; those still in the final table are counted.
        kept1 = len(tags['1'] & final_ids)
        n = {'case2': counts['2'], 'case3': counts['3'], 'case4': counts['4'],
             'case1_initial': counts['1'], 'rv_missing': len(missing),
             'total_cds': len(set().union(*tags.values())),
             'case1_removed': counts['1'] - kept1, 'case1_kept': kept1}
        h, rs = read_tsv(os.path.join(results, 'CDS_CASES', 'summary', f'{sample}-summary.tsv'))
        c = {x: i for i, x in enumerate(h)}
        for r in rs:
            r[c['Number']] = str(n[r[c['Type']]])
        write_tsv(os.path.join(out, 'cases', f'{sample}-summary.tsv'), h, rs)

        # Bakta: CDS rows of the stretches, the summary and the GFF3 header
        bdir = os.path.join(out, 'bakta', sample)
        os.makedirs(bdir, exist_ok=True)
        src = os.path.join(results, 'BAKTA', sample, sample)
        with open(src + '.tsv') as fh, open(os.path.join(bdir, f'{sample}.tsv'), 'w') as o:
            for line in fh:
                if line.startswith('#'):
                    o.write(line)
                    continue
                f = line.rstrip('\n').split('\t')
                if f[1] == 'cds' and f[5] in final_ids:
                    f[8] = ''
                    o.write('\t'.join(f) + '\n')
        shutil.copy(src + '.txt', os.path.join(bdir, f'{sample}.txt'))
        with open(src + '.gff3') as fh, open(os.path.join(bdir, f'{sample}.gff3'), 'w') as o:
            for line in fh:
                if line.startswith('##FASTA'):
                    break
                if line.startswith('#'):
                    o.write(line)

    # The reference: the genes of the stretches
    with open(os.path.join(out, 'reference.faa'), 'w') as o:
        for header in order:
            if HEADER.match(header).group(1) in window:
                o.write('>' + header + '\n' + ''.join(proteins[header]))
    with open(mycobrowser, newline='') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    with open(os.path.join(out, 'mycobrowser.txt'), 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(MYCO_COLUMNS)
        for r in rows:
            if r['Locus'] in window:
                w.writerow([r.get(k, '') for k in MYCO_COLUMNS])
    for name in ('software_versions.tsv', 'params.yml'):
        with open(os.path.join(results, 'pipeline_info', name)) as fh:
            text = fh.read()
        # Where the run was made says nothing the tests need
        text = re.sub(r'/\S*/amap-report-run', '/data/run', text)
        text = re.sub(r'/\S*/AMAP-[\w-]+', '/data/AMAP', text)
        with open(os.path.join(out, name), 'w') as fh:
            fh.write(text)
    print(f'{len(window)} reference genes, {len(SAMPLES)} genomes -> {out}')


if __name__ == '__main__':
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
