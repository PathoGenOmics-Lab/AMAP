#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The QC table counts the pseudogenes AMAP called.

annotation_qc.py took the pseudogene count of `*-AMAP.qc.tsv` from the first
column the master table has of Is_Pseudogene, pseudogene and is_pseudogene.
Every master table has Mycobrowser's Is_Pseudogene, empty when Mycobrowser is
off, so every organism-agnostic run reported 0 pseudogenes (E. coli Sakai
against K-12: 0, where AMAP called 129), and an MTBC run counted the sample's
genes named after one of the 13 genes Mycobrowser marks as pseudogenes in H37Rv
(CDC1551: 13, where AMAP called 123 and finds 6 of the 13 intact). The count
is now AMAP's own call, the master table's `pseudogene` column, which the
GFF3's pseudo=true and the GenBank's /pseudo carry.

The tests run the steps after the correction as Nextflow renders them, chained
as tests/test_sample_ids.py chains them, on its genomes with FIX_CDS calling
some of their genes pseudogenes and Mycobrowser marking another one, with
Mycobrowser on and off; and read the report's QC table, built by its own R
code from what REPORT_INPUTS hands over, where R is installed.
"""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'bin'))

from tests import test_sample_ids as chain  # noqa: E402
from tests.test_sample_ids import rscript, tools  # noqa: E402,F401  (fixtures)

SAMPLES = [('L2_Beijing_1', 'Beijing_1', 'L2'), ('L4_H37Rv_test', 'H37Rv_test', 'L4')]

# Each genome has Rv0001, Rv0002, Rv0003 and 2 + k copies of a gene of no
# reference, k its place. FIX_CDS calls Rv0002 a pseudogene, and every copy of
# that gene but the first, so the k-th sample has 2 + k.
REASONS = {'Rv0002': 'frameshift or lost stop codon: its frame runs on 40 residues '
                     'past the end of the reference protein',
           'novel': '1 internal stop codon(s)'}


def pseudogenes(k):
    """How many pseudogenes AMAP calls in the k-th sample."""
    return 2 + k


# tests/test_sample_ids.py's writers, kept before the runs patch them
_write_sample, _write_mycobrowser = chain.write_sample, chain.write_mycobrowser


def write_sample(d, name, copies):
    """tests/test_sample_ids.py's sample, with FIX_CDS's pseudogenes and their
    reasons in its table."""
    _write_sample(d, name, copies)
    table = d / f'{name}.AMAP.tsv'
    header, *rows = table.read_text().splitlines()
    novel = 0
    for i, row in enumerate(rows):
        cells = row.split('\t')
        key = cells[4] or 'novel'
        novel += key == 'novel'
        if key == 'Rv0002' or (key == 'novel' and novel > 1):
            cells[5], cells[8] = 'yes', REASONS[key]
            rows[i] = '\t'.join(cells)
    table.write_text('\n'.join([header] + rows) + '\n')
    return d


def write_mycobrowser(path):
    """Mycobrowser's table, which marks Rv0001 as a pseudogene of H37Rv, a gene
    AMAP finds intact in every sample here, as it finds 6 of the 13 in CDC1551."""
    _write_mycobrowser(path)
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines):
        cells = line.split('\t')
        if cells[8] == 'Rv0001':
            cells[14] = 'Yes'
            lines[i] = '\t'.join(cells)
    path.write_text('\n'.join(lines) + '\n')


def skipping_mycobrowser(run_task):
    """run_task, with GET_OUTPUTS run under --skip_mycobrowser, as a run with a
    reference from another genus is."""
    def run(work, nf_file, process, values, inputs, tools):
        if process == 'GET_OUTPUTS':
            values = dict(values, **{'params.skip_mycobrowser': 'true'})
        return run_task(work, nf_file, process, values, inputs, tools)
    return run


RUNS = ['mycobrowser', 'no_mycobrowser']


@pytest.fixture(scope='module', params=RUNS)
def run(request, tmp_path_factory, tools):
    """(the tasks of a run, whether Mycobrowser was on) for both runs."""
    mycobrowser = request.param == 'mycobrowser'
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(chain, 'write_sample', write_sample)
        mp.setattr(chain, 'write_mycobrowser', write_mycobrowser)
        if not mycobrowser:
            mp.setattr(chain, 'run_task', skipping_mycobrowser(chain.run_task))
        tasks = chain.run_pipeline(tmp_path_factory.mktemp(request.param), SAMPLES, tools)
    return tasks, mycobrowser


@pytest.fixture(scope='module')
def report(run, rscript):
    """The report's QC table, as its own R code reads REPORT_INPUTS' table."""
    tasks, _ = run
    with pytest.MonkeyPatch.context() as mp:
        # only the QC section: without Mycobrowser the PE/PPE one builds no table
        mp.setattr(chain, 'R_CHUNKS', ['setup', 'qc_table'])
        mp.setattr(chain, 'R_TABLES', {'qc_data': 'qc'})
        return chain.report_tables(tasks, rscript)['qc']


def qc_table(tasks, name):
    return chain.read_tsv(tasks['outputs'][name] / f'{name}-AMAP.qc.tsv')


class TestGetOutputs:

    def test_the_master_table_has_both_columns(self, run):
        # as every master table has: Mycobrowser's, empty when it is off, and AMAP's
        tasks, mycobrowser = run
        for k, (name, _, _) in enumerate(SAMPLES):
            master = chain.read_tsv(tasks['outputs'][name] / f'{name}-AMAPmastertable.tsv')
            marked = master.loc[master['Is_Pseudogene'] != '', 'Locus']
            assert dict(zip(marked, master.loc[marked.index, 'Is_Pseudogene'])) == \
                ({'Rv0001': 'Yes', 'Rv0002': 'No', 'Rv0003': 'No'} if mycobrowser else {})
            assert (master['pseudogene'] == 'yes').sum() == pseudogenes(k)

    def test_the_qc_table_counts_the_pseudogenes_amap_called(self, run):
        # 1 each with Mycobrowser on (Rv0001, intact in every sample) and 0 without, before
        tasks, _ = run
        for k, (name, _, _) in enumerate(SAMPLES):
            assert qc_table(tasks, name)['pseudogene_count'].tolist() == [pseudogenes(k)]

    def test_as_many_as_the_gff3_and_the_genbank_mark(self, run):
        tasks, _ = run
        for name, _, _ in SAMPLES:
            out = tasks['outputs'][name]
            features = [line.split('\t') for line in
                        (out / f'{name}-AMAP.gff3').read_text().split('##FASTA')[0].splitlines()
                        if line and not line.startswith('#')]
            gff3 = sum(f[2] == 'gene' and 'pseudo=true' in f[8].split(';') for f in features)
            genbank = (out / f'{name}-AMAP.gbk').read_text().count('/pseudo=""')
            assert qc_table(tasks, name)['pseudogene_count'].tolist() == [gff3] == [genbank]


class TestReport:

    def test_report_inputs_hands_over_each_sample_s_count(self, run):
        tasks, _ = run
        qc = chain.read_tsv(tasks['report'] / 'input_qc.txt')
        assert dict(zip(qc['sample'], qc['pseudogene_count'])) == \
            {name: pseudogenes(k) for k, (name, _, _) in enumerate(SAMPLES)}

    def test_the_report_s_qc_table_shows_it(self, report):
        assert dict(zip(report['sample'], report['pseudogene_count'])) == \
            {name: pseudogenes(k) for k, (name, _, _) in enumerate(SAMPLES)}
