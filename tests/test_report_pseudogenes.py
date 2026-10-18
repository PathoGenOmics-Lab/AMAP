#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The report's pseudogene chart counts the pseudogenes AMAP called, sample by
sample.

Its Pseudogene Detection section kept the master table's genes flagged in
Mycobrowser's Is_Pseudogene or in AMAP's `pseudogene` column and counted them
by lineage, as "Pseudogenes per Lineage (Mycobrowser annotation)". Is_Pseudogene
says which genes are pseudogenes in H37Rv, not in the sample: CDC1551 came out
at 129, where AMAP calls 123 and its QC table says 123, the other 6 being genes
Mycobrowser marks that are intact in CDC1551. Without Mycobrowser the chart
showed AMAP's calls under Mycobrowser's name. Counting the flagged rows left
out a lineage with none, and adding up a lineage's samples made one with more
samples look as if it had more pseudogenes. REPORT_INPUTS now hands the report
the master tables with the sample each row is from, and the chart counts
AMAP's own call in each sample, a sample with none at 0.

The tests run the steps after the correction as tests/test_sample_ids.py
chains them, with FIX_CDS calling a different number of pseudogenes in each
sample, none in one, and Mycobrowser marking a gene that every sample has
intact, as tests/test_qc_pseudogenes.py marks it, with Mycobrowser on and off;
then the report's own R code for the section on what REPORT_INPUTS made, where
R is installed.
"""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'bin'))

from tests import test_qc_pseudogenes as qc  # noqa: E402
from tests import test_sample_ids as chain  # noqa: E402
from tests.test_sample_ids import rscript, tools  # noqa: E402,F401  (fixtures)

# (the name a genome's files carry, its sample ID, its lineage). The k-th has
# Rv0001, Rv0002, Rv0003 and 2 + k copies of a gene of no reference.
SAMPLES = [('L2_Beijing_1', 'Beijing_1', 'L2'),
           ('L2_Beijing_2', 'Beijing_2', 'L2'),
           ('L4_H37Rv_test', 'H37Rv_test', 'L4'),
           ('NA_K12_self', 'K12_self', 'NA')]          # typing found no lineage

# How many genes FIX_CDS calls pseudogenes in each. L2's two samples have 4
# together, more than the 3 of lineage NA's one sample, which has more than
# either of them; L4's one sample has none.
PSEUDOGENES = {'L2_Beijing_1': 2, 'L2_Beijing_2': 2, 'L4_H37Rv_test': 0, 'NA_K12_self': 3}
REASONS = {'Rv0002': 'frameshift or lost stop codon: its frame runs on 40 residues '
                     'past the end of the reference protein',
           'Rv0003': 'no stop codon in its frame within 300 bp',
           'novel': '1 internal stop codon(s)'}

# tests/test_sample_ids.py's writer, kept before the runs patch it
_write_sample = chain.write_sample


def write_sample(d, name, copies):
    """tests/test_sample_ids.py's sample, with FIX_CDS calling as many of its
    genes pseudogenes as PSEUDOGENES says, never Rv0001, which Mycobrowser
    marks, nor the first copy of the gene of no reference, with the reason in
    the note."""
    _write_sample(d, name, copies)
    table = d / f'{name}.AMAP.tsv'
    header, *rows = table.read_text().splitlines()
    flagged = novel = 0
    for i, row in enumerate(rows):
        cells = row.split('\t')
        key = cells[4] or 'novel'
        novel += key == 'novel'
        if flagged < PSEUDOGENES[name] and key != 'Rv0001' and novel != 1:
            cells[5], cells[8] = 'yes', REASONS[key]
            rows[i] = '\t'.join(cells)
            flagged += 1
    assert flagged == PSEUDOGENES[name]
    table.write_text('\n'.join([header] + rows) + '\n')
    return d


RUNS = ['mycobrowser', 'no_mycobrowser']


@pytest.fixture(scope='module', params=RUNS)
def run(request, tmp_path_factory, tools):
    """(the tasks of a run, whether Mycobrowser was on) for both runs."""
    mycobrowser = request.param == 'mycobrowser'
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(chain, 'write_sample', write_sample)
        mp.setattr(chain, 'write_mycobrowser', qc.write_mycobrowser)
        if not mycobrowser:
            mp.setattr(chain, 'run_task', qc.skipping_mycobrowser(chain.run_task))
        tasks = chain.run_pipeline(tmp_path_factory.mktemp(request.param), SAMPLES, tools)
    return tasks, mycobrowser


def section(tasks, rscript, tables, replace=None):
    """Run the report's pseudogene section on a run's report inputs, as
    tests/test_sample_ids.py runs the report's chunks, and read back `tables`
    ({R expression: name})."""
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(chain, 'R_CHUNKS', ['setup', 'pseudogenes'])
        mp.setattr(chain, 'R_TABLES', tables)
        return chain.report_tables(tasks, rscript, replace=replace)


@pytest.fixture(scope='module')
def chart(run, rscript):
    """The count the section draws for each sample, and the chart's labels."""
    tasks, _ = run
    return section(tasks, rscript, {
        'pseudo_count': 'counts',
        'data.frame(title = p$labels$title, x = p$labels$x, y = p$labels$y)': 'labels'})


def per_sample(counts):
    """{(lineage, sample): pseudogenes} from the counts the section draws."""
    assert {'Lineage', 'Sample', 'Pseudogenes'} <= set(counts.columns), \
        f'the chart counts {list(counts.columns)}, not each sample'
    return dict(zip(zip(counts['Lineage'], counts['Sample']), counts['Pseudogenes']))


def amap_calls():
    """{(lineage, sample): pseudogenes} as FIX_CDS called them."""
    return {(lineage, name): PSEUDOGENES[name] for name, _, lineage in SAMPLES}


class TestReportInputs:

    def test_each_master_table_row_says_its_sample(self, run):
        # the name its file carries, as the QC table's rows say it
        tasks, _ = run
        header, *rows = (tasks['report'] / 'input_mastertables.txt').read_text().splitlines()
        for name, _, _ in SAMPLES:
            own = (tasks['outputs'][name] / f'{name}-AMAPmastertable.tsv').read_text().splitlines()
            assert header == 'sample\t' + own[0]
            assert [r.split('\t', 1)[1] for r in rows if r.split('\t', 1)[0] == name] == own[1:]
        assert len(rows) == sum(
            len((tasks['outputs'][name] / f'{name}-AMAPmastertable.tsv').read_text().splitlines()) - 1
            for name, _, _ in SAMPLES)


class TestChart:

    def test_counts_the_pseudogenes_amap_called_in_each_sample(self, chart):
        # Main counted by lineage, Mycobrowser's Rv0001 in every sample with
        # it on: L2 6, NA 4, L4 1; and without it L2 4, NA 3, and no L4
        assert per_sample(chart['counts']) == amap_calls()

    def test_a_sample_with_none_counts_zero(self, chart):
        counts = per_sample(chart['counts'])
        assert counts[('L4', 'L4_H37Rv_test')] == 0

    def test_as_many_as_the_qc_table_says(self, run, chart):
        tasks, _ = run
        qc_table = chain.read_tsv(tasks['report'] / 'input_qc.txt')
        assert {sample: n for (_, sample), n in per_sample(chart['counts']).items()} == \
            dict(zip(qc_table['sample'], qc_table['pseudogene_count']))

    def test_the_title_names_amap_s_calls_per_sample(self, chart):
        labels = chart['labels'].iloc[0]
        assert 'Mycobrowser' not in labels['title']
        assert 'AMAP' in labels['title'] and 'per Sample' in labels['title']
        assert (labels['x'], labels['y']) == ('Lineage', 'Pseudogenes per sample')

    def test_every_sample_at_zero_when_none_has_any(self, run, rscript):
        # no chart then, but the table still lists every sample, at 0
        tasks, _ = run
        header, *rows = (tasks['report'] / 'input_mastertables.txt').read_text().splitlines()
        col = header.split('\t').index('pseudogene')
        intact = []
        for row in rows:
            cells = row.split('\t')
            cells[col] = 'no'
            intact.append('\t'.join(cells))
        tables = section(tasks, rscript, {'pseudo_count': 'counts'},
                         replace={'input_mastertables.txt': '\n'.join([header] + intact) + '\n'})
        assert per_sample(tables['counts']) == {key: 0 for key in amap_calls()}
