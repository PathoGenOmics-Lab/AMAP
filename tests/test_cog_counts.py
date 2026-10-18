#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A gene eggNOG gives no COG category counts once in the COG counts.

count_annot.py counts a gene once in each letter of its COG category, "JK" in
J and in K. A gene with none, which eggNOG writes as '-', was filled in as
'Unassigned' and split the same way, so it counted once in A, D, E, G, I and U
and twice in N and S: a `-profile docker,test` run with a stand-in for eggNOG
that gave none of its 4036 genes a category reported 4036 genes in each of A,
D, E, G, I and U and 8072 in N and S. Such a gene now counts once, as 'NaN',
which the report reads as its category "NA", shown as "None". Percentages are
still of all the genome's genes.

The report's COG figures took their categories from a list that ran A to V,
without R, and None: a gene of R, W, X, Y or Z was plotted as NA in the
overview and in no figure of its own, so the figure X - Mobilome came out
empty. They now come from the report's own list of figures.

The tests count columns shaped like a master table's, then run the steps after
the correction as Nextflow renders them, chained as tests/test_sample_ids.py
chains them, with eggNOG giving some genes no category and then no gene one,
and the report's own R code on what REPORT_INPUTS hands over, where R is
installed.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'bin'))

from count_annot import process_cog_categories  # noqa: E402
from tests import test_sample_ids as chain  # noqa: E402
from tests.test_sample_ids import rscript, tools  # noqa: E402,F401  (fixtures)


def counts(categories):
    """{category: (genes, percentage of the genome's genes)} of a genome
    whose genes have these COG categories."""
    df = pd.DataFrame({'lineage': 'L4', 'genome': 'H37Rv', 'COG_category': categories})
    cog = process_cog_categories(df)
    return {c: (n, round(pct, 2)) for c, n, pct in
            zip(cog['COG_category'], cog['count'], cog['percentage_of_genes'])}


class TestCounts:

    def test_a_gene_without_a_category_counts_once(self):
        # A, D, E, G, I, U and X 1, and K, N and S 2, before
        assert counts(['K', 'K', 'X', np.nan]) == {'K': (2, 50.0), 'NaN': (1, 25.0), 'X': (1, 25.0)}

    def test_no_gene_with_a_category(self):
        # eggNOG's '-' everywhere, which count_annot.py reads as missing: a
        # column with no string in it, float64 under pandas 2
        assert counts(pd.Series([np.nan] * 4)) == {'NaN': (4, 100.0)}

    def test_a_gene_in_two_categories_counts_in_both(self):
        # so the percentages of the genome's genes add up to more than 100
        assert counts(['JK', 'K', 'S']) == {'J': (1, 33.33), 'K': (2, 66.67), 'S': (1, 33.33)}

    def test_only_letters_are_categories(self):
        # lowercase is its letter, a letter counts once a gene, and a gene
        # with no letter has no category
        assert counts(['jk', 'KK', '-', '']) == {'J': (1, 25.0), 'K': (2, 50.0), 'NaN': (2, 50.0)}


# ---------------------------------------------------------------------------
# GET_COUNT and the report on what eggNOG gave the genes of
# tests/test_sample_ids.py's genomes: Rv0001, Rv0002, Rv0003 and the gene of
# no reference in 2 + k copies, k the sample's place
# ---------------------------------------------------------------------------

SAMPLES = [('L2_Beijing_1', 'Beijing_1', 'L2'), ('L4_H37Rv_test', 'H37Rv_test', 'L4')]

# eggNOG's category for each gene, and what each sample's genes count as then,
# the k-th sample having 2 + k copies of the gene of no reference
RUNS = {
    # Rv0002 without a category, and Rv0001 in two
    'some_genes': ({'Rv0001': 'JK', 'Rv0002': '-', 'Rv0003': 'K', 'novel': 'X'},
                   lambda k: {'J': 1, 'K': 2, 'NaN': 1, 'X': 2 + k}),
    # no gene with one, as in that test-profile run
    'no_gene': (dict.fromkeys(['Rv0001', 'Rv0002', 'Rv0003', 'novel'], '-'),
                lambda k: {'NaN': 5 + k}),
}


@pytest.fixture(scope='module', params=sorted(RUNS))
def run(request, tmp_path_factory, tools):
    """(the tasks of a run, what each sample's genes count as) for both runs."""
    cog, expected = RUNS[request.param]
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(chain, 'COG', cog)
        tasks = chain.run_pipeline(tmp_path_factory.mktemp(request.param), SAMPLES, tools)
    return tasks, expected


# What each of the report's per-COG figures plots, by the report's own
# create_scatter_plot, as one more table for report_tables to write out
FIGURES = ('do.call(rbind, lapply(names(cog_names), function(code) {'
           ' d <- create_scatter_plot(datos_cog1, code, "Count", "Number of CDS",'
           ' cog_names[code])$data;'
           ' data.frame(figure=rep(cog_names[[code]], nrow(d)), genome=d$genome,'
           ' Count=d$Count) }))')

# The figure each category of the runs is drawn in
FIGURE = {'J': 'J - Translation, ribosomal structure', 'K': 'K - Transcription',
          'X': 'X - Mobilome: prophages, transposons', 'NaN': 'None'}

# The category the overview plots each row under, R's NA written as "lost",
# since write.table writes it as "NA", the category None's own
OVERVIEW = ('data.frame(genome=datos_cog1$genome, Count=datos_cog1$Count,'
            ' category=ifelse(is.na(datos_cog1$Functional_Category), "lost",'
            ' as.character(datos_cog1$Functional_Category)))')


@pytest.fixture(scope='module')
def report(run, rscript):
    """(the report's COG table, {(figure, genome): count} of its per-COG
    figures, the rows of its overview, what each sample's genes count as)
    for a run."""
    tasks, expected = run
    with pytest.MonkeyPatch.context() as mp:
        mp.setitem(chain.R_TABLES, FIGURES, 'cog_figures')
        mp.setitem(chain.R_TABLES, OVERVIEW, 'cog_overview')
        tables = chain.report_tables(tasks, rscript)
    plotted = tables['cog_figures']
    figures = dict(zip(zip(plotted['figure'], plotted['genome']), plotted['Count']))
    return tables['cog'], figures, tables['cog_overview'], expected


class TestGetCount:

    def test_category_counts(self, run):
        tasks, expected = run
        for k, (name, sample, lineage) in enumerate(SAMPLES):
            cog = chain.read_tsv(tasks['count'][name] / f'{name}-AMAP.COG.tsv')
            genes = 5 + k
            assert {c: (n, pct) for c, n, pct in
                    zip(cog['COG_category'], cog['count'], cog['percentage_of_genes'])} == \
                {c: (n, pytest.approx(100 * n / genes)) for c, n in expected(k).items()}
            assert set(zip(cog['lineage'], cog['genome'])) == {(lineage, sample)}


class TestReport:

    def test_category_counts(self, report):
        cog, _, _, expected = report
        for k, (_, sample, _) in enumerate(SAMPLES):
            rows = cog[cog['genome'] == sample]
            # 'NaN' is the report's category "NA"
            assert dict(zip(rows['Functional_Category'], rows['Count'])) == \
                {'NA' if c == 'NaN' else c: n for c, n in expected(k).items()}
            assert rows['Percentage'].sum() == pytest.approx(1.0)

    def test_the_figure_none_counts_them(self, report):
        _, figures, _, expected = report
        for k, (_, sample, _) in enumerate(SAMPLES):
            assert figures.get(('None', sample)) == expected(k)['NaN']

    def test_each_category_in_a_figure_of_its_own(self, report):
        # the figure X - Mobilome plotted none of the X genes, before
        _, figures, _, expected = report
        for k, (_, sample, _) in enumerate(SAMPLES):
            assert {figure: n for (figure, genome), n in figures.items() if genome == sample} == \
                {FIGURE[c]: n for c, n in expected(k).items()}

    def test_the_overview_plots_each_category_as_itself(self, report):
        # and the overview plotted them as NA
        _, _, overview, expected = report
        for k, (_, sample, _) in enumerate(SAMPLES):
            rows = overview[overview['genome'] == sample]
            assert dict(zip(rows['category'], rows['Count'])) == \
                {'NA' if c == 'NaN' else c: n for c, n in expected(k).items()}
