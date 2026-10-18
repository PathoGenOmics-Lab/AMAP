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

count_annot.py writes a row for each category a sample has genes in and none
for the others, and the report plotted the rows as they came: a sample with no
gene in a category was missing from that category's figures, from the Dunn
tests run on them and from the statistics table, where it should count 0, so
its lineage's box and mean were those of the samples that have some. The
report now counts each sample 0 in every category another sample of the run
has genes in, in the COG and in the Mycobrowser tables.

The tests count columns shaped like a master table's, then run the steps after
the correction as Nextflow renders them, chained as tests/test_sample_ids.py
chains them, with eggNOG giving some genes no category, then no gene one, then
some samples no gene in a category another sample has genes in, and the
report's own R code on what REPORT_INPUTS hands over, where R is installed.
"""

import sys
from pathlib import Path
from statistics import median

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

# Two samples of L2, so that one of them can lack a category the other has
SAMPLES = [('L2_Beijing_1', 'Beijing_1', 'L2'), ('L2_Beijing_2', 'Beijing_2', 'L2'),
           ('L4_H37Rv_test', 'H37Rv_test', 'L4')]

# eggNOG's category for each gene, the same for every sample or one per sample,
# and what each sample's genes count as then, the k-th sample having 2 + k
# copies of the gene of no reference
RUNS = {
    # Rv0002 without a category, and Rv0001 in two
    'some_genes': ({'Rv0001': 'JK', 'Rv0002': '-', 'Rv0003': 'K', 'novel': 'X'},
                   lambda k: {'J': 1, 'K': 2, 'NaN': 1, 'X': 2 + k}),
    # no gene with one, as in that test-profile run
    'no_gene': (dict.fromkeys(['Rv0001', 'Rv0002', 'Rv0003', 'novel'], '-'),
                lambda k: {'NaN': 5 + k}),
    # Beijing_2 with no gene in J, which Beijing_1 has, the L2 samples with
    # none in S and H37Rv_test with none in None
    'samples_lack_some': ({'Rv0001': ('JK', 'K', 'JK'), 'Rv0002': ('-', '-', 'S'),
                           'Rv0003': 'K', 'novel': 'X'},
                          lambda k: [{'J': 1, 'K': 2, 'NaN': 1, 'X': 2},
                                     {'K': 2, 'NaN': 1, 'X': 3},
                                     {'J': 1, 'K': 2, 'S': 1, 'X': 4}][k]),
}


def categories(cog, k):
    """eggNOG's category for each gene of the k-th sample."""
    return {gene: c[k] if isinstance(c, tuple) else c for gene, c in cog.items()}


def shown(expected):
    """What the report counts each sample's genes as: every category some
    sample of the run has genes in, 0 where that sample has none."""
    each = [expected(k) for k in range(len(SAMPLES))]
    return [{c: genes.get(c, 0) for c in set().union(*each)} for genes in each]


def as_in_report(genes):
    """The categories as the report names them, 'NaN' being its "NA"."""
    return {'NA' if c == 'NaN' else c: n for c, n in genes.items()}


def by_lineage(counts):
    """{(lineage, category): (samples, their mean share of it)} of each
    sample's counts, a share being of all the sample's counts as in the
    report."""
    shares = {}
    for (_, _, lineage), genes in zip(SAMPLES, counts):
        for c, n in genes.items():
            shares.setdefault((lineage, c), []).append(n / sum(genes.values()))
    return {key: (len(s), pytest.approx(sum(s) / len(s), abs=1e-4)) for key, s in shares.items()}


def statistics(table):
    """{(lineage, category): (Count, Mean_Pct)} of a statistics table."""
    return dict(zip(zip(table['Lineages'], table['Functional_Category']),
                    zip(table['Count'], table['Mean_Pct'])))


@pytest.fixture(scope='module', params=sorted(RUNS))
def run(request, tmp_path_factory, tools):
    """(the tasks of a run, what each sample's genes count as) for each run."""
    cog, expected = RUNS[request.param]
    names = [name for name, _, _ in SAMPLES]
    write_sample = chain.write_sample

    def write_sample_with_its_categories(d, name, copies):
        with pytest.MonkeyPatch.context() as mp:
            mp.setattr(chain, 'COG', categories(cog, names.index(name)))
            return write_sample(d, name, copies)

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(chain, 'write_sample', write_sample_with_its_categories)
        tasks = chain.run_pipeline(tmp_path_factory.mktemp(request.param), SAMPLES, tools)
    return tasks, expected


# What each of the report's per-COG figures plots, by the report's own
# create_scatter_plot, as one more table for report_tables to write out. The
# figure's Dunn tests between lineages run on the same rows.
FIGURES = ('do.call(rbind, lapply(names(cog_names), function(code) {'
           ' d <- create_scatter_plot(datos_cog1, code, "Count", "Number of CDS",'
           ' cog_names[code])$data;'
           ' data.frame(figure=rep(cog_names[[code]], nrow(d)), genome=d$genome,'
           ' Count=d$Count) }))')

# The figure each category of the runs is drawn in
FIGURE = {'J': 'J - Translation, ribosomal structure', 'K': 'K - Transcription',
          'S': 'S - Function unknown', 'X': 'X - Mobilome: prophages, transposons',
          'NaN': 'None'}

# The category the overview plots each row under, R's NA written as "lost",
# since write.table writes it as "NA", the category None's own
OVERVIEW = ('data.frame(genome=datos_cog1$genome, Count=datos_cog1$Count,'
            ' category=ifelse(is.na(datos_cog1$Functional_Category), "lost",'
            ' as.character(datos_cog1$Functional_Category)))')


@pytest.fixture(scope='module')
def report(run, rscript):
    """(the report's COG table, {(figure, genome): count} of its per-COG
    figures, the rows of its overview, its statistics table, what each
    sample's genes count as) for a run."""
    tasks, expected = run
    with pytest.MonkeyPatch.context() as mp:
        mp.setitem(chain.R_TABLES, FIGURES, 'cog_figures')
        mp.setitem(chain.R_TABLES, OVERVIEW, 'cog_overview')
        mp.setitem(chain.R_TABLES, 'stats_cog', 'cog_statistics')
        tables = chain.report_tables(tasks, rscript)
    plotted = tables['cog_figures']
    figures = dict(zip(zip(plotted['figure'], plotted['genome']), plotted['Count']))
    return tables['cog'], figures, tables['cog_overview'], tables['cog_statistics'], expected


class TestGetCount:

    def test_category_counts(self, run):
        # a row for each category the sample has genes in, and none for the
        # others: the report fills those in
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
        # 0 in each category another sample has genes in
        cog, _, _, _, expected = report
        for (_, sample, _), genes in zip(SAMPLES, shown(expected)):
            rows = cog[cog['genome'] == sample]
            assert dict(zip(rows['Functional_Category'], rows['Count'])) == as_in_report(genes)
            assert rows['Percentage'].sum() == pytest.approx(1.0)

    def test_the_figure_none_counts_them(self, report):
        # H37Rv_test, every gene of which has a category, at 0
        _, figures, _, _, expected = report
        for (_, sample, _), genes in zip(SAMPLES, shown(expected)):
            assert figures.get(('None', sample)) == genes.get('NaN')

    def test_each_category_in_a_figure_of_its_own(self, report):
        # the figure X - Mobilome plotted none of the X genes, before. A
        # sample with no gene in a category is plotted at 0, where Beijing_2
        # was missing from L2's box in J, and a category no sample has genes
        # in has a figure with no point.
        _, figures, _, _, expected = report
        for (_, sample, _), genes in zip(SAMPLES, shown(expected)):
            assert {figure: n for (figure, genome), n in figures.items() if genome == sample} == \
                {FIGURE[c]: n for c, n in genes.items()}

    def test_the_overview_plots_each_category_as_itself(self, report):
        # and the overview plotted them as NA
        _, _, overview, _, expected = report
        for (_, sample, _), genes in zip(SAMPLES, shown(expected)):
            rows = overview[overview['genome'] == sample]
            assert dict(zip(rows['category'], rows['Count'])) == as_in_report(genes)

    def test_the_statistics_count_every_sample(self, report):
        # all the lineage's samples, and a mean share with 0 for those with no
        # gene in the category: L2's J counted Beijing_1 alone, at twice the
        # mean, and L2's S and L4's None had no row
        _, _, _, table, expected = report
        assert statistics(table) == by_lineage(as_in_report(genes) for genes in shown(expected))


# ---------------------------------------------------------------------------
# The report's Mycobrowser tables. Every run gives the reference genes the
# same categories, one each, and the gene of no reference none, which the
# report calls "NA": a sample with no gene in one is the run's counts without
# its row, as REPORT_INPUTS hands over GET_COUNT's tables.
# ---------------------------------------------------------------------------

LACKING = ('Beijing_2', 'conserved hypotheticals')

# What each per-category figure plots, the categories as the report lists them
MYCO_FIGURES = ('do.call(rbind, lapply(setdiff(unique(na.omit(datos_myco$Functional_Category)), "NA"),'
                ' function(category) { d <- create_scatter_plot(datos_myco, category, "Count",'
                ' "Number of CDS", category)$data; data.frame(figure=rep(category, nrow(d)),'
                ' genome=d$genome, Count=d$Count) }))')


def myco_counts():
    """What each sample's genes count as in Mycobrowser's categories, with
    the sample in LACKING at 0 in its category."""
    each = []
    for k, (_, sample, _) in enumerate(SAMPLES):
        genes = dict.fromkeys(chain.CATEGORY.values(), 1)
        genes['NA'] = 2 + k
        if sample == LACKING[0]:
            genes[LACKING[1]] = 0
        each.append(genes)
    return each


@pytest.fixture(scope='module')
def myco(run, rscript):
    """The report's Mycobrowser tables, its per-category figures as
    {(figure, genome): count}, when Beijing_2 has no gene in 'conserved
    hypotheticals'."""
    tasks, _ = run
    lines = (tasks['report'] / 'input_myco.txt').read_text().splitlines(keepends=True)
    kept = [line for line in lines if line.split('\t')[1:3] != list(LACKING)]
    assert len(kept) == len(lines) - 1
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(chain, 'R_CHUNKS', chain.R_CHUNKS + ['myco_stats'])
        mp.setitem(chain.R_TABLES, MYCO_FIGURES, 'myco_figures')
        mp.setitem(chain.R_TABLES, 'stats_myco', 'myco_statistics')
        mp.setitem(chain.R_TABLES, 'merged_data', 'myco_ratios')
        tables = chain.report_tables(tasks, rscript, replace={'input_myco.txt': ''.join(kept)})
    plotted = tables['myco_figures']
    tables['myco_figures'] = dict(zip(zip(plotted['figure'], plotted['genome']), plotted['Count']))
    return tables


class TestMycobrowser:

    def test_category_counts(self, myco):
        counts = myco['myco']
        for (_, sample, _), genes in zip(SAMPLES, myco_counts()):
            rows = counts[counts['genome'] == sample]
            assert dict(zip(rows['Functional_Category'], rows['Count'])) == genes
            assert rows['Percentage'].sum() == pytest.approx(1.0)

    def test_each_category_in_a_figure_of_its_own(self, myco):
        # Beijing_2 at 0 in 'conserved hypotheticals', where it was missing
        # from L2's box. Genes with no category have no figure.
        figures = myco['myco_figures']
        for (_, sample, _), genes in zip(SAMPLES, myco_counts()):
            assert {figure: n for (figure, genome), n in figures.items() if genome == sample} == \
                {c: n for c, n in genes.items() if c != 'NA'}

    def test_the_statistics_count_every_sample(self, myco):
        assert statistics(myco['myco_statistics']) == by_lineage(myco_counts())

    def test_the_ratio_to_the_median_of_every_sample(self, myco):
        # Beijing_2 at 0 in 'conserved hypotheticals', against the median of
        # all three samples, where it was missing from the comparison
        each = myco_counts()
        ratios = myco['myco_ratios']
        assert dict(zip(zip(ratios['genome'], ratios['Functional_Category']),
                        ratios['Ratio_Percentage'])) == \
            {(sample, c): pytest.approx(100 * n / round(median(g[c] for g in each)))
             for (_, sample, _), genes in zip(SAMPLES, each) for c, n in genes.items()}
