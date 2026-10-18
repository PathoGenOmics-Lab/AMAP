#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
How the duplications are reported: DUPLICATIONS/homologous_communities.tsv and
the report's section on it.

- Modularity1 and Modularity2 were computed with each community as its only
  group, which makes them 0 whatever the community: they are gone.
- Elements1 and Elements2 count a community's proteins by COG category. With
  numpy 2 they read {'L': np.int64(4)}, and rename_fastas.py appends a
  category only to the proteins eggNOG gave one, leaving the others as FIX_CDS
  wrote them, <name>;<strand>;<start>-<end>, so those were counted in
  categories '+' and '-'. They are now plain counts, the largest first, with
  '-', eggNOG's mark for none, for the proteins without a category.
- The report's duplication summary printed its markdown with cat() in a chunk
  without results='asis', so the report showed
  "## **48** homologous community pairs ..." in a code box.

The tests run GET_DUPLICATIONS and GET_ALLDUP as Nextflow renders them,
chained as the pipeline chains them, on tests/test_sample_ids.py's genomes
with eggNOG giving one copy of the duplicated gene a category, another '-' and
the last no row at all; then knit the report's chunk on GET_ALLDUP's table,
where R is installed.
"""

import ast
import itertools
import re
import subprocess
import sys
from io import StringIO
from pathlib import Path

import networkx as nx
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'bin'))

import network_genomes  # noqa: E402
from tests import test_sample_ids as chain  # noqa: E402
from tests.test_sample_ids import rscript, tools  # noqa: E402,F401  (fixtures)

REPORT = ROOT / 'bin' / 'report' / 'report_AMAP.Rmd'

# homologous_communities.tsv's columns, for each of the two communities of a row
COLUMNS = [f'{name}{k}' for k in '12' for name in (
    'Genome', 'Community', 'Size', 'Elements', 'AvgDegree', 'Density',
    'AvgClustering', 'AvgPathLength', 'Nodes')]


# ---------------------------------------------------------------------------
# Elements and the metrics, from network_genomes.py itself
# ---------------------------------------------------------------------------

class TestElements:

    def test_the_category_rename_fastas_appends(self):
        assert network_genomes.cog_category('L4_CDC1551_Rv2885c;L') == 'L'
        # a protein in two categories is in both, as eggNOG writes it
        assert network_genomes.cog_category('L4_CDC1551_Rv0001;KT') == 'KT'

    def test_a_protein_without_a_category_is_not_counted_by_its_strand(self):
        # rename_fastas.py leaves it as FIX_CDS wrote it: '-' and '+' before
        assert network_genomes.cog_category('L4_CDC1551_Rv2885c;-;3188319-3189702') == '-'
        assert network_genomes.cog_category('L4_CDC1551_02906;+;3094920-3096300') == '-'
        assert network_genomes.cog_category('L4_CDC1551_02906') == '-'

    def test_plain_counts_the_largest_first(self):
        nodes = ['S_00001;L', 'S_00002;+;10-300', 'S_00003;K', 'S_00004;L',
                 'S_00005;-;400-700', 'S_00006;L', 'S_00007;S']
        counts = network_genomes.count_categories(set(nodes))
        assert str(counts) == "{'L': 3, '-': 2, 'K': 1, 'S': 1}"
        assert all(type(n) is int for n in counts.values())

    def test_the_order_does_not_depend_on_the_order_of_the_proteins(self):
        # the proteins come as a set, whose order changes from one run to the
        # next; ties came out in that order
        nodes = ['S_00001;L', 'S_00002;K', 'S_00003;+;10-300', 'S_00004;S']
        written = {str(network_genomes.count_categories(p)) for p in itertools.permutations(nodes)}
        assert written == {"{'-': 1, 'K': 1, 'L': 1, 'S': 1}"}


class TestMetrics:

    def test_no_modularity(self):
        # of a community taken as its only group, it is 0 whatever its shape
        for graph in (nx.complete_graph(4), nx.path_graph(4), nx.star_graph(3)):
            metrics = network_genomes.calculate_subgraph_metrics(graph)
            assert set(metrics) == {'avg_degree', 'density', 'avg_clustering', 'avg_path_length'}

    def test_a_row_per_pair_of_communities(self):
        nodes = {'A': ['A_00001;L', 'A_00002;-;1-300', 'A_00003;L'],
                 'B': ['B_00001;L', 'B_00002;L']}
        networks = {g: nx.complete_graph(n) for g, n in nodes.items()}
        communities = {g: {'Community_0': set(n)} for g, n in nodes.items()}
        homologous = {('A', 'Community_0'): [('B', 'Community_0')],
                      ('B', 'Community_0'): [('A', 'Community_0')]}
        report = network_genomes.create_summary_report(communities, homologous, networks)
        assert report.columns.tolist() == COLUMNS
        assert len(report) == 1
        assert str(report.loc[0, 'Elements1']) == "{'L': 2, '-': 1}"
        assert str(report.loc[0, 'Elements2']) == "{'L': 2}"


# ---------------------------------------------------------------------------
# GET_DUPLICATIONS and GET_ALLDUP, as Nextflow renders them, on
# tests/test_sample_ids.py's genomes: Rv0001, Rv0002, Rv0003 and copies of a
# gene of no reference, the duplication
# ---------------------------------------------------------------------------

# eggNOG's category for each copy of the duplicated gene, in the genome's
# order (None: eggNOG has no row for it), and the Elements they make
SAMPLES = {
    'L2_Beijing_1': (['X', '-', None], "{'-': 2, 'X': 1}"),
    'L4_H37Rv_test': (['X', 'X'], "{'X': 2}"),
}


def write_eggnog(path, categories):
    """eggNOG's annotation as chain.write_sample writes it, with each copy of
    the duplicated gene given its own category from `categories`."""
    lines, copies = [], iter(categories)
    for line in path.read_text().splitlines(True):
        fields = line.rstrip('\n').split('\t')
        # the copies are the rows of category X, the duplicated gene's
        if not line.startswith('#') and fields[6] == chain.COG['novel']:
            category = next(copies)
            if category is None:
                continue
            fields[6] = category
            line = '\t'.join(fields) + '\n'
        lines.append(line)
    path.write_text(''.join(lines))


@pytest.fixture(scope='module')
def communities(tmp_path_factory, tools):
    """GET_ALLDUP's homologous_communities.tsv, as text, for SAMPLES."""
    with pytest.MonkeyPatch.context() as mp:
        # Both scripts write their network with pyvis' show(), which opens it
        # in a web browser wherever there is one
        mp.setenv('BROWSER', 'true')
        return _run_duplications(tmp_path_factory.mktemp('duplications'), tools)


def _run_duplications(base, tools):
    """GET_DUPLICATIONS for each of SAMPLES, then GET_ALLDUP, in base."""
    dup_fastas = []
    for name, (categories, _) in SAMPLES.items():
        d = chain.write_sample(base / 'inputs' / name, name, copies=len(categories))
        write_eggnog(d / f'{name}.emapper.annotations', categories)
        files = {'sampleId': name, 'amap_tsv': f'{name}.AMAP.tsv', 'amap_faa': f'{name}.AMAP.faa',
                 'amap_fna': f'{name}.AMAP.fna', 'bakta_tsv': f'{name}.tsv',
                 'fasta': f'{name}.fasta', 'emapper': f'{name}.emapper.annotations'}
        task = chain.run_task(base / 'GET_DUPLICATIONS' / name, 'get.nf', 'GET_DUPLICATIONS',
                              dict(files, **{'task.cpus': 1}),
                              [d / f for key, f in files.items() if key != 'sampleId'], tools)
        dup_fastas += [p for p in task.iterdir() if p.suffix == '.dup']
    task = chain.run_task(base / 'GET_ALLDUP', 'get.nf', 'GET_ALLDUP',
                          {'task.memory': '8 GB', 'task.cpus': 1}, dup_fastas, tools)
    return (task / 'homologous_communities.tsv').read_text()


class TestGetAllDup:

    def test_columns(self, communities):
        assert communities.splitlines()[0].split('\t') == COLUMNS

    def test_elements(self, communities):
        table = pd.read_csv(StringIO(communities), sep='\t')
        assert len(table) == 1
        written = {genome: elements for k in '12'
                   for genome, elements in zip(table[f'Genome{k}'], table[f'Elements{k}'])}
        assert written == {name: elements for name, (_, elements) in SAMPLES.items()}
        # a Python dict of plain counts, as numpy 1 wrote it
        for elements in written.values():
            assert all(type(n) is int for n in ast.literal_eval(elements).values())

    def test_no_numpy_and_no_strand(self, communities):
        assert 'np.' not in communities
        assert "'+'" not in communities


# ---------------------------------------------------------------------------
# The report's duplication summary
# ---------------------------------------------------------------------------

def _chunks():
    """(header, code) of every R chunk of the report."""
    return re.findall(r'^```\{(r [^}]*)\}\n(.*?)^```', REPORT.read_text(), re.S | re.M)


class TestReportChunks:

    def test_every_chunk_that_prints_markdown_hands_it_on_as_markdown(self):
        # without results='asis', knitr puts what cat() prints in a code box,
        # every line behind "## "
        markdown = re.compile(r'''cat\(.*?["'](\*\*|#)''')
        printing = [header for header, code in _chunks() if markdown.search(code)]
        assert 'r dup_summary' in [h.split(',')[0] for h in printing]
        for header in printing:
            assert re.search(r'''results\s*=\s*["']asis["']''', header), header


def _report_chunk(label):
    """A chunk of the report, its header line and fences included."""
    return re.search(rf'^```\{{r {label}\b.*?^```\n', REPORT.read_text(), re.S | re.M).group(0)


def _dup_setup():
    """What the report's setup chunk does with the duplication table."""
    text = REPORT.read_text()
    return re.search(r'^# Duplications\n.*?^\}\n', text, re.S | re.M).group(0)


def knit_summary(work, table, rscript):
    """The markdown the report's dup_summary chunk gives on `table`, the text
    of homologous_communities.tsv."""
    (work / 'homologous_communities.tsv').write_text(table)
    (work / 'dup.Rmd').write_text(
        '```{r setup, include=FALSE}\n'
        'suppressMessages({library(ggplot2); library(ggpubr); library(dplyr)})\n'
        'params <- list(dup_file = "homologous_communities.tsv")\n'
        f'{_dup_setup()}'
        '```\n\n'
        f'{_report_chunk("dup_summary")}')
    done = subprocess.run([rscript, '-e', 'knitr::knit("dup.Rmd", "dup.md", quiet = TRUE)'],
                          cwd=work, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr[-3000:]
    return (work / 'dup.md').read_text()


class TestReportSummary:

    def test_the_summary_is_text(self, communities, rscript, tmp_path):
        md = knit_summary(tmp_path, communities, rscript)
        assert re.search(r'^\*\*1\*\* homologous community pairs detected across all genomes\.$', md, re.M), md
        assert '## **' not in md
        # the histogram of community sizes is still drawn
        assert 'figure/dup_summary-1.png' in md, md

    def test_no_duplications(self, rscript, tmp_path):
        # GET_ALLDUP's file when no sample has a duplicated gene
        md = knit_summary(tmp_path, '\n', rscript)
        assert re.search(r'^No duplication data available\.$', md, re.M), md
        assert '## No duplication' not in md
