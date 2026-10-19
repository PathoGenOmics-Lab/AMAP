#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GET_DUPLICATIONS and GET_ALLDUP write the same files from one run to the
next, and write their networks without opening them.

- homologous_communities.tsv listed its pairs of communities in the order of a
  set, which changes with Python's hash seed from one run to the next, and
  which community of a pair came first changed with it.
- AvgClustering changed in its last digit: it is an average over the
  community's proteins, and networkx went through a few proteins of a large
  network in the order of their set.
- get_duplications.py listed each group's Members in the order of a set, the
  order GET_ALLDUP is handed the sequences in, and drew its networks from a
  random layout and a random partition.
- Both scripts wrote their networks with pyvis' show(), which also opens the
  file in a web browser: on a desktop, every run opened a tab per network.
- Both loaded the networks' drawing library, vis-network, from a CDN, so
  without the internet the pages drew nothing. It now goes inside each page.

The scripts run with the arguments the two steps give them, on Diamond tables
made up here, under three hash seeds, with BROWSER set to a stand-in that
records what it is asked to open.
"""

import ast
import hashlib
import itertools
import os
import random
import re
import subprocess
import sys
from pathlib import Path

import networkx as nx
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / 'bin'
sys.path.insert(0, str(BIN))

import network_genomes  # noqa: E402

SEEDS = ['1', '2', '3']
GENOMES = ['L2_Beijing_1', 'L4_H37Rv_test', 'La1_Mbovis']
SAMPLE = 'L4_H37Rv_test'   # the genome get_duplications.py runs on

# A family of seven whose clustering coefficients are 2/3, 0.7 and 0.8: their
# average, a sum of floats, comes out 0.7238095238095237, ...238 or ...239
# depending on the order of its nodes
FAMILY = [e for e in itertools.combinations('abcdefg', 2)
          if e not in {('a', 'b'), ('a', 'c'), ('d', 'e'), ('f', 'g')}]
TRIANGLE = list(itertools.combinations('xyz', 2))
# In each genome: the family and five families of three, so that the family
# holds fewer than half the nodes of its genome's network, which is when
# networkx goes through a subgraph's nodes in the order of their set
FAMILIES = {'F': FAMILY, **{f'T{k}': TRIANGLE for k in range(5)}}

OUTPUTS = ['homologous_communities.tsv', 'homologous_communities.html',
           f'duplications_{SAMPLE}_groups.tsv', f'duplications_{SAMPLE}_network.png',
           f'duplications_{SAMPLE}_network_interactive.html']

RECORDER = '#!/bin/sh\necho "$@" >> "{log}"\n'


def hits(pairs):
    """Diamond's table as the two steps write it: its header, then each pair
    of identical proteins both ways, sorted."""
    rows = sorted('\t'.join(map(str, [q, s, 100.0, 100, 0, 0, 1, 100, 1, 100, 1e-50, 200.0, 100, 100, 100]))
                  for a, b in pairs for q, s in ((a, b), (b, a)))
    return ''.join(f'{line}\n' for line in ['# DIAMOND'] + rows)


def write_inputs(d):
    """GET_DUPLICATIONS' Diamond table for SAMPLE, and GET_ALLDUP's for all
    of GENOMES with its genome table: within each genome every family, and
    between two genomes each family to the same family of the other."""
    within, genome_of = {}, {}
    for genome in GENOMES:
        within[genome] = [(f'{genome}_{fam}{u}', f'{genome}_{fam}{v}')
                          for fam, links in FAMILIES.items() for u, v in links]
        genome_of.update({n: genome for pair in within[genome] for n in pair})
    between = [(f'{g1}_{fam}{links[0][0]}', f'{g2}_{fam}{links[0][0]}')
               for g1, g2 in itertools.combinations(GENOMES, 2) for fam, links in FAMILIES.items()]
    (d / f'duplications_{SAMPLE}').write_text(hits(within[SAMPLE]))
    (d / 'all_vs_all.diamond').write_text(hits([p for g in GENOMES for p in within[g]] + between))
    (d / 'genomes.tsv').write_text(''.join(f'{n}\t{g}\n' for n, g in sorted(genome_of.items())))


def run_scripts(d, seed, browser):
    """get_duplications.py and network_genomes.py in d, as the two steps run
    them, under a hash seed, with `browser` as the web browser."""
    env = dict(os.environ, PYTHONHASHSEED=seed, BROWSER=str(browser),
               PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
    for script, args in [('get_duplications.py', [f'duplications_{SAMPLE}']),
                         ('network_genomes.py', ['all_vs_all.diamond', '--genomes', 'genomes.tsv'])]:
        done = subprocess.run([sys.executable, str(BIN / script)] + args, cwd=d, env=env,
                              capture_output=True, text=True)
        assert done.returncode == 0, f'{script} failed:\n{done.stdout[-1500:]}\n{done.stderr[-3000:]}'


def stand_in_browser(d):
    """A web browser in d that only writes what it is asked to open into
    d/opened.txt, and that file."""
    browser, log = d / 'browser', d / 'opened.txt'
    browser.write_text(RECORDER.format(log=log))
    browser.chmod(0o755)
    return browser, log


@pytest.fixture(scope='module')
def browser(tmp_path_factory):
    """The web browser of the runs, and the file it writes into."""
    return stand_in_browser(tmp_path_factory.mktemp('browser'))


@pytest.fixture(scope='module')
def runs(tmp_path_factory, browser):
    """The folder each hash seed's run wrote its files in."""
    out = {}
    for seed in SEEDS:
        d = out[seed] = tmp_path_factory.mktemp(f'seed{seed}')
        write_inputs(d)
        run_scripts(d, seed, browser[0])
    return out


# ---------------------------------------------------------------------------
# The same bytes under every hash seed
# ---------------------------------------------------------------------------

class TestSameFiles:

    @pytest.mark.parametrize('name', OUTPUTS)
    def test_every_hash_seed_writes_the_same_file(self, runs, name):
        written = {seed: hashlib.md5((d / name).read_bytes()).hexdigest()[:8] for seed, d in runs.items()}
        assert len(set(written.values())) == 1, f'{name}, md5 by hash seed: {written}'

    def test_pairs_in_the_order_of_the_communities(self, runs):
        table = pd.read_csv(runs['1'] / 'homologous_communities.tsv', sep='\t')
        number = lambda name: int(name.rsplit('_', 1)[1])  # noqa: E731
        pairs = [((g1, number(c1)), (g2, number(c2))) for g1, c1, g2, c2 in
                 zip(table['Genome1'], table['Community1'], table['Genome2'], table['Community2'])]
        assert pairs == sorted(pairs)
        # each pair once, the first of its two communities on the left, so
        # L4_H37Rv_test comes second next to L2_Beijing_1 and first next to
        # La1_Mbovis, where the order of the set put it either way
        assert all(one < two for one, two in pairs)
        assert len(set(pairs)) == len(pairs) == len(FAMILIES) * 3   # three pairs of genomes

    def test_members_in_alphabetical_order(self, runs):
        groups = pd.read_csv(runs['1'] / f'duplications_{SAMPLE}_groups.tsv', sep='\t')
        members = [ast.literal_eval(m) for m in groups['Members']]
        assert all(m == sorted(m) for m in members)
        assert sorted(map(len, members)) == [3, 3, 3, 3, 3, 7]


class TestMetrics:

    def test_the_order_of_the_nodes_changes_nothing(self):
        # a subgraph view of a few nodes of a large network goes through them
        # in the order of their set, which changes from one run to the next
        nodes, rng, written = list('abcdefg'), random.Random(0), set()
        for _ in range(50):
            rng.shuffle(nodes)
            graph = nx.Graph()
            graph.add_nodes_from(nodes)
            graph.add_edges_from(FAMILY)
            written.add(tuple(network_genomes.calculate_subgraph_metrics(graph).items()))
        assert len(written) == 1, written


# ---------------------------------------------------------------------------
# No web browser
# ---------------------------------------------------------------------------

class TestNoBrowser:

    def test_no_network_is_opened(self, runs, browser):
        log = browser[1]
        assert not log.exists(), f'opened in a web browser: {log.read_text()}'
        for d in runs.values():
            assert (d / 'homologous_communities.html').stat().st_size
            assert (d / f'duplications_{SAMPLE}_network_interactive.html').stat().st_size

    def test_nor_an_empty_one(self, tmp_path):
        # without hits, both scripts still write a network
        (tmp_path / f'duplications_{SAMPLE}').write_text('# DIAMOND\n')
        (tmp_path / 'all_vs_all.diamond').write_text('# DIAMOND\n')
        (tmp_path / 'genomes.tsv').write_text('')
        stand_in, log = stand_in_browser(tmp_path)
        run_scripts(tmp_path, '1', stand_in)
        assert not log.exists(), f'opened in a web browser: {log.read_text()}'
        assert (tmp_path / 'homologous_communities.html').stat().st_size
        assert (tmp_path / f'duplications_{SAMPLE}_network_interactive.html').stat().st_size
        assert _offline(tmp_path / 'homologous_communities.html')
        assert _offline(tmp_path / f'duplications_{SAMPLE}_network_interactive.html')


# ---------------------------------------------------------------------------
# Networks that draw without the internet
# ---------------------------------------------------------------------------

def _offline(page):
    """Whether a network page carries its drawing library instead of asking a
    CDN for it."""
    html = page.read_text(encoding='utf-8')
    return '* vis-network' in html and not re.search(r'(?:src|href)="https?://[^"]*vis-network', html)


class TestOffline:

    def test_the_networks_draw_without_the_internet(self, runs):
        for d in runs.values():
            for name in ('homologous_communities.html', f'duplications_{SAMPLE}_network_interactive.html'):
                assert _offline(d / name), name
