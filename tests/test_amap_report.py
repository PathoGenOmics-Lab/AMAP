#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bin/report/amap_report.py writes the interactive report, amap_report.html.

Most tests run it on tests/data/report/, four genomes of a real AMAP run cut
to two stretches of the reference (tests/data/report/make_fixture.py says how
and from what): H37Rv, CDC1551 in 40 contigs, M. bovis AF2122/97 and CCDC5079.
Their genes are in every state the report draws: M. bovis lacks RD4, a contig
end cuts dnaE2 in the draft, CCDC5079 has two copies of Rv3424c, and the IS6110
copies name seventeen reference genes at once.

The page carries its data compressed, and its JavaScript decodes the genes.
The decoders are checked against the files they were made from twice: in
Python (decode_genes, decode_gone) and, where node is installed, with the
page's own JavaScript.
"""

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
REPORT = ROOT / 'bin' / 'report'
DATA = ROOT / 'tests' / 'data' / 'report'
sys.path.insert(0, str(REPORT))

import amap_report as ar  # noqa: E402

GENOMES = ['L4_H37Rv', 'L4_CDC1551draft', 'A4_AF2122', 'L2_CCDC5079']


def _args(tmp_path, *extra, **kw):
    argv = ['--amap', str(DATA / 'amap'), '--ref-proteome', str(DATA / 'reference.faa'),
            '--out', str(tmp_path / 'amap_report.html')]
    if kw.get('full', True):
        argv += ['--cases', str(DATA / 'cases'), '--bakta', str(DATA / 'bakta'),
                 '--mycobrowser', str(DATA / 'mycobrowser.txt'),
                 '--versions', str(DATA / 'software_versions.tsv'),
                 '--params', str(DATA / 'params.yml'), '--amap-version', '1.0.1']
    return argv + list(extra)


@pytest.fixture(scope='module')
def run(tmp_path_factory):
    """The report of the fixture, and the data in it, read back."""
    tmp = tmp_path_factory.mktemp('report')
    ar.main(_args(tmp))
    page = (tmp / 'amap_report.html').read_text(encoding='utf-8')
    return page, ar.load_embedded(page)


def _states(data, name):
    """{reference gene: (state, detail)} of one genome, from the page."""
    out, idx = {}, 0
    for tok in filter(None, data['status'][name].split(',')):
        m = re.match(r'^(\d+)([A-Z])(.*)$', tok)
        idx += int(m.group(1))
        out[data['ref']['names'][idx]] = (m.group(2), m.group(3))
    return out


def _sample(data, name):
    return next(s for s in data['samples'] if s['name'] == name)


def _ref(data):
    return {'names': data['ref']['names'],
            'start': [x or 0 for x in data['ref']['start']],
            'end': [x or 0 for x in data['ref']['end']]}


# --------------------------------------------------------------------------
# The page
# --------------------------------------------------------------------------

class TestPage:
    def test_one_file_that_needs_no_network(self, run):
        page, _ = run
        # Scripts, styles, images and fonts: none from anywhere else
        assert not re.search(r'<script[^>]+src=', page)
        assert not re.search(r'<link[^>]+href="https?:', page)
        assert not re.search(r'<img[^>]+src="https?:', page)
        assert not re.search(r'url\(\s*["\']?https?:', page)
        assert '@import' not in page
        assert 'fetch(' not in page
        # The only address in the page is a link a reader can follow
        for m in re.finditer(r'(?:src|href)="(https?:[^"]*)"', page):
            start = page.rfind('<', 0, m.start())
            assert page[start:start + 2] == '<a', m.group(0)

    def test_the_logo_and_the_data_are_inside(self, run):
        page, data = run
        assert 'src="data:image/png;base64,' in page
        assert '{{' not in page.replace('{{TITLE}}', '')   # every marker filled in
        assert [s['name'] for s in data['samples']] == ['L2_CCDC5079', 'L4_CDC1551draft', 'L4_H37Rv', 'A4_AF2122']

    def test_small(self, run):
        page, _ = run
        assert len(page.encode('utf-8')) < 400_000

    def test_no_em_dash(self, run):
        page, _ = run
        assert '\u2014' not in page
        for path in REPORT.glob('amap_report.*'):
            assert '\u2014' not in path.read_text(encoding='utf-8'), path.name

    def test_text_cannot_end_the_data_script(self, tmp_path):
        # A product with "</script>" in it, as plain JSON (--plain)
        bakta = tmp_path / 'bakta'
        shutil.copytree(DATA / 'bakta', bakta)
        tsv = bakta / 'L2_CCDC5079' / 'L2_CCDC5079.tsv'
        lines = tsv.read_text().splitlines()
        noref = {r['id'] for r in ar._rows(DATA / 'amap' / 'L2_CCDC5079.AMAP.tsv') if not r['rvs']}
        out = []
        for line in lines:
            f = line.split('\t')
            if not line.startswith('#') and len(f) > 7 and f[5] in noref:
                f[7] = 'evil </script><script>alert(1)</script> protein'
            out.append('\t'.join(f))
        tsv.write_text('\n'.join(out) + '\n')
        argv = _args(tmp_path, '--plain')
        argv[argv.index('--bakta') + 1] = str(bakta)
        ar.main(argv)
        page = (tmp_path / 'amap_report.html').read_text()
        # The HTML parser ends a script at the first </script: only the
        # page's own three may be there
        assert page.count('</script') == 3
        data = ar.load_embedded(page)
        assert any('</script>' in p for p in data['products'])

    def test_title_is_escaped(self, tmp_path):
        ar.main(_args(tmp_path, '--title', 'Run <1> & "2"'))
        page = (tmp_path / 'amap_report.html').read_text()
        assert '<title>Run &lt;1&gt; &amp; &quot;2&quot;</title>' in page


# --------------------------------------------------------------------------
# The genes, encoded and decoded
# --------------------------------------------------------------------------

class TestGenes:
    @pytest.mark.parametrize('name', GENOMES)
    def test_every_gene_comes_back_as_written(self, run, name):
        _, data = run
        enc = data['tables'][name]
        strand = [c == '+' for c in data['ref']['strand']]
        got = ar.decode_genes(name, enc, _ref(data), strand, enc['cn'])
        want = ar.read_amap_table(DATA / 'amap' / f'{name}.AMAP.tsv')
        assert len(got) == len(want) > 100
        for g, w in zip(got, want):
            assert (g['id'], g['contig'], g['start'], g['end'], g['strand'], g['pseudo'], g['partial']) == \
                   (w['id'], w['contig'], w['start'], w['end'], w['strand'], w['pseudo'], w['partial'])
            assert g['rvs'] == w['rvs']

    @pytest.mark.parametrize('name', GENOMES)
    def test_every_call_left_out_comes_back(self, run, name):
        _, data = run
        enc = data['tables'][name]
        strand = [c == '+' for c in data['ref']['strand']]
        genes = ar.decode_genes(name, enc, _ref(data), strand, enc['cn'])
        got = {c['tag']: (c['case'], c['start'], c['end'], c['strand'], c['contig'])
               for c in ar.decode_gone(name, enc['gone'], genes)}
        cases = {c: ar.read_case_table(DATA / 'cases' / f'{name}-case{c}.tsv') for c in '1234'}
        case_of = ar.case_of_calls(cases)
        final = {g['id'] for g in genes}
        want = {c['tag']: (case, c['start'], c['end'], c['strand'], c['contig'])
                for case, calls in cases.items() for c in calls
                if c['tag'] not in final and case_of[c['tag']] == case}
        assert got == want
        assert got, 'every genome here has calls the correction replaced'

    def test_a_reversed_contig_is_cheap(self, run):
        # The draft's contigs run both ways along the reference; their genes
        # still cost little more than those of a closed genome
        _, data = run
        size = {n: len(json.dumps({k: v for k, v in data['tables'][n].items() if k != 'gone'})) for n in GENOMES}
        assert size['L4_CDC1551draft'] < 3 * size['L4_H37Rv']

    def test_rle(self):
        values = ['1', '1', '1', '=4', 'x', 'x', '-1']
        assert ar.rle(values) == '1*3,=4,x*2,-1'
        assert ar.unrle(ar.rle(values)) == values
        assert ar.unrle('') == []

    @pytest.mark.skipif(not shutil.which('node'), reason='needs node')
    @pytest.mark.parametrize('name', GENOMES)
    def test_the_page_decodes_as_python_does(self, run, name, tmp_path):
        _, data = run
        enc = data['tables'][name]
        inp = tmp_path / 'in.json'
        inp.write_text(json.dumps({'name': name, 'enc': enc, 'ref': _ref(data),
                                   'plus': [c == '+' for c in data['ref']['strand']]}))
        script = (
            "const m = require(process.argv[1]);"
            "const d = JSON.parse(require('fs').readFileSync(process.argv[2], 'utf8'));"
            "const g = m.decodeGenes(d.name, d.enc, d.ref, d.plus);"
            "const gone = m.decodeGone(d.name, d.enc.gone, g);"
            "process.stdout.write(JSON.stringify({ids: g.ids, start: Array.from(g.start), end: Array.from(g.end),"
            " plus: Array.from(g.plus), contig: Array.from(g.contig).map((c) => g.contigs[c]),"
            " gone: gone.map((x) => [x.tag, x.case, x.start, x.end, x.strand, x.contig])}));")
        out = subprocess.run(['node', '-e', script, str(REPORT / 'amap_report.js'), str(inp)],
                             capture_output=True, text=True, check=True).stdout
        js = json.loads(out)
        py = ar.decode_genes(name, enc, _ref(data), [c == '+' for c in data['ref']['strand']], enc['cn'])
        assert js['ids'] == [g['id'] for g in py]
        assert js['start'] == [g['start'] for g in py]
        assert js['end'] == [g['end'] for g in py]
        assert js['plus'] == [1 if g['strand'] == '+' else 0 for g in py]
        assert js['contig'] == [g['contig'] for g in py]
        assert js['gone'] == [[c['tag'], c['case'], c['start'], c['end'], c['strand'], c['contig']]
                              for c in ar.decode_gone(name, enc['gone'], py)]


# --------------------------------------------------------------------------
# What the report says of the genomes
# --------------------------------------------------------------------------

class TestStates:
    def test_rd4_is_missing_from_m_bovis(self, run):
        _, data = run
        states = _states(data, 'A4_AF2122')
        rd4 = ['Rv1507c', 'Rv1507A', 'Rv1508c', 'Rv1508A', 'Rv1509', 'Rv1510', 'Rv1511',
               'Rv1512', 'Rv1513', 'Rv1514c', 'Rv1515c', 'Rv1516c']
        assert all(states[g][0] == 'M' for g in rd4)
        assert 'Rv1507c' not in _states(data, 'L4_H37Rv')

    def test_the_deletion_is_a_block(self, run):
        _, data = run
        names = data['ref']['names']
        blocks = {(names[b['first']], names[b['last']]): sorted({r[0] for r in b['runs']}) for b in data['blocks']}
        # RD4 in M. bovis, and PPE57 to PPE59, which CDC1551 lacks too; a
        # single missing gene is not a block
        assert blocks == {('Rv1507c', 'Rv1516c'): ['A4_AF2122'],
                          ('Rv3425', 'Rv3429'): ['A4_AF2122', 'L4_CDC1551draft']}

    def test_one_missing_gene_is_not_a_block(self):
        ref = {'start': [100, 200, 300, 400, 500, None]}
        status = {'a': {1: ('M', ''), 3: ('M', ''), 4: ('M', '')},
                  'b': {0: ('M', ''), 2: ('P', '0'), 4: ('M', ''), 5: ('M', '')}}
        blocks = ar.deletion_blocks(ref, status)
        # Only genes 3 and 4 of genome a run on; gene 5 has no place
        assert blocks == [{'first': 3, 'last': 4, 'runs': [['a', 3, 4]]}]

    def test_a_gene_a_contig_end_cuts_is_partial(self, run):
        _, data = run
        assert _states(data, 'L4_CDC1551draft')['Rv3370c'] == ('T', '3')
        assert _sample(data, 'L4_CDC1551draft')['partial'] == 1

    def test_a_pseudogene_keeps_its_reason(self, run):
        _, data = run
        code, note = _states(data, 'A4_AF2122')['Rv3373']
        assert code == 'P'
        assert data['notes'][int(note)] == ('frameshift or lost stop codon: its frame runs on 83 '
                                            'residues past the end of the reference protein')

    def test_copies_shared_copies_and_length(self, run):
        _, data = run
        assert _states(data, 'L2_CCDC5079')['Rv3424c'] == ('C', '2')
        # An IS6110 copy names seventeen reference genes at once
        assert _states(data, 'L2_CCDC5079')['Rv3474'][0] == 'S'
        assert _states(data, 'A4_AF2122')['Rv1506c'] == ('L', '-426')

    def test_the_reference_strain_is_intact(self, run):
        _, data = run
        assert _states(data, 'L4_H37Rv') == {}
        h37rv = _sample(data, 'L4_H37Rv')
        assert h37rv['found'] == len(data['ref']['names'])
        assert h37rv['found_pct'] == 100.0

    def test_counts_agree_with_the_states(self, run):
        _, data = run
        for name in GENOMES:
            s, states = _sample(data, name), _states(data, name)
            assert s['missing'] == sum(1 for c, _ in states.values() if c == 'M')
            assert s['copies'] == sum(1 for c, _ in states.values() if c == 'C')
            assert s['found'] == len(data['ref']['names']) - s['missing']

    def test_sample_ids_and_lineages_from_the_summaries(self, run):
        _, data = run
        assert {(s['name'], s['id'], s['lineage']) for s in data['samples']} == {
            ('L4_H37Rv', 'H37Rv', 'L4'), ('L4_CDC1551draft', 'CDC1551draft', 'L4'),
            ('A4_AF2122', 'AF2122', 'A4'), ('L2_CCDC5079', 'CCDC5079', 'L2')}

    def test_genome_statistics_from_bakta(self, run):
        _, data = run
        draft = _sample(data, 'L4_CDC1551draft')
        assert (draft['size'], draft['contigs'], draft['n50'], draft['gc']) == (4403837, 40, 231675, 65.6)
        enc = data['tables']['L4_CDC1551draft']
        assert len(enc['cn']) == len(enc['cl']) == 40 and sum(enc['cl']) == 4403837


class TestCorrection:
    @pytest.mark.parametrize('name', GENOMES)
    def test_the_steps_add_up(self, run, name):
        _, data = run
        s = _sample(data, name)
        f = s['flow']
        calls = f['c1'] + f['c2'] + f['c3'] + f['c4']
        left_out = (f['c1'] - f['k1']) + (f['c2'] - f['k2']) + (f['c3'] - f['k3']) + (f['c4'] - f['k4'])
        assert calls - left_out + f['R'] + f['A'] + f['other'] == s['genes']
        assert len(ar.decode_gone(name, data['tables'][name]['gone'],
                                  ar.decode_genes(name, data['tables'][name], _ref(data),
                                                  [c == '+' for c in data['ref']['strand']],
                                                  data['tables'][name]['cn']))) == left_out

    def test_a_call_in_two_cases_counts_once(self, run):
        # CDC1551's Bakta call 03643 matches one placement exactly and another
        # in part, so TYPE_CASES writes it in both case 2 and case 3
        _, data = run
        f = _sample(data, 'L4_CDC1551draft')['flow']
        assert f['twice'] == 1
        assert f['c1'] + f['c2'] + f['c3'] + f['c4'] == _sample(data, 'L4_CDC1551draft')['cases']['total_cds']

    def test_a_call_dropped_after_its_case_is_listed(self):
        # A case 2 call that the last deduplication or the stop codon check
        # took out is a call the final table does not keep, like the others
        genes = [{'id': 'S_00001', 'contig': 'c', 'start': 100, 'end': 399, 'strand': '+', 'rvs': ['Rv1']},
                 {'id': 'S_Rv2', 'contig': 'c', 'start': 500, 'end': 799, 'strand': '+', 'rvs': ['Rv2']}]
        cases = {'2': [{'tag': 'S_00001', 'start': 100, 'end': 399, 'strand': '+', 'contig': 'c'},
                       {'tag': 'S_00002', 'start': 900, 'end': 1199, 'strand': '+', 'contig': 'c'}],
                 '4': [{'tag': 'S_00003', 'start': 530, 'end': 799, 'strand': '+', 'contig': 'c'}]}
        case_of = ar.case_of_calls(cases)
        enc = ar.encode_gone('S', genes, cases, case_of, {'S_00001', 'S_Rv2'})
        gone = ar.decode_gone('S', enc, genes)
        assert [(c['tag'], c['case'], c['row'], c['start'], c['end']) for c in gone] == \
            [('S_00002', '2', -1, 900, 1199), ('S_00003', '4', 1, 530, 799)]

    def test_case_of_calls_prefers_the_better_case(self):
        cases = {'3': [{'tag': 'a'}], '2': [{'tag': 'a'}, {'tag': 'b'}], '4': [{'tag': 'c'}], '1': [{'tag': 'b'}]}
        assert ar.case_of_calls(cases) == {'a': '2', 'b': '2', 'c': '4'}

    def test_origins(self):
        genes = [{'id': 'S_00001', 'rvs': ['Rv1']}, {'id': 'S_Rv2', 'rvs': ['Rv2']},
                 {'id': 'S_Rv3', 'rvs': ['Rv3']}, {'id': 'S_00009', 'rvs': []}]
        assert ar.origins_of('S', genes, {'S_00001': '2', 'S_00009': '1'}, {'Rv3'}) == ['2', 'R', 'A', '1']
        # Without case tables only the kind of identifier is known
        assert ar.origins_of('S', genes, {}, set()) == ['b', 'r', 'r', 'b']


class TestFlags:
    @staticmethod
    def _genomes(values, key='missing', lineage='L4'):
        return [{'name': f'g{i}', 'lineage': lineage, key: v} for i, v in enumerate(values)]

    def test_far_from_the_others(self):
        genomes = self._genomes([10, 12, 9, 11, 10, 13, 80])
        ar.robust_flags(genomes)
        flagged = [g['name'] for g in genomes if g['flags']]
        assert flagged == ['g6']
        f = genomes[-1]['flags'][0]
        assert f['key'] == 'missing' and f['value'] == 80 and f['median'] == 10.5 and f['z'] > 3.5

    def test_only_on_the_side_that_is_bad(self):
        genomes = self._genomes([99.1, 99.2, 99.0, 99.3, 99.1, 99.9], key='found_pct')
        ar.robust_flags(genomes)
        assert not any(g['flags'] for g in genomes)
        genomes = self._genomes([99.1, 99.2, 99.0, 99.3, 99.1, 95.0], key='found_pct')
        ar.robust_flags(genomes)
        assert [g['name'] for g in genomes if g['flags']] == ['g5']

    def test_identical_genomes_flag_a_real_gap_only(self):
        # Every other genome has one contig: the deviation is 0, so any gap
        # would be infinitely far; only one past the floor is flagged
        genomes = self._genomes([1, 1, 1, 1, 2], key='contigs') + self._genomes([40], key='contigs')
        for i, g in enumerate(genomes):
            g['name'] = f'g{i}'
        ar.robust_flags(genomes)
        assert [g['name'] for g in genomes if g['flags']] == ['g5']
        assert genomes[5]['flags'][0]['z'] is None

    def test_three_genomes_are_too_few(self):
        genomes = self._genomes([10, 10, 200])
        ar.robust_flags(genomes)
        assert not any(g['flags'] for g in genomes)

    def test_a_large_lineage_is_its_own_group(self):
        # M. bovis lacks its RDs: among six A4 genomes 78 missing genes is normal
        genomes = self._genomes([10, 11, 9, 12, 10, 11], lineage='L4') + \
            self._genomes([78, 80, 77, 79, 81, 78], lineage='A4')
        for i, g in enumerate(genomes):
            g['name'] = f'g{i}'
        ar.robust_flags(genomes)
        assert not any(g['flags'] for g in genomes)


# --------------------------------------------------------------------------
# Inputs that are not there
# --------------------------------------------------------------------------

class TestWhatIsMissing:
    def test_the_annotation_alone_is_enough(self, tmp_path):
        ar.main(_args(tmp_path, full=False))
        data = ar.load_embedded((tmp_path / 'amap_report.html').read_text())
        modes = data['meta']['modes']
        assert not modes['cases'] and not modes['genomes'] and not modes['mycobrowser'] and not modes['functional']
        assert all(s['lineage'] == 'NA' and 'flow' not in s for s in data['samples'])
        assert all('gone' not in t for t in data['tables'].values())
        # The identifiers still come back, from the kind of identifier alone
        name = 'L4_CDC1551draft'
        genes = ar.decode_genes(name, data['tables'][name], _ref(data),
                                [c == '+' for c in data['ref']['strand']], data['tables'][name]['cn'])
        assert [g['id'] for g in genes] == [r['id'] for r in ar.read_amap_table(DATA / 'amap' / f'{name}.AMAP.tsv')]

    def test_a_genome_without_an_annotation_is_named(self, tmp_path):
        cases = tmp_path / 'cases'
        shutil.copytree(DATA / 'cases', cases)
        summary = (cases / 'L4_H37Rv-summary.tsv').read_text()
        (cases / 'L1_G1-summary.tsv').write_text(summary.replace('L4_H37Rv', 'L1_G1').replace('\tL4\t', '\tL1\t').replace('H37Rv', 'G1'))
        expected = tmp_path / 'samples.txt'
        expected.write_text('H37Rv\nCDC1551draft\nAF2122\nCCDC5079\nG1\nG2\n')
        argv = _args(tmp_path, '--expected', str(expected))
        argv[argv.index('--cases') + 1] = str(cases)
        ar.main(argv)
        data = ar.load_embedded((tmp_path / 'amap_report.html').read_text())
        failed = [s for s in data['samples'] if s.get('failed')]
        assert [(s['name'], s['id'], s['lineage']) for s in failed] == [('L1_G1', 'G1', 'L1')]
        assert data['meta']['expected_missing'] == ['G2']

    def test_no_genome_at_all(self, tmp_path):
        empty = tmp_path / 'empty'
        empty.mkdir()
        argv = _args(tmp_path, full=False)
        argv[argv.index('--amap') + 1] = str(empty)
        ar.main(argv)
        data = ar.load_embedded((tmp_path / 'amap_report.html').read_text())
        assert data['samples'] == [] and data['tables'] == {}

    def test_mycobrowser_names_the_genes(self, run, tmp_path):
        _, data = run
        i = data['ref']['names'].index('Rv3457c')
        assert (data['ref']['name'][i], data['ref']['cat'][i]) == ('rpoA', 'information pathways')
        assert data['ref']['src'][i] == 'm'
        # Without it, the names and products come from the genomes' own annotation
        argv = _args(tmp_path)
        del argv[argv.index('--mycobrowser'):argv.index('--mycobrowser') + 2]
        ar.main(argv)
        bare = ar.load_embedded((tmp_path / 'amap_report.html').read_text())
        assert not bare['meta']['modes']['mycobrowser']
        assert bare['ref']['cat'][i] == '' and bare['ref']['refpseudo'][i] == -1
        assert bare['ref']['product'][i]           # Bakta's product
        assert bare['ref']['src'][i] == 'e'

    def test_params_as_dump_params_writes_them(self):
        # Only the first parameter is indented under `params:` in that file
        params, info = ar.read_params(DATA / 'params.yml')
        assert len(params) == 28
        assert params['skip_functional_annotation'] == 'true'
        assert params['mycobrowser_ref'] == '/reference/Mycobacterium_tuberculosis_H37Rv_txt_v4.txt'
        assert info['nextflow'] == '26.10.2' and info['profile'] == 'docker'

    def test_the_methods_carry_the_run(self, run):
        _, data = run
        meta = data['meta']
        assert meta['version'] == '1.0.1'
        assert dict(meta['versions'])['Bakta'] == '1.12.0'
        assert meta['modes']['skip_functional'] and not meta['reference']['custom']


class TestFunctionalInputs:
    """No eggNOG database ran on the fixture's genomes, so these tables are
    written here, in the shape GET_OUTPUTS, PANGENE and GET_ALLDUP give them."""

    def test_master_tables_give_cog_and_kegg_counts(self, tmp_path):
        master = tmp_path / 'master'
        master.mkdir()
        rows = list(ar._rows(DATA / 'amap' / 'L4_H37Rv.AMAP.tsv'))
        cogs = ['J', 'KT', '-', 'S']
        with open(master / 'L4_H37Rv-AMAPmastertable.tsv', 'w') as fh:
            fh.write('id\tCOG_category\tKEGG_Pathway\tPreferred_name\tDescription\n')
            for k, r in enumerate(rows):
                fh.write(f"{r['id']}\t{cogs[k % 4]}\t{'ko00010,map00010' if k % 2 else '-'}\tgene{k}\tdesc {k}\n")
        argv = _args(tmp_path, '--master', str(master))
        ar.main(argv)
        data = ar.load_embedded((tmp_path / 'amap_report.html').read_text())
        assert data['meta']['modes']['functional']
        cog = data['functional']['cog']['L4_H37Rv']
        n = len(rows)
        assert cog['J'] == len(range(0, n, 4))
        assert cog['K'] == cog['T'] == len(range(1, n, 4))     # a gene with two letters counts in each
        assert cog['-'] == len(range(2, n, 4))
        assert data['functional']['kegg']['L4_H37Rv'] == {'map00010': n // 2}
        assert set(data['functional']['cog']) == {'L4_H37Rv'}

    def test_pangene_matrix(self, tmp_path):
        # Two genomes' Rtabs merged as REPORT_INPUTS merges them: each genome's
        # genes against the other genome's sequence, named <other>_<genome>#0
        rtab = tmp_path / 'input_pangene.txt'
        rtab.write_text(
            'Gene\tA4_AF2122_L4_H37Rv#0\tL4_H37Rv_A4_AF2122#0\n'
            'L4_H37Rv_Rv1507c;+;0-100\t0\t0\n'
            'L4_H37Rv_Rv0001;+;0-100\t1\t0\n'
            'A4_AF2122_Rv0001;+;0-100\t0\t1\n')
        ar.main(_args(tmp_path, '--pangene', str(rtab)))
        data = ar.load_embedded((tmp_path / 'amap_report.html').read_text())
        assert data['pangene'] == {'L4_H37Rv': {'others': 1, 'hist': {'0': 1, '1': 1}},
                                   'A4_AF2122': {'others': 1, 'hist': {'1': 1}}}

    def test_homologous_communities(self, tmp_path):
        tsv = tmp_path / 'homologous_communities.tsv'
        tsv.write_text('Community\tGenomes\tNodes\n1\t2\tL4_H37Rv_00001,A4_AF2122_00002\n')
        ar.main(_args(tmp_path, '--communities', str(tsv)))
        data = ar.load_embedded((tmp_path / 'amap_report.html').read_text())
        assert data['communities'] == [{'Community': '1', 'Genomes': '2', 'Nodes': 'L4_H37Rv_00001,A4_AF2122_00002'}]


class TestReference:
    def test_names_with_underscores_and_no_coordinates(self, tmp_path):
        faa = tmp_path / 'ref.faa'
        faa.write_text('>PE_PGRS1_100-400\nMAAA*\n>b0001_190-255\nMKR\n>orphan\nMA\n')
        genes, warnings = ar.read_proteome(faa)
        assert [(g['name'], g['start'], g['end'], g['aa']) for g in genes] == \
            [('PE_PGRS1', 100, 400, 4), ('b0001', 190, 255, 3), ('orphan', None, None, 2)]
        assert warnings == []

    def test_a_name_given_twice(self, tmp_path):
        faa = tmp_path / 'ref.faa'
        faa.write_text('>g1_1-30\nMA\n>g1_50-80\nMA\n')
        genes, warnings = ar.read_proteome(faa)
        assert len(genes) == 1 and 'g1' in warnings[0]


class TestPipeline:
    NF = ROOT / 'modules' / 'local' / 'amap_report.nf'

    def test_every_option_the_process_passes_exists(self):
        script = self.NF.read_text()
        call = script[script.index('amap_report.py'):script.index('--out amap_report.html')]
        options = set(re.findall(r'(--[a-z][\w-]*)', call)) | {'--pangene', '--communities'}
        parser_options = set(re.findall(r"add_argument\('(?:-\w, )?'?(--[\w-]+)'", (REPORT / 'amap_report.py').read_text()))
        assert options <= parser_options, options - parser_options

    def test_the_folders_it_stages_are_the_ones_it_reads(self):
        script = self.NF.read_text()
        staged = set(re.findall(r"stageAs: '([\w-]+)/\*'", script))
        for folder in ('amap', 'cases', 'bakta', 'master'):
            assert folder in staged and f'--{folder} {folder}' in script
        assert {'run', 'reference', 'pangene', 'duplications', 'mycobrowser'} <= staged
        for name in ('software_versions.tsv', 'params.yml', 'samples.txt'):
            assert f'run/{name}' in script

    def test_main_hands_it_what_it_declares(self):
        main = (ROOT / 'main.nf').read_text()
        call = main[main.index('AMAP_REPORT('):]
        call = call[:call.index('\n        )') + 10]
        inputs = re.findall(r'^    path (\w+)', self.NF.read_text(), re.M)
        args = [a for a in re.split(r',\n', call[call.index('(') + 1:call.rindex(')')].strip()) if a.strip()]
        assert len(args) == len(inputs) == 9
