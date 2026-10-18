#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A reference file given on the host reaches the task that reads it.

BAKTA_ANNOTATOR passed --ref_annotation, or --bakta_proteins, to Bakta as a
path written into its script, and PATHOTYPR, GET_OUTPUTS and BUNDLED_REFERENCE
did the same with the lineage markers, Mycobrowser's table and the miniprot
references. A container sees only what Nextflow mounts for a task, and a task
runs in a directory of its own, so a file on the host, under a relative path
or in a folder nothing else mounts, was not there for it. With K-12 as the
custom reference for Sakai (2026-10-18), Bakta said "user proteins file
(../inputs/ecoli/k12.gbff) not valid!" three times, the sample was ignored and
the run ended with exit 0; given the lineage markers on the host, pathotypr
found no file and every genome went on untyped, as NA_<sampleId>.

A file on the host is now staged into the task as an input, in a directory of
its own, and a default, which lives inside the container, is read by its path.

The module tests render each script the way Nextflow does, its own Groovy
included, and run it in a task directory laid out the way Nextflow lays one
out, with stand-ins for Bakta, pathotypr and the Python it calls. As in
tests/test_validate_inputs.py, a folder the task cannot see, as a container
sees no folder Nextflow does not mount for it, is one that does not exist.
"""

import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'bin'))

SCRIPT = ROOT / 'bin' / 'validate_inputs.py'
NF_FILES = sorted((ROOT / 'modules' / 'local').glob('*.nf')) + [ROOT / 'main.nf']


# ---------------------------------------------------------------------------
# A process, rendered as Nextflow renders it
# ---------------------------------------------------------------------------

def _process(nf, name):
    """The text of process `name` in the file `nf`, relative to the repository."""
    text = (ROOT / nf).read_text()
    found = re.search(rf'^process {name} \{{\n(.*?)^\}}', text, re.S | re.M)
    assert found, f'no process {name} in {nf}'
    return found.group(1)


INPUT = re.compile(r"\b(path|val)(?:\(\s*|\s+)(\w+)(?:,\s*stageAs:\s*'([^']+)')?")


def _inputs(body):
    """{name: (path or val, stageAs pattern or None)} of a process's inputs."""
    block = re.search(r'^    input:\n(.*?)^    output:', body, re.S | re.M)
    if block is None:
        return {}
    block = re.sub(r'//.*', '', block.group(1))
    return {name: (kind, stage or None) for kind, name, stage in INPUT.findall(block)}


def _script(body):
    """The Groovy that runs before the script's string, and the string."""
    found = re.search(r'^    script:\n(.*?)^    """\n(.*?)\n    """\n', body, re.S | re.M)
    return found.group(1), found.group(2)


TOKEN = re.compile(r"\s*(\?:|[?:()]|'[^']*'|[A-Za-z_][\w.]*)")


def _evaluate(expr, values):
    """The little Groovy these scripts interpolate: names, a.b, 'text', the
    Elvis operator, a ? b : c and parentheses, with Groovy's idea of false
    (null, '', an empty list, false). Anything else fails, rather than give
    a different answer from Groovy."""
    tokens = TOKEN.findall(expr)
    if ''.join(tokens).replace(' ', '') != expr.replace(' ', ''):
        raise ValueError(f'cannot evaluate {expr!r}')
    pos = 0

    def peek():
        return tokens[pos] if pos < len(tokens) else None

    def take(expected=None):
        nonlocal pos
        token = tokens[pos]
        if expected is not None and token != expected:
            raise ValueError(f'expected {expected!r} in {expr!r}')
        pos += 1
        return token

    def primary():
        token = take()
        if token == '(':
            value = ternary()
            take(')')
            return value
        if token.startswith("'"):
            return token[1:-1]
        value = values
        for part in token.split('.'):
            value = value[part]
        return value

    def elvis():
        value = primary()
        if peek() == '?:':
            take()
            other = elvis()
            return value if value not in (None, '', [], False) else other
        return value

    def ternary():
        value = elvis()
        if peek() == '?':
            take()
            yes = ternary()
            take(':')
            no = ternary()
            return yes if value not in (None, '', [], False) else no
        return value

    value = ternary()
    if pos != len(tokens):
        raise ValueError(f'cannot evaluate {expr!r}')
    return value


def _text(value):
    """A value as Groovy writes it into a string."""
    if isinstance(value, bool):
        return 'true' if value else 'false'
    if value is None:
        return 'null'
    if value == []:
        return ''
    return str(value)


def _render(template, values):
    """Groovy's triple-quoted string, for what these scripts use: \\\\, \\$, a
    backslash that ends a line and ${...}. `task.ext.code` renders empty."""
    out, i = [], 0
    while i < len(template):
        c = template[i]
        if c == '\\':
            nxt = template[i + 1]
            if nxt not in '\\$\n':
                raise ValueError(f'escape \\{nxt} not handled here')
            out.append('' if nxt == '\n' else nxt)
            i += 2
        elif c == '$':
            if template[i + 1] != '{':
                raise ValueError(f'bare $ at {template[i:i + 20]!r}')
            depth, j = 1, i + 2
            while depth:
                depth += {'{': 1, '}': -1}.get(template[j], 0)
                j += 1
            expr = template[i + 2:j - 1]
            out.append('' if expr.startswith('task.ext.code') else _text(_evaluate(expr, values)))
            i = j
        else:
            out.append(c)
            i += 1
    return ''.join(out)


def _staged_name(pattern, name):
    """Where Nextflow stages one file under a stageAs pattern: a * in the
    directories is its index, 1, a final * its name."""
    parent, _, last = pattern.rpartition('/')
    last = name if last == '*' else last
    return f"{parent.replace('*', '1')}/{last}" if parent else last


def _run(tmp_path, nf, name, inputs, params, tools):
    """Run process `name`'s script in a task directory under the launch
    folder's work/. `inputs` is what the task gets: for a path input, the
    file to stage, or [] for nothing staged; a value otherwise. Inputs the
    process does not declare are left out, as Nextflow would refuse them, and
    a value the process takes that `inputs` does not give is its own name:
    these tests are about the files. `tools` maps a command to the Python
    source of its stand-in. (exit code, stderr, task directory)."""
    body = _process(nf, name)
    declared = _inputs(body)
    launch = tmp_path / 'runs' / 'after'
    work = launch / 'work' / 'ab' / 'cdef0123'
    work.mkdir(parents=True)

    values = {'params': params, 'task': {'cpus': 2}, 'projectDir': str(ROOT)}
    for input_name, (kind, stage) in declared.items():
        given = inputs[input_name] if kind == 'path' else inputs.get(input_name, input_name)
        if kind == 'path' and given != []:
            staged = work / _staged_name(stage or '*', given.name)
            staged.parent.mkdir(parents=True, exist_ok=True)
            staged.symlink_to(given)
            given = staged.relative_to(work).as_posix()
        values[input_name] = given

    groovy, template = _script(body)
    for line in groovy.splitlines():
        line = line.strip()
        if line and not line.startswith('//'):
            found = re.fullmatch(r'def (\w+) = (.+)', line)
            assert found, f'{name}: cannot run the Groovy line {line!r}'
            values[found.group(1)] = _evaluate(found.group(2), values)
    script = _render(textwrap.dedent(template) + '\n', values)
    (work / '.command.sh').write_text(script)

    bin_dir = tmp_path / 'tools'
    bin_dir.mkdir()
    for command, source in tools.items():
        (bin_dir / f'{command}.py').write_text(textwrap.dedent(source))
        (bin_dir / command).write_text(f'#!/bin/bash\nexec "{sys.executable}" "{bin_dir / command}.py" "$@"\n')
        (bin_dir / command).chmod(0o755)
    env = dict(os.environ, PATH=f'{bin_dir}{os.pathsep}{os.environ["PATH"]}')
    # A script with its own #! line runs as it is; any other gets the
    # shell nextflow.config gives every task, bash -euo pipefail
    shell = ['bash'] if script.startswith('#!') else ['bash', '-euo', 'pipefail']
    done = subprocess.run(shell + ['.command.sh'], cwd=work, env=env, capture_output=True, text=True)
    return done.returncode, done.stderr, work


def _write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def _unseen(tmp_path, name):
    """Where a file is on the host, in a folder the task cannot see."""
    return str(tmp_path / 'host' / 'references' / name)


def _container(tmp_path, name, text):
    """A file inside the container, which the task reads where it is."""
    return _write(tmp_path / 'container' / 'reference' / name, text)


# ---------------------------------------------------------------------------
# BAKTA_ANNOTATOR
# ---------------------------------------------------------------------------

# Bakta 1.12's own check of --proteins (bakta/config.py, check_user_proteins)
FAKE_BAKTA = '''
    import os, sys
    from pathlib import Path
    args = sys.argv[1:]
    def option(name):
        return args[args.index(name) + 1]
    if not os.path.isdir(option('--db')):
        sys.exit('ERROR: wrong database path! --db=' + option('--db'))
    proteins = option('--proteins')
    path = Path(proteins).resolve()
    if not (path.is_file() and os.access(path, os.R_OK) and path.stat().st_size > 0):
        sys.exit(f'ERROR: user proteins file ({proteins}) not valid!')
    out = Path(option('--output'))
    out.mkdir(exist_ok=True)
    for ext in ('faa', 'tsv'):
        (out / f"{option('--prefix')}.{ext}").write_text('')
    # What it annotated with, for the test to look at
    (out / 'proteins.txt').write_text(path.read_text())
'''

K12 = 'LOCUS       NC_000913               4641652 bp    DNA     circular CON\n'
H37RV = 'LOCUS       AL123456                4411532 bp    DNA     linear   CON\n'


def _bakta(tmp_path, params, proteins):
    """Run BAKTA_ANNOTATOR on Sakai with the container's database."""
    params = dict({'bakta_db': str(tmp_path / 'container' / 'db-light'),
                   'ref_annotation': None,
                   'bakta_proteins': '/reference/sequence.gb'}, **params)
    (tmp_path / 'container' / 'db-light').mkdir(parents=True)
    inputs = {'sampleId': 'Sakai', 'genus': 'Escherichia', 'species': 'Escherichia coli',
              'fasta': _write(tmp_path / 'store' / 'sakai.fna', '>NC_002695.2\nACGT\n'),
              'summary_log': _write(tmp_path / 'store' / 'empty_summary.log', ''),
              'bakta_db': [], 'proteins': proteins}
    code, err, work = _run(tmp_path, 'modules/local/bakta.nf', 'BAKTA_ANNOTATOR',
                           inputs, params, {'bakta': FAKE_BAKTA})
    used = work / 'Sakai' / 'proteins.txt'
    return code, err, (used.read_text() if used.exists() else None)


class TestBakta:

    def test_the_custom_reference_reaches_bakta(self, tmp_path):
        # As seen on 2026-10-18: K-12 as the reference for Sakai, given
        # relative to the launch folder
        k12 = _write(tmp_path / 'runs' / 'inputs' / 'ecoli' / 'k12.gbff', K12)
        code, err, used = _bakta(tmp_path, {'ref_annotation': '../inputs/ecoli/k12.gbff'}, k12)
        assert 'not valid' not in err
        assert code == 0, err
        assert used == K12

    def test_bakta_proteins_on_the_host_reach_bakta(self, tmp_path):
        own = _write(tmp_path / 'store' / 'my_h37rv.gb', H37RV)
        code, err, used = _bakta(tmp_path, {'bakta_proteins': _unseen(tmp_path, 'my_h37rv.gb')}, own)
        assert code == 0, err
        assert used == H37RV

    def test_the_default_is_read_inside_the_container(self, tmp_path):
        # Not on the host, so nothing is staged and the path is the container's
        default = _container(tmp_path, 'sequence.gb', H37RV)
        code, err, used = _bakta(tmp_path, {'bakta_proteins': str(default)}, [])
        assert code == 0, err
        assert used == H37RV

    def test_the_custom_reference_wins_over_bakta_proteins(self, tmp_path):
        _container(tmp_path, 'sequence.gb', H37RV)
        k12 = _write(tmp_path / 'store' / 'k12.gbff', K12)
        code, err, used = _bakta(tmp_path, {'ref_annotation': _unseen(tmp_path, 'k12.gbff'),
                                            'bakta_proteins': str(tmp_path / 'container' / 'reference' / 'sequence.gb')},
                                 k12)
        assert code == 0, err
        assert used == K12

    def test_a_file_named_after_the_sample_cannot_collide(self, tmp_path):
        # Bakta writes Sakai/, so a reference called Sakai staged next to it
        # would be in the way
        reference = _write(tmp_path / 'store' / 'Sakai', K12)
        code, err, used = _bakta(tmp_path, {'ref_annotation': _unseen(tmp_path, 'Sakai')}, reference)
        assert code == 0, err
        assert used == K12


# ---------------------------------------------------------------------------
# PATHOTYPR
# ---------------------------------------------------------------------------

FAKE_PATHOTYPR = '''
    import os, sys
    from pathlib import Path
    args = sys.argv[1:]
    def option(name):
        return args[args.index(name) + 1]
    markers, reference = option('-m'), option('-r')
    if not os.path.isfile(markers):
        sys.exit(f'Error: The TSV positions file {markers} does not exist')
    if not os.path.isfile(reference):
        sys.exit(f'Error: The reference file {reference} does not exist')
    genome = open(option('-i')).readline()[1:].split()[0]
    prefix = option('--output-prefix')
    Path(f'{prefix}_summary.tsv').write_text(
        f'genome\\tlineage:count\\tmajor_lineage\\n{genome}\\tL4:61 L4;L4.9:6\\tL4\\n')
    # What it typed with, for the test to look at
    Path('markers.txt').write_text(Path(markers).read_text() + Path(reference).read_text())
'''

MARKERS = 'lineage\tposition\tallele\nL4\t1234\tA\n'
ANCESTOR = '>MTB_anc\nACGTACGT\n'


def _pathotypr(tmp_path, params, snps, ref):
    """Run PATHOTYPR on one genome; (exit code, stderr, task directory)."""
    params = dict({'pathotyper_snps': '/reference/snps_lineages.tsv',
                   'pathotyper_ref': '/reference/MTB_ancestor_reference.fasta'}, **params)
    inputs = {'sampleId': 'G1', 'fasta': _write(tmp_path / 'store' / 'G1.fasta', '>G1\nACGT\n'),
              'pathotyper_snps': snps, 'pathotyper_ref': ref}
    return _run(tmp_path, 'modules/local/pathotypr.nf', 'PATHOTYPR', inputs, params,
                {'pathotypr': FAKE_PATHOTYPR})


def _typed(work):
    """The genome PATHOTYPR wrote, by the name it gave it."""
    return sorted(p.name for p in work.glob('*.fasta') if not p.is_symlink())


class TestPathotypr:

    def test_markers_on_the_host_reach_pathotypr(self, tmp_path):
        snps = _write(tmp_path / 'store' / 'snps_lineages.tsv', MARKERS)
        ref = _write(tmp_path / 'store' / 'MTB_ancestor_reference.fasta', ANCESTOR)
        code, err, work = _pathotypr(tmp_path, {
            'pathotyper_snps': _unseen(tmp_path, 'snps_lineages.tsv'),
            'pathotyper_ref': _unseen(tmp_path, 'MTB_ancestor_reference.fasta')}, snps, ref)
        assert 'does not exist' not in err
        assert _typed(work) == ['L4_G1.fasta']
        assert (work / 'markers.txt').read_text() == MARKERS + ANCESTOR

    def test_the_defaults_are_read_inside_the_container(self, tmp_path):
        snps = _container(tmp_path, 'snps_lineages.tsv', MARKERS)
        ref = _container(tmp_path, 'MTB_ancestor_reference.fasta', ANCESTOR)
        code, err, work = _pathotypr(tmp_path, {'pathotyper_snps': str(snps),
                                                'pathotyper_ref': str(ref)}, [], [])
        assert _typed(work) == ['L4_G1.fasta']
        assert (work / 'markers.txt').read_text() == MARKERS + ANCESTOR

    def test_one_on_the_host_and_one_default(self, tmp_path):
        snps = _write(tmp_path / 'store' / 'my_markers.tsv', MARKERS)
        ref = _container(tmp_path, 'MTB_ancestor_reference.fasta', ANCESTOR)
        code, err, work = _pathotypr(tmp_path, {'pathotyper_snps': _unseen(tmp_path, 'my_markers.tsv'),
                                                'pathotyper_ref': str(ref)}, snps, [])
        assert _typed(work) == ['L4_G1.fasta']


# ---------------------------------------------------------------------------
# GET_OUTPUTS
# ---------------------------------------------------------------------------

# Every Python script GET_OUTPUTS runs, writing down how it was called
FAKE_PYTHON = '''
    import sys
    with open('python3.calls', 'a') as calls:
        calls.write('\\t'.join(sys.argv[1:]) + '\\n')
'''

MYCOBROWSER = 'Locus\tName\tFunctional_Category\nRv0001\tdnaA\tinformation pathways\n'


def _get_outputs(tmp_path, params, mycobrowser_ref):
    """Run GET_OUTPUTS on one sample; (exit code, stderr, the Mycobrowser
    argument generate_tsv_gff_gbk.py got, task directory)."""
    params = dict({'mycobrowser_ref': '/reference/Mycobacterium_tuberculosis_H37Rv_txt_v4.txt',
                   'skip_mycobrowser': False}, **params)
    store = tmp_path / 'store'
    inputs = {'sampleId': 'S1',
              'amap_tsv': _write(store / 'S1.AMAP.tsv', 'id\tstart\nRv0001\t1\n'),
              'amap_faa': _write(store / 'S1.AMAP.faa', '>Rv0001\nM\n'),
              'amap_fna': _write(store / 'S1.AMAP.fna', '>Rv0001\nATG\n'),
              'bakta_tsv': _write(store / 'S1.tsv', ''),
              'fasta': _write(store / 'S1.fasta', '>S1\nACGT\n'),
              'emapper': _write(store / 'S1.emapper.annotations', ''),
              'mycobrowser_ref': mycobrowser_ref}
    code, err, work = _run(tmp_path, 'modules/local/get.nf', 'GET_OUTPUTS', inputs, params,
                           {'python3': FAKE_PYTHON})
    calls = [line.split('\t') for line in (work / 'python3.calls').read_text().splitlines()]
    generate = next(call for call in calls if call[0].endswith('generate_tsv_gff_gbk.py'))
    return code, err, generate[4], work


class TestGetOutputs:

    def test_mycobrowser_on_the_host_reaches_the_master_table(self, tmp_path):
        table = _write(tmp_path / 'store' / 'mycobrowser_v5.txt', MYCOBROWSER)
        code, err, myco, work = _get_outputs(
            tmp_path, {'mycobrowser_ref': _unseen(tmp_path, 'mycobrowser_v5.txt')}, table)
        assert 'skipping Mycobrowser' not in err
        assert code == 0, err
        assert myco != 'NONE'
        assert (work / myco).read_text() == MYCOBROWSER

    def test_the_default_is_read_inside_the_container(self, tmp_path):
        default = _container(tmp_path, 'Mycobacterium_tuberculosis_H37Rv_txt_v4.txt', MYCOBROWSER)
        code, err, myco, _ = _get_outputs(tmp_path, {'mycobrowser_ref': str(default)}, [])
        assert code == 0, err
        assert myco == str(default)

    def test_skip_mycobrowser_still_skips_it(self, tmp_path):
        table = _write(tmp_path / 'store' / 'mycobrowser_v5.txt', MYCOBROWSER)
        code, err, myco, _ = _get_outputs(
            tmp_path, {'mycobrowser_ref': _unseen(tmp_path, 'mycobrowser_v5.txt'),
                       'skip_mycobrowser': True}, table)
        assert code == 0, err
        assert myco == 'NONE'


# ---------------------------------------------------------------------------
# BUNDLED_REFERENCE
# ---------------------------------------------------------------------------

PROTEOME = '>Rv0001\nMTDDPGSGFTTVWNAVVSELNGDPKVDDGP\n'
GENOME = '>NC_000962.3\nTTGACCGATGACCCCGGTTC\n'
BAD = '>bad_1\nMKV\n'


def _bundled(tmp_path, params, inputs):
    defaults = {'miniprot_ref_aa': str(_container(tmp_path, 'aa_H37Rv.translated.fas', PROTEOME)),
                'miniprot_ref_h37rv': str(_container(tmp_path, 'new_refH37Rv.fasta', GENOME)),
                'miniprot_ref_bad_bakta': str(_container(tmp_path, 'bad.bakta.translated.fas', BAD))}
    inputs = dict({'miniprot_ref_aa': [], 'miniprot_ref_h37rv': [], 'miniprot_ref_bad_bakta': []}, **inputs)
    code, err, work = _run(tmp_path, 'modules/local/reference.nf', 'BUNDLED_REFERENCE', inputs,
                           dict(defaults, **params), {})
    out = {name: (work / name).read_text() if (work / name).is_file() else None
           for name in ('ref_proteome.faa', 'ref_genome.fasta', 'bad_bakta.faa')}
    return code, err, out


class TestBundledReference:

    def test_the_defaults_are_copied_out_of_the_container(self, tmp_path):
        code, err, out = _bundled(tmp_path, {}, {})
        assert code == 0, err
        assert out == {'ref_proteome.faa': PROTEOME, 'ref_genome.fasta': GENOME, 'bad_bakta.faa': BAD}

    def test_one_of_them_on_the_host(self, tmp_path):
        # A newer proteome of your own, with the genome and the list of bad
        # calls that ship in the container
        own = '>Rv0001\nMTDDPGSGFTTVWNAVVSELNGDPKVDDGPSSDANLSAPLTPQQRAWLNLVQPLTIVEGFALLSVPSSFVQNEIERHLRAPITDALSRRLGHQIQLGVRIAPPATDEADDTTVPPSENPATTSPDTTTDNDEIDDSAAARGDNQHSWNPSHGVSV\n'
        proteome = _write(tmp_path / 'store' / 'my_proteome.faa', own)
        code, err, out = _bundled(tmp_path, {'miniprot_ref_aa': _unseen(tmp_path, 'my_proteome.faa')},
                                  {'miniprot_ref_aa': proteome})
        assert 'not found' not in err
        assert code == 0, err
        assert out == {'ref_proteome.faa': own, 'ref_genome.fasta': GENOME, 'bad_bakta.faa': BAD}

    def test_one_named_like_an_output_is_copied_not_written_over(self, tmp_path):
        genome = _write(tmp_path / 'store' / 'ref_genome.fasta', GENOME)
        code, err, out = _bundled(tmp_path, {'miniprot_ref_h37rv': _unseen(tmp_path, 'ref_genome.fasta')},
                                  {'miniprot_ref_h37rv': genome})
        assert code == 0, err
        assert out['ref_genome.fasta'] == GENOME
        assert genome.read_text() == GENOME


# ---------------------------------------------------------------------------
# No task reads a reference file by the path it was given
# ---------------------------------------------------------------------------

def _reference_params():
    """The params that name a reference file: the custom reference, and those
    whose default lives in the container, under /reference."""
    text = (ROOT / 'main.nf').read_text()
    return {'ref_annotation'} | set(re.findall(r'^params\.(\w+) = "/reference/', text, re.M))


def _scripts():
    """(file, process, script block) of every process."""
    found = []
    for nf in NF_FILES:
        for m in re.finditer(r'^process (\w+) \{\n(.*?)^\}', nf.read_text(), re.S | re.M):
            script = re.split(r'^[ \t]*script:[ \t]*\n', m.group(2), maxsplit=1, flags=re.M)[1]
            found.append((nf.relative_to(ROOT).as_posix(), m.group(1), script))
    return found


# TEST_GENOME copies the genome the container ships, which is what -profile
# test annotates
READS_THE_CONTAINER = {'TEST_GENOME'}


def test_finds_the_reference_params():
    assert {'ref_annotation', 'bakta_proteins', 'pathotyper_snps', 'pathotyper_ref', 'mycobrowser_ref',
            'miniprot_ref_aa', 'miniprot_ref_h37rv', 'miniprot_ref_bad_bakta'} <= _reference_params()


@pytest.mark.parametrize('nf,name,script', _scripts(), ids=[name for _, name, _ in _scripts()])
def test_reference_params_are_only_a_fallback(nf, name, script):
    # Each may be read only where the task got nothing staged: after `?:`,
    # behind the input that carries the file when it is on the host
    if name in READS_THE_CONTAINER:
        return
    direct = sorted({m.group(1) for m in re.finditer(r'\bparams\.(\w+)', script)
                     if m.group(1) in _reference_params()
                     and not script[:m.start()].rstrip().rstrip('(').rstrip().endswith('?:')})
    assert not direct, (
        f"{name} ({nf}) reads {', '.join('params.' + p for p in direct)} by path, which a "
        f"container cannot see when the file is on the host; take it as a path input, "
        f"staged from the host, with the param as its fallback (see main.nf, stagedIfOnHost)")


@pytest.mark.parametrize('nf,name', [('modules/local/bakta.nf', 'BAKTA_ANNOTATOR'),
                                     ('modules/local/pathotypr.nf', 'PATHOTYPR'),
                                     ('modules/local/get.nf', 'GET_OUTPUTS'),
                                     ('modules/local/reference.nf', 'BUNDLED_REFERENCE')])
def test_each_file_is_staged_in_a_directory_of_its_own(nf, name):
    # So that it cannot collide with what the task writes, named after the
    # sample, nor with another input of the same name
    staged = {n: stage for n, (kind, stage) in _inputs(_process(nf, name)).items()
              if kind == 'path' and stage}
    assert staged, f'{name} stages no reference file'
    dirs = [stage.split('/')[0] for stage in staged.values()]
    assert all('/' in stage for stage in staged.values()), staged
    assert len(set(dirs)) == len(dirs), staged


# ---------------------------------------------------------------------------
# validate_inputs.py --reference
# ---------------------------------------------------------------------------

def _genome(n_bp=600_000):
    """A FASTA every check passes without a warning: 70% GC, no N."""
    seq = 'GGCCAATGCC' * (n_bp // 10)
    return '>chromosome\n' + '\n'.join(seq[i:i + 80] for i in range(0, len(seq), 80)) + '\n'


def _check(tmp_path, *args):
    """validate_inputs.py on a sheet whose genome passes, with these options."""
    _write(tmp_path / 'G1.fasta', _genome())
    _write(tmp_path / 'samples.tsv', 'sampleId\tfasta\nG1\tG1.fasta\n')
    done = subprocess.run([sys.executable, str(SCRIPT), '--tsv', 'samples.tsv', *args],
                          cwd=tmp_path, capture_output=True, text=True)
    return done.returncode, done.stdout + done.stderr


class TestValidateReferences:

    def test_a_reference_in_the_container_passes(self, tmp_path):
        snps = _container(tmp_path, 'snps_lineages.tsv', MARKERS)
        code, out = _check(tmp_path, '--reference', 'pathotyper_snps', str(snps))
        assert code == 0, out
        assert f'✓ --pathotyper_snps: {snps}' in out

    def test_one_that_is_nowhere_stops_the_run(self, tmp_path):
        # A typo, or a default without the container
        typo = _unseen(tmp_path, 'snps_lineage.tsv')
        code, out = _check(tmp_path, '--reference', 'pathotyper_snps', typo)
        assert code == 2
        assert f'--pathotyper_snps: {typo} not found, neither on the host nor inside the container' in out

    def test_a_directory_stops_the_run(self, tmp_path):
        (tmp_path / 'reference').mkdir()
        code, out = _check(tmp_path, '--reference', 'bakta_proteins', str(tmp_path / 'reference'))
        assert code == 2
        assert f"--bakta_proteins: {tmp_path / 'reference'} is not a file" in out

    def test_an_empty_one_stops_the_run(self, tmp_path):
        empty = _container(tmp_path, 'sequence.gb', '')
        code, out = _check(tmp_path, '--reference', 'bakta_proteins', str(empty))
        assert code == 2
        assert f'--bakta_proteins: {empty} is empty' in out

    def test_each_is_checked(self, tmp_path):
        snps = _container(tmp_path, 'snps_lineages.tsv', MARKERS)
        typo = _unseen(tmp_path, 'MTB_ancestor.fasta')
        code, out = _check(tmp_path, '--reference', 'pathotyper_snps', str(snps),
                           '--reference', 'pathotyper_ref', typo)
        assert code == 2
        assert f'✓ --pathotyper_snps: {snps}' in out
        assert f'✗ --pathotyper_ref: {typo}' in out
        assert 'Total: 1 error(s)' in out

    def test_a_missing_mycobrowser_table_only_warns(self, tmp_path):
        # GET_OUTPUTS goes on without it and leaves the categories out
        typo = _unseen(tmp_path, 'mycobrowser.txt')
        code, out = _check(tmp_path, '--reference', 'mycobrowser_ref', typo)
        assert code == 0, out
        assert (f'⚠ --mycobrowser_ref: {typo} not found, neither on the host nor inside the '
                f'container; the master table will have no Mycobrowser categories') in out
        code, out = _check(tmp_path / 'strict', '--reference', 'mycobrowser_ref', typo, '--strict')
        assert code == 1

    def test_without_any_there_is_nothing_to_check(self, tmp_path):
        code, out = _check(tmp_path)
        assert code == 0, out
        assert 'Reference Files' not in out

    def test_the_module_passes_each_one_on(self):
        body = _process('modules/local/validation.nf', 'VALIDATE_INPUTS')
        assert _inputs(body).get('references') == ('val', None)
        assert re.search(r"references\.each \{ param, given -> options << '--reference ' \+ "
                         r"quoted\.call\(param\) \+ ' ' \+ quoted\.call\(given\) \}", body)
