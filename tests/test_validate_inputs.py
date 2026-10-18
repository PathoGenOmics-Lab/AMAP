#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VALIDATE_INPUTS reads every genome the sample sheet names, wherever it lives.

The sheet gives each genome's path on the host, absolute or relative to where
Nextflow was launched. VALIDATE_INPUTS staged the sheet alone and opened those
paths from inside its task, where a relative path points somewhere else and,
under -profile docker, only the folder that the work directory and the sheet
share is mounted. With the sheet and the work directory in runs/after and the
genomes in runs/inputs, all of them were "FASTA file not found" and the run
stopped before anything ran. The genomes are now inputs of VALIDATE_INPUTS,
which reads each one where Nextflow staged it.

The module tests run VALIDATE_INPUTS's script the way Nextflow renders it, in
a task directory laid out the way Nextflow lays one out, so they also fail if
modules/local/validation.nf goes back to opening the sheet's paths. A folder
the task cannot see, as a container sees no folder Nextflow does not mount
for it, is one that does not exist.
"""

import io
import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'bin'))

SCRIPT = ROOT / 'bin' / 'validate_inputs.py'
MODULE = ROOT / 'modules' / 'local' / 'validation.nf'


def _genome(n_bp=600_000):
    """A FASTA every check passes without a warning: 70% GC, no N."""
    seq = 'GGCCAATGCC' * (n_bp // 10)
    return '>chromosome\n' + '\n'.join(seq[i:i + 80] for i in range(0, len(seq), 80)) + '\n'


def _write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def _sheet(path, rows):
    return _write(path, 'sampleId\tfasta\n' + ''.join(f'{s}\t{f}\n' for s, f in rows))


def _unseen(tmp_path, name):
    """Where a genome is on the host, in a folder the task cannot see."""
    return str(tmp_path / 'host' / 'genomes' / name)


def _validate(cwd, *args, stdin=None):
    """Run validate_inputs.py; (exit code, what it printed)."""
    done = subprocess.run([sys.executable, str(SCRIPT), *args], cwd=cwd, input=stdin,
                          capture_output=True, text=True)
    return done.returncode, done.stdout + done.stderr


# ---------------------------------------------------------------------------
# validate_inputs.py --genomes
# ---------------------------------------------------------------------------

class TestGenomeTable:

    def test_each_genome_is_read_from_the_table(self, tmp_path):
        _sheet(tmp_path / 'samples.tsv', [('G1', _unseen(tmp_path, 'G1.fasta')),
                                          ('G2', _unseen(tmp_path, 'G2.fasta'))])
        _write(tmp_path / 'genomes' / '1' / 'G1.fasta', _genome())
        _write(tmp_path / 'genomes' / '2' / 'G2.fasta', _genome())
        _sheet(tmp_path / 'staged.tsv', [('G1', 'genomes/1/G1.fasta'),
                                         ('G2', 'genomes/2/G2.fasta')])
        code, out = _validate(tmp_path, '--tsv', 'samples.tsv', '--genomes', 'staged.tsv')
        assert code == 0, out
        assert 'Results: 2 OK, 0 warnings, 0 failed' in out

    def test_the_table_can_come_on_standard_input(self, tmp_path):
        _sheet(tmp_path / 'samples.tsv', [('G1', _unseen(tmp_path, 'G1.fasta'))])
        _write(tmp_path / 'genomes' / '1' / 'G1.fasta', _genome())
        code, out = _validate(tmp_path, '--tsv', 'samples.tsv', '--genomes', '-',
                              stdin='sampleId\tfasta\nG1\tgenomes/1/G1.fasta\n')
        assert code == 0, out
        assert 'Results: 1 OK, 0 warnings, 0 failed' in out

    def test_messages_name_the_path_in_the_sheet(self, tmp_path):
        on_host = _unseen(tmp_path, 'G1.fasta')
        _sheet(tmp_path / 'samples.tsv', [('G1', on_host)])
        _write(tmp_path / 'genomes' / '1' / 'G1.fasta', '')
        _sheet(tmp_path / 'staged.tsv', [('G1', 'genomes/1/G1.fasta')])
        code, out = _validate(tmp_path, '--tsv', 'samples.tsv', '--genomes', 'staged.tsv')
        assert code == 2
        assert f'FASTA file is empty: {on_host}' in out
        assert 'genomes/1' not in out

    def test_a_sample_not_in_the_table_is_read_from_the_sheet(self, tmp_path):
        _write(tmp_path / 'G2.fasta', _genome())
        _sheet(tmp_path / 'samples.tsv', [('G1', _unseen(tmp_path, 'G1.fasta')),
                                          ('G2', 'G2.fasta')])
        _write(tmp_path / 'genomes' / '1' / 'G1.fasta', _genome())
        _sheet(tmp_path / 'staged.tsv', [('G1', 'genomes/1/G1.fasta')])
        code, out = _validate(tmp_path, '--tsv', 'samples.tsv', '--genomes', 'staged.tsv')
        assert code == 0, out
        assert 'Results: 2 OK, 0 warnings, 0 failed' in out

    def test_without_a_table_the_sheet_is_read_as_before(self, tmp_path):
        _write(tmp_path / 'G1.fasta', _genome())
        _sheet(tmp_path / 'samples.tsv', [('G1', 'G1.fasta'), ('G2', 'missing.fasta')])
        code, out = _validate(tmp_path, '--tsv', 'samples.tsv')
        assert code == 2
        assert 'FASTA file not found: missing.fasta' in out
        assert 'Results: 1 OK, 0 warnings, 1 failed' in out

    def test_quotes_and_spaces_are_taken_as_they_are(self):
        from validate_inputs import read_genome_table
        table = 'sampleId\tfasta\nG1\tgenomes/1/my "best" genome.fasta\nG\'2\tgenomes/2/it\'s.fa\n'
        assert read_genome_table(io.StringIO(table)) == {
            'G1': 'genomes/1/my "best" genome.fasta', "G'2": "genomes/2/it's.fa"}


# ---------------------------------------------------------------------------
# VALIDATE_INPUTS's own script
# ---------------------------------------------------------------------------

def _module():
    """VALIDATE_INPUTS's script block as Nextflow strips it, and the pattern its
    genomes are staged with (None when it stages none)."""
    text = MODULE.read_text()
    script = re.search(r'^    script:\n.*?^    """\n(.*?)\n    """\n', text, re.S | re.M)
    stage = re.search(r"path\(genomes, stageAs: '([^']+)'\)", text)
    return textwrap.dedent(script.group(1)) + '\n', stage.group(1) if stage else None


def _render(script, values):
    """Groovy's triple-quoted string, for the little this script uses: \\\\, \\$,
    a backslash that ends a line and ${...}. Anything else fails, rather than
    render differently from Nextflow."""
    out, i = [], 0
    while i < len(script):
        c = script[i]
        if c == '\\':
            nxt = script[i + 1]
            if nxt not in '\\$\n':
                raise ValueError(f'escape \\{nxt} not handled here')
            out.append('' if nxt == '\n' else nxt)
            i += 2
        elif c == '$':
            if script[i + 1] != '{':
                raise ValueError(f'bare $ at {script[i:i + 20]!r}')
            depth, j = 1, i + 2
            while depth:
                depth += {'{': 1, '}': -1}.get(script[j], 0)
                j += 1
            expr = script[i + 2:j - 1]
            out.append('' if expr.startswith('task.ext.code') else values[expr])
            i = j
        else:
            out.append(c)
            i += 1
    return ''.join(out)


def _staged_name(pattern, index, name):
    """Where Nextflow stages the index-th file (from 1) of a list under a
    stageAs pattern: a * in the directories is the index, a final * the name."""
    parent, _, last = pattern.rpartition('/')
    last = name if last == '*' else last
    return f"{parent.replace('*', str(index))}/{last}" if parent else last


def _quoted(value):
    """A shell word in single quotes, as the module's Groovy writes each one."""
    return "'" + str(value).replace("'", "'\\''") + "'"


def _run_module(tmp_path, rows, sources):
    """Run VALIDATE_INPUTS's script in a task directory under the launch
    folder's work/, the sheet and each genome staged as a link the way
    Nextflow stages them. `rows` is the sheet, (sampleId, path as written);
    `sources` the genome each row's path stands for, which the task reaches
    only where the module stages it. (exit code, the log, task directory)."""
    script, stage = _module()
    launch = tmp_path / 'runs' / 'after'
    sheet = _sheet(launch / 'samples.tsv', rows)
    work = launch / 'work' / 'ab' / 'cdef0123'
    work.mkdir(parents=True)
    (work / sheet.name).symlink_to(sheet)

    genome_rows = []
    if stage is not None:
        for index, ((sample, _), source) in enumerate(zip(rows, sources), start=1):
            staged = work / _staged_name(stage, index, source.name)
            staged.parent.mkdir(parents=True, exist_ok=True)
            staged.symlink_to(source)
            genome_rows.append(f'{_quoted(sample)} {_quoted(staged.relative_to(work).as_posix())}')

    tools = tmp_path / 'tools'
    tools.mkdir()
    (tools / 'python3').write_text(f'#!/bin/bash\nexec "{sys.executable}" "$@"\n')
    (tools / 'python3').chmod(0o755)

    rendered = _render(script, {
        'tsv_file': sheet.name,
        'genome_rows': ' '.join(genome_rows),
        'bakta_arg': '',
        'eggnog_arg': '',
        'projectDir': str(ROOT),
    })
    (work / '.command.sh').write_text(rendered)
    env = dict(os.environ, PATH=f'{tools}{os.pathsep}{os.environ["PATH"]}')
    # nextflow.config runs every task with bash -euo pipefail
    done = subprocess.run(['bash', '-euo', 'pipefail', '.command.sh'], cwd=work, env=env,
                          capture_output=True, text=True)
    log = work / 'validation.log'
    return done.returncode, (log.read_text() if log.exists() else done.stderr), work


class TestModule:

    def test_genomes_in_a_folder_the_task_does_not_see(self, tmp_path):
        # As seen on 2026-10-18: sheet and work directory in runs/after,
        # genomes in a folder of their own
        sources = [_write(tmp_path / 'store' / f'G{n}.fasta', _genome()) for n in (1, 2, 3, 4)]
        rows = [(f'G{n}', _unseen(tmp_path, f'G{n}.fasta')) for n in (1, 2, 3, 4)]
        code, log, _ = _run_module(tmp_path, rows, sources)
        assert 'FASTA file not found' not in log
        assert 'Results: 4 OK, 0 warnings, 0 failed' in log
        assert code == 0

    def test_relative_paths_are_relative_to_the_launch_folder(self, tmp_path):
        # As documented, and as the rest of the pipeline reads them. The task
        # runs in work/ab/cdef0123, where the same path points elsewhere.
        sources = [_write(tmp_path / 'runs' / 'inputs' / 'G1.fasta', _genome()),
                   _write(tmp_path / 'runs' / 'after' / 'assemblies' / 'G2.fa', _genome())]
        rows = [('G1', '../inputs/G1.fasta'), ('G2', './assemblies/G2.fa')]
        code, log, _ = _run_module(tmp_path, rows, sources)
        assert 'Results: 2 OK, 0 warnings, 0 failed' in log
        assert code == 0

    def test_two_genomes_with_the_same_file_name(self, tmp_path):
        sources = [_write(tmp_path / 'store' / 'closed' / 'genome.fasta', _genome()),
                   _write(tmp_path / 'store' / 'draft' / 'genome.fasta', '>contig_1\nACGTACGT\n')]
        rows = [('CLOSED', _unseen(tmp_path, 'closed/genome.fasta')),
                ('DRAFT', _unseen(tmp_path, 'draft/genome.fasta'))]
        code, log, _ = _run_module(tmp_path, rows, sources)
        assert code == 2
        assert '[DRAFT] Genome too small: 8 bp' in log
        assert '[CLOSED]' not in log
        assert 'Results: 1 OK, 0 warnings, 1 failed' in log

    def test_a_bad_genome_still_stops_the_run_and_is_named_as_in_the_sheet(self, tmp_path):
        sources = [_write(tmp_path / 'runs' / 'inputs' / 'G1.fasta', _genome()),
                   _write(tmp_path / 'runs' / 'inputs' / 'G2.fasta', '')]
        rows = [('G1', '../inputs/G1.fasta'), ('G2', '../inputs/G2.fasta')]
        code, log, _ = _run_module(tmp_path, rows, sources)
        assert code == 2
        assert 'FASTA file is empty: ../inputs/G2.fasta' in log

    def test_quotes_and_spaces_in_sample_ids_and_file_names(self, tmp_path):
        sources = [_write(tmp_path / 'store' / "it's a genome.fasta", _genome()),
                   _write(tmp_path / 'store' / 'small one.fasta', '>c\nACGT\n')]
        rows = [("G'1", _unseen(tmp_path, "it's a genome.fasta")),
                ('G 2', _unseen(tmp_path, 'small one.fasta'))]
        code, log, _ = _run_module(tmp_path, rows, sources)
        assert code == 2
        assert '[G 2] Genome too small: 4 bp' in log
        assert 'Results: 1 OK, 0 warnings, 1 failed' in log

    def test_a_sheet_of_500_samples(self, tmp_path):
        # One genome named 500 times, staged 500 times as Nextflow does, but
        # each sample's own path: the last one is too small to pass. Under
        # 500 kb the others pass with a warning, and read five times faster.
        good = _write(tmp_path / 'store' / 'G.fasta', _genome(110_000))
        small = _write(tmp_path / 'store' / 'small.fasta', '>c\nACGT\n')
        sources = [good] * 499 + [small]
        rows = [(f'S{n:03d}', _unseen(tmp_path, f'S{n:03d}.fasta')) for n in range(1, 501)]
        code, log, work = _run_module(tmp_path, rows, sources)
        assert code == 2
        assert 'Samples: 500' in log
        assert '[S500] Genome too small: 4 bp' in log
        assert 'Results: 0 OK, 499 warnings, 1 failed' in log
        assert len(list((work / 'genomes').iterdir())) == 500

    def test_the_sheet_is_left_as_it_was(self, tmp_path):
        sources = [_write(tmp_path / 'runs' / 'inputs' / 'G1.fasta', _genome())]
        rows = [('G1', '../inputs/G1.fasta')]
        _run_module(tmp_path, rows, sources)
        sheet = tmp_path / 'runs' / 'after' / 'samples.tsv'
        assert sheet.read_text() == 'sampleId\tfasta\nG1\t../inputs/G1.fasta\n'

    def test_module_stages_the_genomes_one_directory_each(self):
        # Without a directory each, two genomes called genome.fasta collide
        _, stage = _module()
        assert stage is not None, 'VALIDATE_INPUTS does not stage the genomes'
        assert _staged_name(stage, 1, 'genome.fasta') != _staged_name(stage, 2, 'genome.fasta')
