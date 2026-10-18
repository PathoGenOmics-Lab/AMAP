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

The databases had the same blind spot. The eggNOG path, relative to where
Nextflow was launched, was looked for in the task directory, and one that is
not built yet, which EGGNOG_DB builds, was an error, so every run with
functional annotation stopped at validation. A Bakta database on the host was
"path does not exist", and the container's own drew a warning about a file no
Bakta database has.
"""

import io
import json
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
    """VALIDATE_INPUTS's script block as Nextflow strips it, and the stageAs
    pattern of each input that has one, by name."""
    text = MODULE.read_text()
    script = re.search(r'^    script:\n.*?^    """\n(.*?)\n    """\n', text, re.S | re.M)
    stages = dict(re.findall(r"path[( ](\w+), stageAs: '([^']+)'", text))
    return textwrap.dedent(script.group(1)) + '\n', stages


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


def _options(databases, strict):
    """The options the module's Groovy writes: each database by the path given,
    and where it is staged if it was, then --strict for --strict_validation."""
    options = []
    for name in ('bakta', 'eggnog'):
        if name in databases:
            given, staged = databases[name]
            options.append(f'--{name}-db {_quoted(given)}')
            if staged:
                options.append(f'--{name}-db-staged {_quoted(staged)}')
    if strict:
        options.append('--strict')
    return ' '.join(options)


def _run_module(tmp_path, rows, sources, databases=None, strict=False):
    """Run VALIDATE_INPUTS's script in a task directory under the launch
    folder's work/, the sheet and each genome staged as a link the way
    Nextflow stages them. `rows` is the sheet, (sampleId, path as written);
    `sources` the genome each row's path stands for, which the task reaches
    only where the module stages it. `databases` maps bakta and eggnog to
    (path given, the database on the host or None), staged likewise.
    (exit code, the log, task directory)."""
    script, stages = _module()
    stage = stages.get('genomes')
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

    staged_databases = {}
    for name, (given, source) in (databases or {}).items():
        staged = None
        if source is not None and stages.get(f'{name}_db'):
            staged = work / _staged_name(stages[f'{name}_db'], 1, source.name)
            staged.parent.mkdir(parents=True, exist_ok=True)
            staged.symlink_to(source)
            staged = staged.relative_to(work).as_posix()
        staged_databases[name] = (given, staged)

    tools = tmp_path / 'tools'
    tools.mkdir()
    (tools / 'python3').write_text(f'#!/bin/bash\nexec "{sys.executable}" "$@"\n')
    (tools / 'python3').chmod(0o755)

    rendered = _render(script, {
        'tsv_file': sheet.name,
        'genome_rows': ' '.join(genome_rows),
        "options.join(' ')": _options(staged_databases, strict),
        # What the module wrote before it staged the databases
        'bakta_arg': f"--bakta-db {_quoted(databases['bakta'][0])}" if 'bakta' in (databases or {}) else '',
        'eggnog_arg': f"--eggnog-db {_quoted(databases['eggnog'][0])}" if 'eggnog' in (databases or {}) else '',
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
        stage = _module()[1].get('genomes')
        assert stage is not None, 'VALIDATE_INPUTS does not stage the genomes'
        assert _staged_name(stage, 1, 'genome.fasta') != _staged_name(stage, 2, 'genome.fasta')


# ---------------------------------------------------------------------------
# The databases
# ---------------------------------------------------------------------------

LIGHT = {'date': '2025-02-24', 'major': 6, 'minor': 0, 'type': 'light',
         'software-min': {'major': 1, 'minor': 11}}


def _bakta_db(path, files=('bakta.db', 'pfam.h3f', 'pscc.dmnd')):
    """A stand-in for a Bakta database, laid out as the container's light one:
    version.json and some of its files, all empty."""
    path.mkdir(parents=True)
    (path / 'version.json').write_text(json.dumps(LIGHT))
    for name in files:
        (path / name).write_text('')
    return path


def _eggnog_db(path, dmnd_bytes=None, annotations=True, built=True):
    """A stand-in for eggNOG's database: bacteria.dmnd of that many bytes (a
    sparse file, so the real size costs nothing), eggnog.db and the
    eggnog_db.ok that EGGNOG_DB's storeDir leaves."""
    path.mkdir(parents=True)
    if dmnd_bytes is not None:
        with open(path / 'bacteria.dmnd', 'wb') as f:
            f.truncate(dmnd_bytes)
    if annotations:
        (path / 'eggnog.db').write_text('')
    if built:
        (path / 'eggnog_db.ok').write_text('')
    return path


BUILT = 100_000_001


def _check(tmp_path, *args):
    """validate_inputs.py on a sheet whose genome passes, with these options."""
    _write(tmp_path / 'G1.fasta', _genome())
    _sheet(tmp_path / 'samples.tsv', [('G1', 'G1.fasta')])
    return _validate(tmp_path, '--tsv', 'samples.tsv', *args)


class TestBaktaDb:

    def test_the_container_database_passes_without_a_warning(self, tmp_path):
        # Its files are bakta.db and version.json; no Bakta database has a
        # db.sqlite3, which drew a warning on every run
        db = _bakta_db(tmp_path / 'db-light')
        code, out = _check(tmp_path, '--bakta-db', str(db))
        assert code == 0, out
        assert f'✓ Bakta DB: {db} (light 6.0)' in out
        assert '⚠' not in out

    def test_one_on_the_host_is_read_where_staged(self, tmp_path):
        on_host = _unseen(tmp_path, 'db-light')
        _bakta_db(tmp_path / 'bakta-db' / 'db-light')
        code, out = _check(tmp_path, '--bakta-db', on_host, '--bakta-db-staged', 'bakta-db/db-light')
        assert code == 0, out
        assert f'✓ Bakta DB: {on_host} (light 6.0)' in out

    def test_one_that_is_not_there_stops_the_run(self, tmp_path):
        on_host = _unseen(tmp_path, 'db-light')
        code, out = _check(tmp_path, '--bakta-db', on_host)
        assert code == 2
        assert f'Bakta DB not found: {on_host}' in out

    def test_a_folder_that_is_not_a_bakta_database_stops_the_run(self, tmp_path):
        (tmp_path / 'genomes').mkdir()
        code, out = _check(tmp_path, '--bakta-db', str(tmp_path / 'genomes'))
        assert code == 2
        assert 'has no version.json, so it is not a Bakta database' in out

    def test_the_folder_it_was_downloaded_to_points_at_the_one_inside(self, tmp_path):
        _bakta_db(tmp_path / 'bakta_db' / 'db-light')
        code, out = _check(tmp_path, '--bakta-db', str(tmp_path / 'bakta_db'))
        assert code == 2
        assert f"point --bakta_db at {tmp_path / 'bakta_db' / 'db-light'}" in out

    def test_one_without_bakta_db_stops_the_run(self, tmp_path):
        db = _bakta_db(tmp_path / 'db-light', files=('pfam.h3f',))
        code, out = _check(tmp_path, '--bakta-db', str(db))
        assert code == 2
        assert f'Bakta DB: {db} has no bakta.db' in out


class TestEggnogDb:

    def test_one_not_built_yet_is_no_error(self, tmp_path):
        # EGGNOG_DB builds it, so a first run must get past validation
        on_host = _unseen(tmp_path, 'eggnog_db')
        code, out = _check(tmp_path, '--eggnog-db', on_host)
        assert code == 0, out
        assert f'✓ eggNOG DB: {on_host} (not built yet: EGGNOG_DB will download and build it' in out

    def test_a_built_one_is_read_where_staged(self, tmp_path):
        on_host = _unseen(tmp_path, 'eggnog_db')
        _eggnog_db(tmp_path / 'eggnog-db' / 'eggnog_db', dmnd_bytes=BUILT)
        code, out = _check(tmp_path, '--eggnog-db', on_host, '--eggnog-db-staged', 'eggnog-db/eggnog_db')
        assert code == 0, out
        assert f'✓ eggNOG DB: {on_host}\n' in out
        assert '⚠' not in out

    def test_eggnog_db_ok_with_nothing_built_stops_the_run(self, tmp_path):
        # What EGGNOG_DB leaves on the host under Docker: Nextflow will never
        # run it again, and eggNOG-mapper would find no database
        db = _eggnog_db(tmp_path / 'eggnog_db', annotations=False)
        code, out = _check(tmp_path, '--eggnog-db', str(db))
        assert code == 2
        assert ('has eggnog_db.ok, so EGGNOG_DB will not run again, but there is no '
                'bacteria.dmnd; delete eggnog_db.ok to have it built again') in out

    def test_a_short_bacteria_dmnd_is_built_again(self, tmp_path):
        db = _eggnog_db(tmp_path / 'eggnog_db', dmnd_bytes=1000, built=False)
        code, out = _check(tmp_path, '--eggnog-db', str(db))
        assert code == 0, out
        assert ('bacteria.dmnd in {} is only 1,000 bytes; EGGNOG_DB will take it for a '
                'truncated download and build it again').format(db) in out

    def test_a_short_bacteria_dmnd_under_eggnog_db_ok_stops_the_run(self, tmp_path):
        db = _eggnog_db(tmp_path / 'eggnog_db', dmnd_bytes=1000)
        code, out = _check(tmp_path, '--eggnog-db', str(db))
        assert code == 2
        assert 'but bacteria.dmnd is only 1,000 bytes; delete eggnog_db.ok to have it built again' in out

    def test_no_eggnog_db_stops_the_run(self, tmp_path):
        db = _eggnog_db(tmp_path / 'eggnog_db', dmnd_bytes=BUILT, annotations=False, built=False)
        code, out = _check(tmp_path, '--eggnog-db', str(db))
        assert code == 2
        assert 'has no eggnog.db, which eggNOG-mapper needs' in out
        assert 'delete bacteria.dmnd to have it built again' in out

    def test_the_size_it_takes_for_truncated_is_eggnog_dbs(self):
        from validate_inputs import EGGNOG_DMND_MIN_BYTES
        text = (ROOT / 'modules' / 'local' / 'eggnog.nf').read_text()
        assert re.findall(r'"\\\$DB_SIZE" -lt (\d+)', text) == [str(EGGNOG_DMND_MIN_BYTES)]


class TestModuleDatabases:

    def test_module_stages_each_database_in_a_folder_of_its_own(self):
        stages = _module()[1]
        assert stages.get('bakta_db') and stages.get('eggnog_db'), \
            'VALIDATE_INPUTS does not stage the databases'
        names = {_staged_name(stages[n], 1, 'db') for n in ('bakta_db', 'eggnog_db')}
        assert len(names) == 2 and not any(n == 'db' for n in names)

    def test_databases_on_the_host_are_read_where_staged(self, tmp_path):
        sources = [_write(tmp_path / 'store' / 'G1.fasta', _genome())]
        rows = [('G1', _unseen(tmp_path, 'G1.fasta'))]
        bakta = _unseen(tmp_path, 'db-light')
        eggnog = _unseen(tmp_path, 'eggnog_db')
        databases = {'bakta': (bakta, _bakta_db(tmp_path / 'store' / 'db-light')),
                     'eggnog': (eggnog, _eggnog_db(tmp_path / 'store' / 'eggnog_db', dmnd_bytes=BUILT))}
        code, log, _ = _run_module(tmp_path, rows, sources, databases)
        assert f'✓ Bakta DB: {bakta} (light 6.0)' in log
        assert f'✓ eggNOG DB: {eggnog}\n' in log
        assert code == 0

    def test_strict_validation_stops_on_a_warning(self, tmp_path):
        # Under 500 kb a genome passes with a warning
        sources = [_write(tmp_path / 'store' / 'G1.fasta', _genome(110_000))]
        rows = [('G1', _unseen(tmp_path, 'G1.fasta'))]
        code, log, _ = _run_module(tmp_path, rows, sources)
        assert code == 0
        code, log, _ = _run_module(tmp_path / 'strict', rows, sources, strict=True)
        assert code == 1
        assert 'FAILED (strict mode): 1 warning(s)' in log
