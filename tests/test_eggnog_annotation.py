#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANNOT_EGGNOGMAP_BATCH finds eggNOG's database whatever its folder is called.

Nextflow stages a folder that is a path input under the folder's own name, and
the script asked emapper.py for eggnog_db/bacteria.dmnd and --data_dir
eggnog_db. So only a --eggnog_db folder called eggnog_db worked: with the
/shared/databases/eggnog of docs/deployment/hpc.md or the tutorial's
~/databases/eggnog, EGGNOG_DB built the database and the annotation stopped at
"DIAMOND database eggnog_db/bacteria.dmnd not present", with or without a
container. The tag and the first log line also counted the FAA files with
size(), which for a single sample is the size of its FAA in bytes
("batch_531_samples" for a FAA of 531 bytes).

The tests run the process's script the way Nextflow renders it, in a task
directory with each input staged as Nextflow stages it, a link under the name
its declaration gives it or else under its own, and a stand-in for emapper.py
that looks for the database as eggNOG-mapper 2.1.13 does.
"""

import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from tests.test_eggnog_db import BUILT, _sparse
from tests.test_validate_inputs import _render, _staged_name

ROOT = Path(__file__).resolve().parent.parent
MODULE = ROOT / 'modules' / 'local' / 'eggnog.nf'

# What EGGNOG_DB leaves in --eggnog_db, its storeDir: the bacteria-only index,
# the annotation databases it downloads, and the marker storeDir keeps
DATABASE = ('bacteria.dmnd', 'eggnog.db', 'eggnog.taxa.db', 'eggnog.taxa.db.traverse.pkl', 'eggnog_db.ok')

# eggNOG-mapper 2.1.13's checks, with its messages: the DIAMOND database at
# --dmnd_db and eggnog.db in --data_dir have to exist. Then one hit per query,
# in the columns of a real .emapper.annotations file.
EMAPPER = '''#!@PYTHON@
import os, sys
args = sys.argv[1:]
opt = lambda k: args[args.index(k) + 1]
with open(@LOG@, 'a') as f:
    f.write(' '.join(args) + '\\n')
dmnd, data_dir = opt('--dmnd_db'), opt('--data_dir')
if not os.path.exists(dmnd):
    sys.exit('DIAMOND database %s not present. Use download_eggnog_data.py to fetch it' % dmnd)
if not os.path.exists(os.path.join(data_dir, 'eggnog.db')):
    sys.exit('Annotation database data/eggnog.db not present. Use download_eggnog_data.py to fetch it')
cols = ['#query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 'max_annot_lvl', 'COG_category',
        'Description', 'Preferred_name', 'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module',
        'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs']
with open(opt('-i')) as fq, open(opt('-o') + '.emapper.annotations', 'w') as fo:
    fo.write('## emapper-2.1.13\\n' + '\\t'.join(cols) + '\\n')
    for line in fq:
        if line.startswith('>'):
            fo.write('\\t'.join([line[1:].split()[0], '83332.Rv0001', '1e-100', '500.0'] + ['-'] * (len(cols) - 4)) + '\\n')
'''

# Every sample's proteins as FIX_CDS writes them, <sampleId>_<gene>;<strand>;<start>-<end>:
# dnaA is the same in every sample, so eggNOG-mapper sees it once
DNAA = 'MTDDPGSGFTTVWNAVVSELNGDPKVDDGPSSDANLSAPLTPQQRAWLNLVQPLTIVEGFALLSVPSSFVQNEIERHLRAPITDALSRRLG'


def _grep_reads_perl_regexps():
    done = subprocess.run(['grep', '-oP', r'\d+'], input='2 GB\n', capture_output=True, text=True)
    return done.stdout.strip() == '2'


# The script works its memory out with grep -P, which BSD grep (macOS) lacks
pytestmark = pytest.mark.skipif(not _grep_reads_perl_regexps(),
                                reason="the script's MEM_GB line needs GNU grep -P; run the tests in the container")


def _annot():
    """ANNOT_EGGNOGMAP_BATCH's tag, the options of each of its inputs, and its
    script block as Nextflow strips it."""
    text = MODULE.read_text()
    process = re.search(r'^process ANNOT_EGGNOGMAP_BATCH \{\n(.*?)^\}', text, re.S | re.M).group(1)
    tag = re.search(r'^\s*tag "(.*)"\s*$', process, re.M).group(1)
    block = re.search(r'^\s*input:\n(.*?)^\s*output:', process, re.S | re.M).group(1)
    inputs = {name: dict(re.findall(r"(\w+):\s*'([^']*)'", options))
              for name, options in re.findall(r'^\s*(?:path|val)\s+(\w+)\s*,?([^\n]*?)\s*(?://.*)?$', block, re.M)}
    script = re.search(r'^\s*script:\n\s*"""\n(.*?)\n\s*"""\n', process, re.S | re.M).group(1)
    return tag, inputs, textwrap.dedent(script) + '\n'


def _stage(work, options, sources):
    """Link each source into the task directory as Nextflow stages a path
    input: under its stageAs, or else under the file's own name. Returns the
    names, blank-separated as the script gets them."""
    names = []
    for index, source in enumerate(sources, 1):
        name = _staged_name(options['stageAs'], index, source.name) if 'stageAs' in options else source.name
        (work / name).parent.mkdir(parents=True, exist_ok=True)
        (work / name).symlink_to(source)
        names.append(name)
    return ' '.join(names)


def _size(options, files):
    """What size() gives on a path input. Nextflow hands one file over as
    that file, whose size() is its size in bytes, unless the input's arity
    takes more than one: then, as for several files, a list, which size()
    counts."""
    upper = options.get('arity', '1').split('..')[-1]
    if len(files) == 1 and upper in ('0', '1'):
        return str(files[0].stat().st_size)
    return str(len(files))


def _run(tmp_path, folder='eggnog_db', samples=('L4_S1', 'L4_S2')):
    """Run the script for `samples` against the database EGGNOG_DB built in a
    --eggnog_db folder called `folder`. (exit code, stdout + stderr, the tag,
    the folder on the host, the task directory, emapper.py's calls)."""
    tag, inputs, script = _annot()
    host_db = tmp_path / 'shared' / 'databases' / folder
    host_db.mkdir(parents=True)
    for name in DATABASE:
        (host_db / name).write_text('')
    _sparse(host_db / 'bacteria.dmnd', BUILT)

    faas = []
    for n, sid in enumerate(samples):
        faa = tmp_path / 'run' / 'results' / sid / f'{sid}.AMAP.faa'
        faa.parent.mkdir(parents=True)
        faa.write_text(f'>{sid}_Rv0001;+;0-1524\n{DNAA}\n'
                       f'>{sid}_Rv0002;+;2051-3260\nMAEPV{"KLMNPQRSTV"[n]}LRGSN\n')
        faas.append(faa)

    work = tmp_path / 'run' / 'work' / 'ab' / 'cdef0123'
    work.mkdir(parents=True)
    tools = tmp_path / 'tools'
    tools.mkdir()
    log = tmp_path / 'emapper.calls'
    (tools / 'emapper.py').write_text(EMAPPER.replace('@PYTHON@', sys.executable).replace('@LOG@', repr(str(log))))
    (tools / 'emapper.py').chmod(0o755)

    values = {
        'faa_files': _stage(work, inputs['faa_files'], faas),
        'faa_files.size()': _size(inputs['faa_files'], faas),
        'sample_ids.size()': str(len(samples)),
        'eggnog_db': _stage(work, inputs['eggnog_db'], [host_db]),
        'task.memory': '2 GB',
        'task.cpus': '1',
    }
    (work / '.command.sh').write_text(_render(script, values))
    env = dict(os.environ, PATH=f'{tools}{os.pathsep}{os.environ["PATH"]}')
    # nextflow.config runs every task with bash -euo pipefail
    done = subprocess.run(['bash', '-euo', 'pipefail', '.command.sh'], cwd=work, env=env,
                          capture_output=True, text=True)
    calls = log.read_text().splitlines() if log.exists() else []
    return done.returncode, done.stdout + done.stderr, _render(tag, values), host_db, work, calls


def _annotated(work, sid):
    """The queries in a sample's split annotations."""
    lines = (work / 'split' / f'{sid}.emapper.annotations').read_text().splitlines()
    return [line.split('\t')[0] for line in lines if not line.startswith('#')]


class TestAnnotEggnogmapBatch:

    @pytest.mark.parametrize('folder', ['eggnog_db', 'eggnog', 'eggnog-5.0.2'])
    def test_the_database_is_found_whatever_its_folder_is_called(self, tmp_path, folder):
        # 'eggnog' as in docs/deployment/hpc.md's --eggnog_db /shared/databases/eggnog
        code, out, _, _, work, calls = _run(tmp_path, folder=folder)
        assert 'not present' not in out, out
        assert code == 0, out
        assert len(calls) == 1, 'eggNOG-mapper ran more than once for one batch'
        # Every protein comes back under its own sample, dnaA too, which
        # eggNOG-mapper saw only once
        for sid in ('L4_S1', 'L4_S2'):
            assert _annotated(work, sid) == [f'{sid}_Rv0001;+;0-1524', f'{sid}_Rv0002;+;2051-3260']

    def test_a_folder_named_like_the_one_the_script_writes_in_is_left_as_it_was(self, tmp_path):
        # The script writes each sample's annotations in split/: a database
        # staged under its own name, split, would take them in
        code, out, _, host_db, work, _ = _run(tmp_path, folder='split')
        assert code == 0, out
        assert sorted(p.name for p in host_db.iterdir()) == sorted(DATABASE)
        assert not (work / 'split').is_symlink()
        assert _annotated(work, 'L4_S1') == ['L4_S1_Rv0001;+;0-1524', 'L4_S1_Rv0002;+;2051-3260']

    @pytest.mark.parametrize('n', [1, 3])
    def test_the_tag_and_the_log_count_the_samples(self, tmp_path, n):
        code, out, tag, *_ = _run(tmp_path, samples=[f'L4_S{i}' for i in range(1, n + 1)])
        assert code == 0, out
        assert tag == f'batch_{n}_samples'
        assert f'Preparing batch annotation for {n} samples...' in out
