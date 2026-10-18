#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EGGNOG_DB builds eggNOG's database in the folder main.nf hands it.

Its script used to work the folder out for itself, from --eggnog_db and the
launch directory. Without a container that is the folder on the host, but a
container sees only the folders of a task's inputs, and EGGNOG_DB had none.
Under -profile docker, unless the launch directory happened to share a parent
folder with the pipeline's bin/ (Nextflow mounts the two together), the
database was downloaded and built inside the container and thrown away with
it, and storeDir moved only eggnog_db.ok to the host. Every later run skipped
EGGNOG_DB on that marker, and ANNOT_EGGNOGMAP_BATCH found no bacteria.dmnd. A
database already on the host was not seen either, and was built again. main.nf
now makes the folder on the host and hands it over as an input, which Nextflow
mounts read-write.

The tests run EGGNOG_DB's script the way Nextflow renders it, in a task
directory with the folder staged as Nextflow stages an input, a link to the
one on the host, and stand-ins for eggNOG-mapper's download and build. A
folder the script works out for itself is rendered as one inside the
container, which the host never sees.
"""

import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from tests.test_validate_inputs import _render

ROOT = Path(__file__).resolve().parent.parent
MODULE = ROOT / 'modules' / 'local' / 'eggnog.nf'

# Just over the 100 MB under which EGGNOG_DB takes bacteria.dmnd for a
# truncated download (validate_inputs.EGGNOG_DMND_MIN_BYTES)
BUILT = 100_000_001

# Stand-ins for eggNOG-mapper's download and build: each writes the files the
# real one leaves in --data_dir, empty, and logs its call outside the
# database, so that a test can tell whether it ran.
STAND_INS = {
    'download_eggnog_data.py': '''#!{python}
import os, sys
args = sys.argv[1:]
data_dir = args[args.index('--data_dir') + 1]
with open({log!r}, 'a') as f:
    f.write('download\\n')
for name in ('eggnog.db', 'eggnog.taxa.db', 'eggnog.taxa.db.traverse.pkl'):
    open(os.path.join(data_dir, name), 'w').close()
''',
    'create_dbs.py': '''#!{python}
import os, sys
args = sys.argv[1:]
data_dir = args[args.index('--data_dir') + 1]
with open({log!r}, 'a') as f:
    f.write('build\\n')
for name in ('e5.proteomes.faa', 'e5.taxid_info.tsv'):
    open(os.path.join(data_dir, name), 'w').close()
with open(os.path.join(data_dir, args[args.index('--dbname') + 1] + '.dmnd'), 'wb') as f:
    f.truncate({size})
''',
}


def _eggnog_db():
    """EGGNOG_DB's script block as Nextflow strips it, and the names of its
    path inputs."""
    text = MODULE.read_text()
    process = re.search(r'^process EGGNOG_DB \{\n(.*?)^\}', text, re.S | re.M).group(1)
    script = re.search(r'^    script:\n    """\n(.*?)\n    """\n', process, re.S | re.M)
    block = re.search(r'^    input:\n(.*?)^    output:', process, re.S | re.M)
    inputs = re.findall(r'^    path (\w+)\s*$', block.group(1), re.M) if block else []
    return textwrap.dedent(script.group(1)) + '\n', inputs


def _sparse(path, size):
    """A file of that many bytes that costs nothing on disk."""
    with open(path, 'wb') as f:
        f.truncate(size)


def _run(tmp_path, name='eggnog_db', dmnd_bytes=None, annotations=False):
    """Run EGGNOG_DB's script in a task directory under the launch folder's
    work/, with --eggnog_db given as `name`, relative to the launch folder,
    whose folder main.nf has made on the host: empty, or holding a
    bacteria.dmnd of `dmnd_bytes` and eggnog.db if `annotations`.
    (exit code, stdout + stderr, the folder on the host, the task directory,
    the stand-ins called)."""
    script, inputs = _eggnog_db()
    launch = tmp_path / 'host' / 'run'
    host_db = launch / name
    host_db.mkdir(parents=True)
    if dmnd_bytes is not None:
        _sparse(host_db / 'bacteria.dmnd', dmnd_bytes)
    if annotations:
        (host_db / 'eggnog.db').write_text('')

    work = launch / 'work' / 'ab' / 'cdef0123'
    work.mkdir(parents=True)
    # A path input is staged as a link named after what it points to
    for _ in inputs:
        (work / host_db.name).symlink_to(host_db)

    tools = tmp_path / 'tools'
    tools.mkdir()
    log = tmp_path / 'calls.log'
    for tool, stand_in in STAND_INS.items():
        (tools / tool).write_text(stand_in.format(python=sys.executable, log=str(log), size=BUILT))
        (tools / tool).chmod(0o755)

    rendered = _render(script, {
        'eggnog_db': host_db.name,
        # What the script used to work the folder out from. Inside a
        # container, the launch folder Nextflow does not mount is a folder of
        # the container's own.
        'params.eggnog_db': name,
        'launchDir': str(tmp_path / 'container' / 'run'),
    })
    (work / '.command.sh').write_text(rendered)
    env = dict(os.environ, PATH=f'{tools}{os.pathsep}{os.environ["PATH"]}')
    # nextflow.config runs every task with bash -euo pipefail
    done = subprocess.run(['bash', '-euo', 'pipefail', '.command.sh'], cwd=work, env=env,
                          capture_output=True, text=True)
    calls = log.read_text().split() if log.exists() else []
    return done.returncode, done.stdout + done.stderr, host_db, work, calls


class TestEggnogDb:

    @pytest.mark.parametrize('name', ['eggnog_db', 'eggnog'])
    def test_the_database_is_built_in_the_folder_handed_over(self, tmp_path, name):
        # 'eggnog' as in docs/deployment/hpc.md's --eggnog_db /shared/databases/eggnog
        code, out, host_db, work, calls = _run(tmp_path, name=name)
        assert code == 0, out
        assert calls == ['download', 'build']
        assert (host_db / 'bacteria.dmnd').is_file(), \
            f'bacteria.dmnd was not built in the folder EGGNOG_DB was handed: {out}'
        assert (host_db / 'bacteria.dmnd').stat().st_size == BUILT
        assert (host_db / 'eggnog.db').is_file()
        # The build-only files are gone, and the marker storeDir moves is made
        assert not (host_db / 'e5.proteomes.faa').exists()
        assert not (host_db / 'e5.taxid_info.tsv').exists()
        assert (work / 'eggnog_db.ok').is_file()

    def test_a_database_already_on_the_host_is_used_as_it_is(self, tmp_path):
        # Built by an earlier run without a container, say, which leaves
        # eggnog_db.ok there too; without it, Nextflow runs EGGNOG_DB
        code, out, host_db, work, calls = _run(tmp_path, dmnd_bytes=BUILT, annotations=True)
        assert code == 0, out
        assert calls == [], f'EGGNOG_DB downloaded and built a database the host has: {out}'
        assert f'bacteria.dmnd found ({BUILT} bytes)' in out
        assert (host_db / 'bacteria.dmnd').stat().st_size == BUILT
        assert (work / 'eggnog_db.ok').is_file()

    def test_a_short_bacteria_dmnd_on_the_host_is_built_again(self, tmp_path):
        code, out, host_db, work, calls = _run(tmp_path, dmnd_bytes=1000)
        assert code == 0, out
        assert 'bacteria.dmnd is only 1000 bytes' in out
        assert calls == ['download', 'build']
        assert (host_db / 'bacteria.dmnd').stat().st_size == BUILT

    def test_main_nf_makes_the_folder_and_both_processes_get_it(self):
        text = (ROOT / 'main.nf').read_text()
        call = re.search(r'\bEGGNOG_DB\(channel\.value\((\w+)\)\)', text)
        assert call, 'main.nf does not hand EGGNOG_DB the folder to build the database in'
        folder = call.group(1)
        assert re.search(rf'\b{folder}\??\.mkdirs\(\)', text[:call.start()]), \
            f'main.nf hands EGGNOG_DB {folder} without making it on the host first'
        assert re.search(rf'EGGNOG_DB\.out\.db_path\.map \{{ {folder} \}}', text), \
            f'ANNOT_EGGNOGMAP_BATCH does not get {folder}, the folder EGGNOG_DB built in'
