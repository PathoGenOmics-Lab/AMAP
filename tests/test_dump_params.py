"""DUMP_PARAMS writes down the parameters of the run it is part of, as YAML.

Two defects lived here. Nextflow strips a script's indentation once the script
is filled in, and the parameters went in as one string indented less than the
rest, so Nextflow stripped less: all but the first parameter fell outside
`params:`, and the heredoc's delimiter kept two spaces, so the file ended in a
line `  YAML`. And the script reads `params` as a whole, which Nextflow leaves
out of the task's hash, so with -resume every later run got the first run's
file back. The tests with Nextflow run the process itself; CI has no Nextflow,
so the guard on the cache runs everywhere.
"""
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
MODULE = ROOT / 'modules' / 'local' / 'versions.nf'

NEXTFLOW = shutil.which('nextflow')
needs_nextflow = pytest.mark.skipif(NEXTFLOW is None, reason='Nextflow is not on PATH')

# Values that broke the old file or the shell that writes it
TRICKY = {
    'multi': 'line one\nline two',
    'quote': 'a "b" c',
    'dollar': '$HOME `id`',
    'back': 'C:\\path',
    'colon': 'a: b # c',
}


def _process(name):
    body = re.search(rf'^process {name} \{{\n(.*?)^\}}', MODULE.read_text(), re.S | re.M)
    assert body, f'no process {name} in {MODULE}'
    return body.group(1)


def test_never_from_the_cache():
    # -resume has to write the parameters of the run, not bring back the first one's
    assert re.search(r'^\s*cache false\s*$', _process('DUMP_PARAMS'), re.M)


def _read(path):
    """The header lines and the parameters of a params.yml, each value read as
    the JSON string it is written as."""
    header, params, in_params = [], {}, False
    for line in path.read_text().splitlines():
        if line == 'params:':
            in_params = True
        elif not in_params:
            header.append(line)
        else:
            m = re.fullmatch(r'  ([\w.-]+): (".*")', line)
            assert m, f'not a parameter under params: {line!r}'
            params[m.group(1)] = json.loads(m.group(2))
    return header, params


def _run(work, args, parser):
    result = subprocess.run(
        [NEXTFLOW, '-q', '-log', str(work / 'nextflow.log'), 'run', 'main.nf',
         '-params-file', 'p.json', '--outdir', 'out'] + args,
        cwd=work, env=dict(os.environ, NXF_SYNTAX_PARSER=parser),
        capture_output=True, text=True, timeout=600,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return _read(work / 'out' / 'pipeline_info' / 'params.yml')


@needs_nextflow
@pytest.mark.parametrize('parser', ['v1', 'v2'])
def test_params_yml_is_the_run_s(tmp_path, parser):
    module = (ROOT / 'modules' / 'local' / 'versions').as_posix()
    (tmp_path / 'main.nf').write_text(f"include {{ DUMP_PARAMS }} from '{module}'\nworkflow {{ DUMP_PARAMS() }}\n")
    (tmp_path / 'p.json').write_text(json.dumps(TRICKY))

    header, params = _run(tmp_path, ['--num', '8'], parser)
    assert header[0] == '# AMAP Pipeline Parameters'
    assert all(line.startswith('# ') for line in header)
    assert params == dict(TRICKY, outdir='out', num='8')

    # The same run again with one parameter changed
    _, params = _run(tmp_path, ['--num', '9', '-resume'], parser)
    assert params == dict(TRICKY, outdir='out', num='9')
