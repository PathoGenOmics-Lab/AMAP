"""Bakta's jobs ask for SLURM's QoS only under SLURM.

nextflow.config gave BAKTA_ANNOTATOR `clusterOptions = '--qos=medium'` for
every run. It is an option of SLURM's sbatch, so under -profile sge, pbs or lsf
it went to a submit command that has no such option, and as a process selector
it took the place of the `-l h_vmem` the sge profile gives every task. The QoS
now lives in the slurm profile.

Without Nextflow, as in CI, the tests read nextflow.config; with Nextflow on
PATH, `nextflow config` resolves each profile, in both orders with a container
and under both of its parsers.
"""
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
CONFIG = (ROOT / 'nextflow.config').read_text()
QOS = "'--qos=medium'"
SGE_MEMORY = '{ "-l h_vmem=${task.memory.toGiga()}G" }'


def _block(text, header):
    """The body of the first `header {` block in text, braces matched."""
    m = re.search(rf'^\s*{header}\s*\{{', text, re.M)
    assert m, f'no {header} block'
    depth, i = 1, m.end()
    while depth:
        depth += {'{': 1, '}': -1}.get(text[i], 0)
        i += 1
    return text[m.end():i - 1], (m.start(), i)


def _without(text, span):
    return text[:span[0]] + text[span[1]:]


def test_no_cluster_option_outside_the_profiles():
    _, span = _block(CONFIG, 'profiles')
    assert 'clusterOptions' not in _without(CONFIG, span)


def test_bakta_asks_for_the_qos_under_slurm():
    profiles, _ = _block(CONFIG, 'profiles')
    slurm, _ = _block(profiles, 'slurm')
    bakta, _ = _block(slurm, "withName: 'BAKTA_ANNOTATOR'")
    assert re.search(rf'^\s*clusterOptions\s*=\s*{re.escape(QOS)}\s*$', bakta, re.M)


@pytest.mark.parametrize('name', ['sge', 'pbs', 'lsf'])
def test_other_schedulers_get_no_qos(name):
    profiles, _ = _block(CONFIG, 'profiles')
    body, _ = _block(profiles, name)
    assert '--qos' not in body and 'BAKTA_ANNOTATOR' not in body


# ── nextflow config ──

NEXTFLOW = shutil.which('nextflow')
needs_nextflow = pytest.mark.skipif(NEXTFLOW is None, reason='Nextflow is not on PATH; the tests above read the config instead')

# -profile, then Bakta's own clusterOptions and every task's
CASES = [
    ('slurm', QOS, None),
    ('singularity,slurm', QOS, None),
    ('slurm,singularity', QOS, None),
    ('sge', None, SGE_MEMORY),
    ('docker,sge', None, SGE_MEMORY),
    ('pbs', None, None),
    ('lsf', None, None),
    ('docker', None, None),
]


def _flat(profiles, parser):
    with tempfile.TemporaryDirectory() as cwd:
        result = subprocess.run(
            [NEXTFLOW, '-log', os.path.join(cwd, 'nextflow.log'), '-q', 'config', '-flat',
             '-profile', profiles, str(ROOT)],
            cwd=cwd, env=dict(os.environ, NXF_SYNTAX_PARSER=parser),
            capture_output=True, text=True, timeout=600,
        )
    assert result.returncode == 0, result.stdout + result.stderr
    bakta = task = None
    for line in result.stdout.splitlines():
        name, sep, value = line.partition(' = ')
        if not sep:
            continue
        if name == "process.'withName:BAKTA_ANNOTATOR'.clusterOptions":
            bakta = value.strip()
        elif name == 'process.clusterOptions':
            task = value.strip()
    return bakta, task


@needs_nextflow
@pytest.mark.parametrize('parser', ['v1', 'v2'])
@pytest.mark.parametrize('profiles,bakta,task', CASES)
def test_cluster_options_from_nextflow_config(profiles, bakta, task, parser):
    assert _flat(profiles, parser) == (bakta, task)
