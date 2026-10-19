"""environment.yml pins the versions the container ships.

environment.yml says its pinned versions match the image paururo/amap:latest,
but it pinned eggnog-mapper 2.1.12 while docker/Dockerfile installs 2.1.13,
which is what the image has. Every tool environment.yml pins has to be in the
Dockerfile at that version, whichever way the Dockerfile installs it: from
conda (`bakta=1.12.0`), from a git tag (`--branch v0.18`) or from a release
(`/download/v1.1/`).
"""
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
ENVIRONMENT = (ROOT / 'environment.yml').read_text()
DOCKERFILE = (ROOT / 'docker' / 'Dockerfile').read_text()


def _pins():
    """name -> version for every exact pin in environment.yml's conda list."""
    pins = {}
    for line in ENVIRONMENT.splitlines():
        m = re.match(r'\s{2}- ([A-Za-z0-9_.-]+)=([0-9][\w.]*)\s*(#.*)?$', line)
        if m:
            pins[m.group(1)] = m.group(2)
    return pins


def test_reads_the_pins():
    pins = _pins()
    assert pins['bakta'] == '1.12.0' and 'miniprot' in pins and 'pangene' in pins


def _in_dockerfile(name, version):
    """The versions the Dockerfile gives this tool, however it installs it."""
    found = set(re.findall(rf'(?<![\w-]){re.escape(name)}=([0-9][\w.]*)', DOCKERFILE))
    for line in DOCKERFILE.splitlines():
        if re.search(rf'/{re.escape(name)}\b', line):
            found |= set(re.findall(r'(?:--branch |/download/)v([0-9][\w.]*)', line))
    return found


@pytest.mark.parametrize('name,version', sorted(_pins().items()))
def test_the_container_has_the_pinned_version(name, version):
    assert _in_dockerfile(name, version) == {version}, \
        f'environment.yml pins {name}={version}, docker/Dockerfile installs {_in_dockerfile(name, version) or "nothing"}'
