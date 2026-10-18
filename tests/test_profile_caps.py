#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The caps on what any one task may ask for do not depend on the order of
-profile.

Nextflow 24.04 applies the profiles a run names in the order nextflow.config
defines them. From 25.10, whose strict parser is the default, it applies them
in the order -profile gives them, and where two set the same thing the last one
wins. The container profiles used to set a laptop's caps, the machine's CPUs and
16 GB, so `-profile slurm,singularity` held a cluster run to them where
`-profile singularity,slurm` got the scheduler's 48 CPUs and 128 GB, and
`-profile test,docker` gave the test 16 GB instead of its 8. Now no profile sets
a cap: each says what kind of run it is, and the caps are set once, after the
profiles block.

Without Nextflow, as in CI, the tests read nextflow.config and assign its params
the way Nextflow does: the params block, then each profile named, where the
profiles block stands, then what follows the block. With Nextflow on PATH,
`nextflow config` also resolves each case in both orders and under both of its
parsers, and so does `nextflow run --help` for a few, since only a run makes
Nextflow 24.04 refuse a config that reads a param it never set.
"""

import itertools
import os
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
CONFIG = (ROOT / 'nextflow.config').read_text()
CAPS = ('max_cpus', 'max_memory', 'max_time')
CONTAINERS = ('standard', 'docker', 'singularity')
OTHERS = ('slurm', 'sge', 'pbs', 'lsf', 'test', 'highmem')
# What Runtime.runtime.availableProcessors() gives, whatever the machine has.
MACHINE = 'the machine'

# -profile, then the CPUs, memory and time it caps a task at. Each is also
# checked with its profiles in reverse order.
CASES = [
    ('docker', MACHINE, '16.GB', '72.h'),
    ('singularity', MACHINE, '16.GB', '72.h'),
    ('standard', MACHINE, '16.GB', '72.h'),
    ('singularity,slurm', 48, '128.GB', '72.h'),
    ('docker,sge', 48, '128.GB', '72.h'),
    ('singularity,pbs', 48, '128.GB', '72.h'),
    ('singularity,lsf', 48, '128.GB', '72.h'),
    ('docker,test', 4, '8.GB', '1.h'),
    ('slurm,test', 4, '8.GB', '1.h'),
    ('slurm,highmem', 48, '256.GB', '120.h'),
    ('singularity,slurm,highmem', 48, '256.GB', '120.h'),
]
PARSERS = ('v1', 'v2')


def _orders(profiles):
    names = profiles.split(',')
    return [profiles] if len(names) == 1 else [profiles, ','.join(reversed(names))]


# ── nextflow.config, read the way Nextflow applies it ──

def _sections():
    """The params block, each profile's body by name in the order the file
    defines them, and what follows the profiles block."""
    params = re.search(r'^params \{\n(.*?)^\}', CONFIG, re.M | re.S).group(1)
    start = re.search(r'^profiles \{\n', CONFIG, re.M).end()
    end = start + re.search(r'^\}', CONFIG[start:], re.M).end()
    profiles = dict(re.findall(r'^    (\w+) \{\n(.*?)^    \}', CONFIG[start:end], re.M | re.S))
    return params, profiles, CONFIG[end:]


PARAMS_BLOCK, PROFILES, AFTER = _sections()


def _assignments(text, prefix):
    """(name, expression) of each param a stretch of the config assigns:
    `name = x` inside the params block, `params.name = x` elsewhere."""
    head = r'params\.' if prefix else ''
    return re.findall(r'^[ \t]*' + head + r'(\w+)[ \t]*=[ \t]*(.+?)[ \t]*(?://.*)?$', text, re.M)


def _split(expr, mark):
    """expr cut at the first `mark` outside quotes and brackets."""
    depth, quote = 0, None
    for i, ch in enumerate(expr):
        if quote:
            quote = None if ch == quote else quote
        elif ch in '\'"':
            quote = ch
        elif ch in '([{':
            depth += 1
        elif ch in ')]}':
            depth -= 1
        elif ch == mark and depth == 0:
            return expr[:i], expr[i + 1:]
    return expr, None


def _value(expr, params):
    """What Nextflow makes of an expression, as far as the caps need it:
    literals, params.NAME, the machine's CPUs, `a ? b : c` and `a ?: b`.
    Anything else is Unread, which fails a test only if a cap depends on it."""
    expr = expr.strip()
    condition, rest = _split(expr, '?')
    if rest is not None:
        first = _value(condition, params)
        if not _known(first):
            return Unread(expr)
        if rest.startswith(':'):
            return first or _value(rest[1:], params)
        then, otherwise = _split(rest, ':')
        if otherwise is None:
            return Unread(expr)
        return _value(then if first else otherwise, params)
    if re.fullmatch(r'-?\d+', expr):
        return int(expr)
    if re.fullmatch(r"'[^']*'", expr):
        return expr[1:-1]
    if expr in ('true', 'false', 'null'):
        return {'true': True, 'false': False, 'null': None}[expr]
    if expr == 'Runtime.runtime.availableProcessors()':
        return MACHINE
    param = re.fullmatch(r'params\.(\w+)', expr)
    if param:
        name = param.group(1)
        assert name in params, (
            f'nextflow.config reads params.{name} before giving it a value, and a run on Nextflow 24.04 '
            'stops there ("Unknown config attribute"); give it one in the params block'
        )
        return params[name]
    return Unread(expr)


class Unread(str):
    """An expression these tests do not read, a timestamp for instance."""


def _known(value):
    return not isinstance(value, Unread)


def read_caps(profiles, parser):
    """The caps nextflow.config gives a run with this -profile. v1 (Nextflow
    24.04, or later ones with NXF_SYNTAX_PARSER=v1) applies the profiles named
    in the order the file defines them, v2 (the strict parser, the default from
    25.10) in the order given; both where the profiles block stands. With no
    -profile, Nextflow applies `standard`."""
    names = profiles.split(',') if profiles else ['standard']
    unknown = [name for name in names if name not in PROFILES]
    assert not unknown, f'nextflow.config defines no profile {unknown}'
    if parser == 'v1':
        names = [name for name in PROFILES if name in names]
    params = {}
    stretches = [(PARAMS_BLOCK, False)] + [(PROFILES[name], True) for name in names] + [(AFTER, True)]
    for text, prefix in stretches:
        for name, expr in _assignments(text, prefix):
            params[name] = _value(expr, params)
    caps = tuple(params.get(cap) for cap in CAPS)
    assert all(_known(cap) for cap in caps), f'the caps depend on {caps}, which these tests cannot read'
    return caps


def test_reads_the_profiles():
    assert set(CONTAINERS + OTHERS) <= set(PROFILES)


@pytest.mark.parametrize('name', sorted(PROFILES))
def test_no_profile_sets_a_cap(name):
    """Set in a profile, a cap would be the last profile's from Nextflow 25.10."""
    set_here = sorted({param for param, _ in _assignments(PROFILES[name], True)} & set(CAPS))
    assert not set_here, (
        f'profile {name} sets {set_here}, which the caps set after the profiles block override; '
        'set there, a cap would depend on the order of -profile from Nextflow 25.10. '
        'Give the profile a param of its own and read it after the profiles block, as highmem does'
    )


def test_caps_are_set_after_the_profiles():
    assigned = [name for name, _ in _assignments(AFTER, True)]
    for cap in CAPS:
        assert assigned.count(cap) == 1, f'{cap} is assigned {assigned.count(cap)} times after the profiles block'


@pytest.mark.parametrize('profiles,cpus,memory,time', CASES, ids=[case[0] for case in CASES])
def test_caps_from_the_config(profiles, cpus, memory, time):
    for order in _orders(profiles):
        for parser in PARSERS:
            assert read_caps(order, parser) == (cpus, memory, time), f'-profile {order}, {parser} parser'


def test_no_profile_given_is_the_machine():
    for parser in PARSERS:
        assert read_caps('', parser) == (MACHINE, '16.GB', '72.h')


def test_order_never_changes_the_caps():
    """Every set of up to three of the profiles, in every order."""
    changed = []
    for n in (2, 3):
        for names in itertools.combinations(PROFILES, n):
            seen = {read_caps(','.join(order), parser) for order in itertools.permutations(names) for parser in PARSERS}
            if len(seen) > 1:
                changed.append(f"{','.join(names)}: {sorted(seen, key=str)}")
    assert not changed, 'the order of -profile changes the caps:\n' + '\n'.join(changed)


@pytest.mark.parametrize('container', CONTAINERS)
def test_a_container_profile_never_changes_the_caps(container):
    """Where the tools come from says nothing about how big a task may be."""
    for n in (0, 1, 2):
        for others in itertools.combinations(OTHERS, n):
            for parser in PARSERS:
                assert read_caps(','.join(others + (container,)), parser) == read_caps(','.join(others), parser), (
                    f'-profile {",".join(others + (container,))}, {parser} parser'
                )


# ── nextflow config and nextflow run ──

NEXTFLOW = shutil.which('nextflow')
needs_nextflow = pytest.mark.skipif(NEXTFLOW is None, reason='Nextflow is not on PATH; the tests above read the config instead')

# A run is slower to start, so these only cover each way the caps can go: the
# machine's, a scheduler's, the test's and highmem's.
RUN_CASES = [case for case in CASES if case[0] in ('docker', 'singularity', 'singularity,slurm', 'docker,test', 'slurm,highmem')]


def _nextflow(args, parser):
    """Nextflow's output for these arguments, run from a directory of its own
    so that nothing lands in the repository."""
    with tempfile.TemporaryDirectory() as cwd:
        result = subprocess.run(
            [NEXTFLOW, '-log', os.path.join(cwd, 'nextflow.log')] + args, cwd=cwd,
            env=dict(os.environ, NXF_SYNTAX_PARSER=parser), capture_output=True, text=True, timeout=600,
        )
    assert result.returncode == 0, f'nextflow {" ".join(args)} ({parser} parser):\n{result.stdout}{result.stderr}'
    return result.stdout


def _profile(profiles):
    return ['-profile', profiles] if profiles else []


def _config_caps(project, profiles, parser):
    """The caps `nextflow config -flat` gives."""
    flat = {}
    for line in _nextflow(['-q', 'config', '-flat'] + _profile(profiles) + [str(project)], parser).splitlines():
        name, sep, value = line.partition(' = ')
        if sep:
            flat[name] = int(value) if value.isdigit() else value.strip("'")
    return tuple(flat.get('params.' + cap) for cap in CAPS)


def _run_caps(profiles, parser):
    """The caps AMAP's --help gives as the defaults of a run."""
    text = _nextflow(['run', str(ROOT)] + _profile(profiles) + ['--help'], parser)
    caps = []
    for cap in CAPS:
        default = re.search(r'--' + cap + r'\s.*?\(Default: ([^)]+)\)', text)
        assert default, f'--help no longer gives the default of --{cap}'
        caps.append(int(default.group(1)) if default.group(1).isdigit() else default.group(1))
    return tuple(caps)


@pytest.fixture(scope='module')
def nextflow_caps():
    """The caps Nextflow gives every case, in both orders and under both
    parsers, by `nextflow config` and by a run, and what this machine's CPUs
    come to for Nextflow."""
    config_jobs = [(order, parser) for case in CASES for order in _orders(case[0]) for parser in PARSERS]
    config_jobs += [('', parser) for parser in PARSERS]
    run_jobs = [(order, parser) for case in RUN_CASES for order in _orders(case[0]) for parser in PARSERS]
    with tempfile.TemporaryDirectory() as probe:
        Path(probe, 'nextflow.config').write_text('params.max_cpus = Runtime.runtime.availableProcessors()\n')
        with ThreadPoolExecutor(max_workers=min(8, os.cpu_count() or 1)) as pool:
            machine = pool.submit(_config_caps, probe, '', 'v2')
            config = pool.map(lambda job: _config_caps(ROOT, *job), config_jobs)
            run = pool.map(lambda job: _run_caps(*job), run_jobs)
            return {
                'config': dict(zip(config_jobs, config)),
                'run': dict(zip(run_jobs, run)),
                'machine': machine.result()[0],
            }


def _want(caps, machine):
    return tuple(machine if value == MACHINE else value for value in caps)


@needs_nextflow
@pytest.mark.parametrize('profiles,cpus,memory,time', CASES, ids=[case[0] for case in CASES])
def test_caps_from_nextflow_config(nextflow_caps, profiles, cpus, memory, time):
    for order in _orders(profiles):
        for parser in PARSERS:
            got = nextflow_caps['config'][(order, parser)]
            assert got == _want((cpus, memory, time), nextflow_caps['machine']), f'-profile {order}, {parser} parser'


@needs_nextflow
@pytest.mark.parametrize('profiles,cpus,memory,time', RUN_CASES, ids=[case[0] for case in RUN_CASES])
def test_caps_in_a_run(nextflow_caps, profiles, cpus, memory, time):
    for order in _orders(profiles):
        for parser in PARSERS:
            got = nextflow_caps['run'][(order, parser)]
            assert got == _want((cpus, memory, time), nextflow_caps['machine']), f'-profile {order}, {parser} parser'


@needs_nextflow
def test_nextflow_agrees_with_the_reading(nextflow_caps):
    for (profiles, parser), got in nextflow_caps['config'].items():
        assert got == _want(read_caps(profiles, parser), nextflow_caps['machine']), (
            f'-profile {profiles or "(none)"}, {parser} parser'
        )
