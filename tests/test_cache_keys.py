#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Every AMAP file a Nextflow process runs must be in that process's cache key.

Nextflow hashes a task's inputs and its script, not the files the script calls
by path: with `python3 ${projectDir}/bin/check_cds.py` alone, an edit to
check_cds.py leaves the cache key as it was, and `-resume` hands back what the
old code produced. So each process that runs AMAP code lists those files in
`ext code` and names `task.ext.code` in its script, which makes Nextflow hash
them (Nextflow 24.04 and later). A script that imports the `classify` package
needs the package listed too.

The other way round, Nextflow adds to a task's cache key every bin/ entry whose
name appears in its script as a bare word, so a stray word ties a task to code
it never runs: `pathotypr classify` tied PATHOTYPR, and Bakta after it, to
bin/classify. A bare word that names a bin/ entry must therefore be code the
process declares.

These tests read the .nf files the way Nextflow does and need no Nextflow.
"""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "bin"
NF_FILES = sorted((ROOT / "modules" / "local").glob("*.nf")) + [ROOT / "main.nf"]

# What Nextflow splits a script on when it looks for bin/ entries
# (TaskHasher.getTaskBinEntries).
NF_WORD_SPLIT = re.compile(r"[ \t\n\r\f()\[\]{};&|<>`]+")
BIN_CALL = re.compile(r"\$\{projectDir\}/bin/([\w./-]+)")
EXT_CODE = re.compile(r'^    ext code: \{ files\("\$\{projectDir\}/bin/([^"]+)"\) \}$', re.M)
# `task.ext.code` itself, as a property: `task.ext.code*.name` would make
# Nextflow look up a variable of that name instead, and hash null.
NAMES_EXT_CODE = re.compile(r"\btask\.ext\.code(?![\w*])")
IMPORTS_CLASSIFY = re.compile(r"^\s*(?:from|import)\s+classify\b", re.M)


def _processes():
    """(file, name, directives, script) of every process, the script being
    what Nextflow hashes: the whole `script:` block."""
    found = []
    for nf in NF_FILES:
        for m in re.finditer(r"^process (\w+) \{\n(.*?)^\}", nf.read_text(), re.S | re.M):
            parts = re.split(r"^[ \t]*script:[ \t]*\n", m.group(2), maxsplit=1, flags=re.M)
            assert len(parts) == 2, f"{nf.name}: process {m.group(1)} has no script block"
            directives, script = parts
            found.append((nf.relative_to(ROOT).as_posix(), m.group(1), directives, script))
    return found


def _alternatives(pattern):
    """Expand one {a,b,...} group, the only brace syntax used here."""
    m = re.search(r"\{([^{}]*)\}", pattern)
    if not m:
        return [pattern]
    return [pattern[:m.start()] + alt + pattern[m.end():] for alt in m.group(1).split(",")]


def _files(pattern):
    """The files a files() pattern relative to bin/ matches, directories left out."""
    return {p for p in BIN.glob(pattern) if p.is_file()}


def _declared(directives):
    patterns = EXT_CODE.findall(directives)
    assert len(patterns) <= 1, "more than one `ext code` line"
    return _alternatives(patterns[0]) if patterns else []


PROCESSES = _processes()
IDS = [f"{name}" for _, name, _, _ in PROCESSES]


def test_finds_the_processes():
    names = {name for _, name, _, _ in PROCESSES}
    assert {"PATHOTYPR", "TYPE_CASES", "FIX_CDS", "GET_OUTPUTS", "FINAL_REPORT"} <= names


@pytest.mark.parametrize("nf,name,directives,script", PROCESSES, ids=IDS)
def test_every_file_it_runs_is_declared(nf, name, directives, script):
    declared = set().union(*(_files(p) for p in _declared(directives)))
    called = {BIN / rel for rel in BIN_CALL.findall(script)}
    missing = sorted(p.relative_to(BIN).as_posix() for p in called - declared)
    assert not missing, (
        f"{name} ({nf}) runs {missing} without listing them in `ext code`, so "
        f"editing them would not re-run it on -resume")


@pytest.mark.parametrize("nf,name,directives,script", PROCESSES, ids=IDS)
def test_the_package_is_declared_where_imported(nf, name, directives, script):
    declared = set().union(*(_files(p) for p in _declared(directives)))
    importers = sorted(rel for rel in BIN_CALL.findall(script)
                       if rel.endswith(".py") and IMPORTS_CLASSIFY.search((BIN / rel).read_text()))
    if importers:
        missing = sorted(p.relative_to(BIN).as_posix() for p in _files("classify/*.py") - declared)
        assert not missing, (
            f"{name} ({nf}) runs {importers}, which import classify, but does not "
            f"list {missing} in `ext code`")


@pytest.mark.parametrize("nf,name,directives,script", PROCESSES, ids=IDS)
def test_declared_code_exists(nf, name, directives, script):
    empty = [p for p in _declared(directives) if not _files(p)]
    assert not empty, f"{name} ({nf}) lists {empty} in `ext code`, which match no file"


@pytest.mark.parametrize("nf,name,directives,script", PROCESSES, ids=IDS)
def test_declared_code_is_hashed(nf, name, directives, script):
    if _declared(directives):
        assert NAMES_EXT_CODE.search(script), (
            f"{name} ({nf}) lists its code in `ext code` but its script never names "
            f"`task.ext.code`, so Nextflow does not hash it")


@pytest.mark.parametrize("nf,name,directives,script", PROCESSES, ids=IDS)
def test_no_bare_word_ties_it_to_other_code(nf, name, directives, script):
    entries = {p.name: p for p in BIN.iterdir() if p.name != "__pycache__"}
    declared = set().union(*(_files(p) for p in _declared(directives)))
    stray = sorted(
        word for word in set(NF_WORD_SPLIT.split(script)) & entries.keys()
        if not any(f == entries[word] or f.is_relative_to(entries[word]) for f in declared))
    assert not stray, (
        f"{name} ({nf}) names {', '.join('bin/' + w for w in stray)} as a bare word, "
        f"which Nextflow adds to its cache key (any directory, and any script with the "
        f"executable bit) though the process does not declare that code; quote the "
        f"word or reword it")


def test_the_split_matches_nextflow():
    # pathotypr's subcommand is quoted for exactly this reason
    assert "classify" in NF_WORD_SPLIT.split("    pathotypr classify \\\\\n")
    assert "classify" not in NF_WORD_SPLIT.split("    pathotypr 'classify' \\\\\n")
    assert "classify" in NF_WORD_SPLIT.split("(classify)")
