#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PANGENE aligns each sample's proteins onto every other sample's genome.

It stages all the genomes in one folder, aligns against the ones named *.fasta
and leaves out the one named like the sample's own genome. PATHOTYPR's renamed
genomes fit that; with --skip_typing there were none, and PANGENE wrote empty
files for every sample. Now PANGENE_GENOME copies each genome of the sample
sheet to <sampleId>.fasta, as the sheet's names may clash (assembly.fasta from
two assemblies) or end in .fna.

These tests run both scripts from modules/local/pangene.nf the way Nextflow
renders them, inputs staged as links, with stand-ins for miniprot and pangene
that record what they were given. How main.nf wires the two together needs
Nextflow itself and is not covered here.
"""

import os
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MODULE = ROOT / 'modules' / 'local' / 'pangene.nf'


def _process(name):
    """A process's script block as Nextflow strips it, its stageAs patterns and
    the file names its outputs declare."""
    text = MODULE.read_text()
    block = re.search(rf'^process {name} \{{\n(.*?)^\}}', text, re.S | re.M).group(1)
    script = re.search(r'^    script:\n    """\n(.*?\n)    """\n', block, re.S | re.M).group(1)
    lines = script.split('\n')
    script = '\n'.join(line[4:] if line.startswith('    ') else line for line in lines)
    outputs = re.search(r'^    output:\n(.*?)\n\n', block, re.S | re.M).group(1)
    return (script, dict(re.findall(r"path\(?\s*(\w+), stageAs: '([^']+)'", block)),
            re.findall(r'path\("([^"]+)"\)', outputs))


def _render(script, values):
    """Groovy's triple-quoted string, for the little these scripts use: \\\\, \\$,
    \\t and ${...}. Anything else fails, rather than render differently from Nextflow."""
    out, i = [], 0
    while i < len(script):
        c = script[i]
        if c == '\\':
            nxt = script[i + 1]
            if nxt not in '\\$t':
                raise ValueError(f'escape \\{nxt} not handled here')
            out.append('\t' if nxt == 't' else nxt)
            i += 2
        elif c == '$':
            if script[i + 1] != '{':
                raise ValueError(f'bare $ at {script[i:i + 20]!r}')
            end = script.index('}', i)
            out.append(values[script[i + 2:end]])
            i = end + 1
        else:
            out.append(c)
            i += 1
    return ''.join(out)


def _stage(work, name, target):
    """A task input as Nextflow stages it: a link in the task directory."""
    assert target.is_file(), f'nothing to stage at {target}'
    link = work / name
    link.parent.mkdir(parents=True, exist_ok=True)
    link.symlink_to(target)
    return name


def _run(work, script, env=None):
    (work / '.command.sh').write_text(script)
    subprocess.run(['bash', '-euo', 'pipefail', '.command.sh'], cwd=work, env=env,
                   check=True, capture_output=True, text=True)


def _copy(tmp_path, sample, genome):
    """Run PANGENE_GENOME on one genome; the copy its output declares."""
    script, stage, (output,) = _process('PANGENE_GENOME')
    work = tmp_path / 'work' / sample
    work.mkdir(parents=True)
    staged = _stage(work, stage.get('fasta', '*').replace('*', genome.name), genome)
    _run(work, _render(script, {'sampleId': sample, 'fasta': staged}))
    return work / _render(output, {'sampleId': sample})


# Stand-ins: miniprot records the genome it is asked to align onto, and the
# rest write something where the script expects output.
STAND_INS = {
    'miniprot': '#!/bin/bash\nfor a in "$@"; do case "$a" in refs/*) echo "$a" >> aligned.txt ;; esac; done\n',
    'pangene': '#!/bin/bash\nprintf "H\\tVN:1.0\\n"\n',
    'k8': '#!/bin/bash\necho stand-in\n',
    'gfatools': '#!/bin/bash\necho stand-in\n',
}


def _pangene(tmp_path, sample, own, genomes):
    """Run PANGENE for one sample, its own genome and every sample's staged as
    main.nf hands them over; the genomes miniprot was asked to align onto."""
    script, stage, _outputs = _process('PANGENE')
    work = tmp_path / 'pangene' / sample
    work.mkdir(parents=True)
    for genome in genomes:
        _stage(work, stage['ref_fastas'].replace('*', genome.name), genome)
    for suffix in ('tsv', 'faa', 'fna'):
        (work / f'{sample}.AMAP.{suffix}').write_text('')
    tools = tmp_path / 'tools'
    tools.mkdir(exist_ok=True)
    for tool, text in STAND_INS.items():
        (tools / tool).write_text(text)
        (tools / tool).chmod(0o755)
    rendered = _render(script, {
        'sampleId': sample,
        'fixfaa': f'{sample}.AMAP.faa',
        'fixfna': f'{sample}.AMAP.fna',
        'fasta': _stage(work, own.name, own),
        # Enough jobs at once for every genome below: the macOS /bin/bash 3.2
        # calls an empty array unbound under set -u, the container's bash 5 does not.
        'task.cpus': '16',
    })
    _run(work, rendered, dict(os.environ, PATH=f'{tools}{os.pathsep}{os.environ["PATH"]}'))
    aligned = work / 'aligned.txt'
    return sorted(aligned.read_text().split()) if aligned.exists() else []


def _genome(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


class TestPangeneGenome:

    def test_copy_is_named_after_the_sample(self, tmp_path):
        genome = _genome(tmp_path / 'data' / 'K12' / 'assembly.fna', '>NC_000913.3\nACGT\n')
        copy = _copy(tmp_path, 'K12', genome)
        assert copy.name == 'K12.fasta'
        assert copy.read_text() == '>NC_000913.3\nACGT\n'
        assert not copy.is_symlink()

    def test_genome_already_named_like_its_copy(self, tmp_path):
        # Staged at the task's top level, the copy would be its own input
        genome = _genome(tmp_path / 'data' / 'S1.fasta', '>S1\nACGT\n')
        copy = _copy(tmp_path, 'S1', genome)
        assert copy.read_text() == '>S1\nACGT\n'
        assert genome.read_text() == '>S1\nACGT\n'


class TestPangene:

    def test_untyped_genomes_are_aligned_against(self, tmp_path):
        """Two assemblies called assembly.fna and one already called O157.fasta:
        each sample is aligned onto the other two, never onto itself."""
        data = tmp_path / 'data'
        sheet = {
            'K12': _genome(data / 'K12' / 'assembly.fna', '>a\nACGT\n'),
            'Sakai': _genome(data / 'Sakai' / 'assembly.fna', '>b\nACGT\n'),
            'O157': _genome(data / 'O157.fasta', '>c\nACGT\n'),
        }
        copies = {sample: _copy(tmp_path, sample, genome) for sample, genome in sheet.items()}
        genomes = list(copies.values())
        assert _pangene(tmp_path, 'K12', copies['K12'], genomes) == ['refs/O157.fasta', 'refs/Sakai.fasta']
        assert _pangene(tmp_path, 'O157', copies['O157'], genomes) == ['refs/K12.fasta', 'refs/Sakai.fasta']

    def test_typed_genomes_as_before(self, tmp_path):
        data = tmp_path / 'typed'
        genomes = [_genome(data / f'{name}.fasta', f'>{name}\nACGT\n') for name in ('L4_CDC1551', 'A4_Mbovis')]
        assert _pangene(tmp_path, 'L4_CDC1551', genomes[0], genomes) == ['refs/A4_Mbovis.fasta']

    def test_one_sample_has_nothing_to_align_against(self, tmp_path):
        copy = _copy(tmp_path, 'K12', _genome(tmp_path / 'data' / 'K12.fna', '>a\nACGT\n'))
        assert _pangene(tmp_path, 'K12', copy, [copy]) == []
        assert (tmp_path / 'pangene' / 'K12' / 'K12.gfa').read_text() == 'H\tVN:1.0\n'
