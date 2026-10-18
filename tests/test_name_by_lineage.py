#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PATHOTYPR names every genome <lineage>_<sampleId>.

pathotypr writes one summary row per record that carries lineage markers,
named after the record, so a draft gets one row per contig with markers. Read
off the first row, the name was a contig's: CDC1551 cut into 40 contigs, sample
CDC1551d, came out as L4_contig_14, and so did every file after it. The
summaries below are pathotypr's own, from that draft and from the closed genome.

The module tests run PATHOTYPR's script the way Nextflow renders it, with a
stand-in for pathotypr that writes a given summary, so they also fail if the
naming in modules/local/pathotypr.nf goes back to reading the summary's rows.
"""

import hashlib
import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'bin'))

from name_by_lineage import read_rows, sample_lineage, write_fasta  # noqa: E402

SCRIPT = ROOT / 'bin' / 'name_by_lineage.py'
MODULE = ROOT / 'modules' / 'local' / 'pathotypr.nf'

HEADER = 'genome\tlineage:count\tmajor_lineage\n'
# The 40-contig CDC1551 draft: only the forward contigs with markers have a row
DRAFT = HEADER + (
    'contig_14\tL4:5 L4;L4.1:3 L4;L4.1;L4.1.1:1\tL4\n'
    'contig_15\tL4:1\tL4\n'
    'contig_16\tL4:5 L4;L4.1:2 L4;L4.1;L4.1.1:1\tL4\n'
    'contig_18\tL4:1 L4;L4.1:1 L4;L4.1;L4.1.1:1\tL4\n'
    'contig_35\tL4:4 L4;L4.1:1 L4;L4.1;L4.1.1:1 L4;L4.1;L4.1.1;L4.1.1.3:1\tL4\n'
    'contig_36\tL4:4 L4;L4.1:1 L4;L4.1;L4.1.1:1 L4;L4.1;L4.1.1;L4.1.1.3:1\tL4\n'
    'contig_4\tL4:8 L4;L4.1:3 L4;L4.1;L4.1.1:3 L4;L4.1;L4.1.1;L4.1.1.3:1\tL4\n'
    'contig_6\tL4:8 L4;L4.1:4 L4;L4.1;L4.1.1:3 L4;L4.1;L4.1.1;L4.1.1.3:2\tL4\n'
    'contig_7\tL4:6 L4;L4.1:1 L4;L4.1;L4.1.1:1 L4;L4.1;L4.1.1;L4.1.1.3:1\tL4\n')
CLOSED_COUNTS = 'L4:79 L4;L4.1:24 L4;L4.1;L4.1.1:17 L4;L4.1;L4.1.1;L4.1.1.3:7'
# pathotypr names the record by its whole header line, description and all
CLOSED_NCBI = HEADER + (
    f'NC_002755.2 Mycobacterium tuberculosis CDC1551, complete sequence\t{CLOSED_COUNTS}\tL4\n')
CLOSED = HEADER + f'CDC1551\t{CLOSED_COUNTS}\tL4\n'
# What PATHOTYPR writes when pathotypr fails
PLACEHOLDER = 'Sample\tLineage\tName\n'


def _summary(tmp_path, text, name='summary.tsv'):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def _name(tmp_path, sample, summary_text, fasta_bytes=b'>c1\nACGT\n>c2\nGGCC\n'):
    """Run the script as PATHOTYPR does; (printed name, files it wrote)."""
    summary = _summary(tmp_path, summary_text)
    genome = tmp_path / 'genome.fa'
    genome.write_bytes(fasta_bytes)
    out = tmp_path / 'out'
    out.mkdir()
    done = subprocess.run([sys.executable, str(SCRIPT), sample, summary, str(genome)],
                          cwd=out, check=True, capture_output=True, text=True)
    return done.stdout.strip(), {p.name: p.read_bytes() for p in out.iterdir()}


class TestSampleLineage:

    def test_draft_is_named_after_the_sample_not_a_contig(self, tmp_path):
        name, files = _name(tmp_path, 'CDC1551d', DRAFT)
        assert name == 'L4_CDC1551d'
        assert list(files) == ['L4_CDC1551d.fasta']

    def test_order_of_the_rows_does_not_matter(self, tmp_path):
        rows = DRAFT.splitlines(keepends=True)
        assert _name(tmp_path, 'CDC1551d', rows[0] + ''.join(reversed(rows[1:])))[0] == 'L4_CDC1551d'

    def test_closed_genome_named_as_before(self, tmp_path):
        assert _name(tmp_path, 'CDC1551', CLOSED)[0] == 'L4_CDC1551'

    def test_a_header_with_a_description_does_not_end_up_in_the_name(self, tmp_path):
        # Read off the summary, this was 'L4_NC_002755.2 Mycobacterium ...',
        # spaces and all, which no file or command downstream survives
        assert _name(tmp_path, 'CDC1551', CLOSED_NCBI)[0] == 'L4_CDC1551'

    @pytest.mark.parametrize('text', [HEADER, PLACEHOLDER, ''])
    def test_no_lineage_gives_NA(self, tmp_path, text):
        name, files = _name(tmp_path, 'G1', text)
        assert name == 'NA_G1'
        assert list(files) == ['NA_G1.fasta']

    def test_missing_summary_gives_no_rows(self, tmp_path):
        assert read_rows(str(tmp_path / 'absent.tsv')) == []

    def test_contigs_that_disagree_go_with_the_most_markers(self, tmp_path):
        rows = read_rows(_summary(tmp_path, HEADER + 'a\tL2:3\tL2\nb\tL4:2 L4;L4.1:2\tL4\nc\tL4:2\tL4\n'))
        assert sample_lineage(rows) == 'L4'

    def test_a_tie_goes_to_the_first_alphabetically_as_in_pathotypr(self, tmp_path):
        rows = read_rows(_summary(tmp_path, HEADER + 'a\tL4:3\tL4\nb\tL2:3\tL2\n'))
        assert sample_lineage(rows) == 'L2'

    def test_a_sublineage_never_outvotes_its_lineage(self, tmp_path):
        rows = read_rows(_summary(tmp_path, HEADER + 'a\tL4:2 L4;L4.1:2\tL4\nb\tL4:1 L4;L4.1:1\tL4\n'))
        assert sample_lineage(rows) == 'L4'

    def test_one_row_keeps_the_call_pathotypr_made(self, tmp_path):
        rows = read_rows(_summary(tmp_path, HEADER + 'a\tL4:5 L2:1\tL2\n'))
        assert sample_lineage(rows) == 'L2'

    def test_unclassified_is_no_lineage(self, tmp_path):
        assert sample_lineage(read_rows(_summary(tmp_path, HEADER + 'a\t\tUnclassified\n'))) is None


class TestWriteFasta:

    def test_one_record_is_renamed_and_its_sequence_left_alone(self, tmp_path):
        src = tmp_path / 'in.fa'
        src.write_bytes(b'>NC_002755.2 Mycobacterium tuberculosis CDC1551\nACGTN\nAC\n')
        write_fasta(str(src), str(tmp_path / 'out.fa'), 'L4_CDC1551')
        assert (tmp_path / 'out.fa').read_bytes() == b'>L4_CDC1551\nACGTN\nAC\n'

    def test_the_contigs_of_a_draft_keep_their_names(self, tmp_path):
        draft = b'>L4_CDC1551d\nACGT\n>contig_2 len=4\nGGCC\n>contig_14\nTTAA\n'
        src = tmp_path / 'in.fa'
        src.write_bytes(draft)
        write_fasta(str(src), str(tmp_path / 'out.fa'), 'L4_CDC1551d')
        assert (tmp_path / 'out.fa').read_bytes() == draft

    def test_lines_end_the_way_the_old_awk_ended_them(self, tmp_path):
        src = tmp_path / 'in.fa'
        src.write_bytes(b'>a\r\nAC\r\nGT')
        write_fasta(str(src), str(tmp_path / 'out.fa'), 'L4_x')
        assert (tmp_path / 'out.fa').read_bytes() == b'>L4_x\nAC\r\nGT\n'


# ---------------------------------------------------------------------------
# PATHOTYPR's own script
# ---------------------------------------------------------------------------

def _module():
    """PATHOTYPR's script block as Nextflow strips it, and where its genome is staged."""
    text = MODULE.read_text()
    # the Groovy lines before the string (the fallbacks for staged files) are not part of it
    script = re.search(r'^    script:\n(?:    (?!""")[^\n]*\n)*?    """\n(.*?)\n    """\n', text, re.S | re.M)
    stage = re.search(r"path\(fasta, stageAs: '([^']+)'\)", text)
    return textwrap.dedent(script.group(1)) + '\n', stage.group(1) if stage else '*'


def _render(script, values):
    """Groovy's triple-quoted string, for the little this script uses: \\\\, \\$
    and ${...}. Anything else fails, rather than render differently from Nextflow."""
    out, i = [], 0
    while i < len(script):
        c = script[i]
        if c == '\\':
            nxt = script[i + 1]
            if nxt not in '\\$':
                raise ValueError(f'escape \\{nxt} not handled here')
            out.append(nxt)
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


STAND_IN = """#!/bin/bash
# pathotypr: copies $FAKE_SUMMARY to <prefix>_summary.tsv, or fails without one
while [ $# -gt 0 ]; do
    case "$1" in
        -i) input="$2"; shift ;;
        --output-prefix) prefix="$2"; shift ;;
    esac
    shift
done
[ -s "$input" ] || { echo "no genome at $input" >&2; exit 1; }
[ -n "$FAKE_SUMMARY" ] || exit 1
: > "$prefix.tsv"
cp "$FAKE_SUMMARY" "${prefix}_summary.tsv"
"""


def _run_module(tmp_path, sample, genome_name, genome, summary_text=None):
    """Run PATHOTYPR's script in a task directory of its own, its genome staged
    as a link the way Nextflow stages it. (SAMPLE_NAME, the task directory)."""
    script, stage = _module()
    original = tmp_path / 'data' / genome_name
    original.parent.mkdir()
    original.write_bytes(genome)
    work = tmp_path / 'work'
    staged = work / stage.replace('*', genome_name)
    staged.parent.mkdir(parents=True, exist_ok=True)
    staged.symlink_to(original)

    tools = tmp_path / 'tools'
    tools.mkdir()
    (tools / 'pathotypr').write_text(STAND_IN)
    (tools / 'python3').write_text(f'#!/bin/bash\nexec "{sys.executable}" "$@"\n')
    for tool in tools.iterdir():
        tool.chmod(0o755)

    rendered = _render(script, {
        'sampleId': sample,
        'fasta': staged.relative_to(work).as_posix(),
        'markers': 'snps_lineages.tsv',
        'markers_ref': 'MTB_ancestor_reference.fasta',
        'projectDir': str(ROOT),
    })
    # Nextflow reads env() outputs after the script, as this does
    (work / '.command.sh').write_text(rendered + 'printf %s "$SAMPLE_NAME" > .sample_name\n')
    env = dict(os.environ, PATH=f'{tools}{os.pathsep}{os.environ["PATH"]}')
    env.pop('FAKE_SUMMARY', None)
    if summary_text is not None:
        env['FAKE_SUMMARY'] = _summary(tmp_path, summary_text, 'pathotypr_summary.tsv')
    subprocess.run(['bash', '.command.sh'], cwd=work, env=env, check=True, capture_output=True)
    return (work / '.sample_name').read_text(), work


def _outputs(work, pattern):
    """What PATHOTYPR's output glob picks up: Nextflow leaves out its inputs."""
    return sorted(p.name for p in work.glob(pattern) if not p.is_symlink())


DRAFT_FASTA = b'>L4_CDC1551d\nACGT\n>contig_2\nGGCC\n>contig_14\nTTAA\n'


class TestModule:

    def test_draft_is_named_after_the_sample(self, tmp_path):
        name, work = _run_module(tmp_path, 'CDC1551d', 'CDC1551d.fasta', DRAFT_FASTA, DRAFT)
        assert name == 'L4_CDC1551d'
        assert _outputs(work, '*.fasta') == ['L4_CDC1551d.fasta']
        assert (work / 'L4_CDC1551d.fasta').read_bytes() == DRAFT_FASTA
        assert _outputs(work, '*_summary.tsv') == ['CDC1551d_summary.tsv']

    def test_closed_genome_named_as_before(self, tmp_path):
        name, work = _run_module(tmp_path, 'CDC1551', 'CDC1551.fasta', b'>CDC1551\nACGT\n', CLOSED)
        assert name == 'L4_CDC1551'
        assert (work / 'L4_CDC1551.fasta').read_bytes() == b'>L4_CDC1551\nACGT\n'

    def test_closed_genome_with_a_description(self, tmp_path):
        genome = b'>NC_002755.2 Mycobacterium tuberculosis CDC1551, complete sequence\nACGT\n'
        name, work = _run_module(tmp_path, 'CDC1551', 'cdc1551.fna', genome, CLOSED_NCBI)
        assert name == 'L4_CDC1551'
        assert _outputs(work, '*.fasta') == ['L4_CDC1551.fasta']

    def test_failed_typing_gives_NA(self, tmp_path):
        name, work = _run_module(tmp_path, 'K12', 'K12.fasta', b'>K12\nACGT\n')
        assert name == 'NA_K12'
        assert _outputs(work, '*.fasta') == ['NA_K12.fasta']
        assert (work / 'K12_summary.tsv').read_text() == PLACEHOLDER

    def test_an_input_named_like_the_output_is_left_alone(self, tmp_path):
        # Sample G26689 in a file already called L4_G26689.fasta: written next
        # to its input, the renamed genome would go through the link into it
        genome = b'>G26689\nACGT\n'
        name, work = _run_module(tmp_path, 'G26689', 'L4_G26689.fasta', genome,
                                 HEADER + f'G26689\t{CLOSED_COUNTS}\tL4\n')
        assert name == 'L4_G26689'
        assert hashlib.md5((tmp_path / 'data' / 'L4_G26689.fasta').read_bytes()).hexdigest() \
            == hashlib.md5(genome).hexdigest()
        assert _outputs(work, '*.fasta') == ['L4_G26689.fasta']
        assert (work / 'L4_G26689.fasta').read_bytes() == b'>L4_G26689\nACGT\n'
