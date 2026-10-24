#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PANGENE builds one pangenome for the whole run.

pangene builds its graph from one set of proteins aligned to every genome.
PANGENE's set is the reference proteome and the genes of each sample that no
reference gene matches (bin/pangenome.py); it aligns the set to every sample's
genome, the sample's own included, runs pangene once and names each column of
the matrix after its genome. It used to align each sample's own proteins to
the other samples' genomes only, one pangene run per sample, and REPORT_INPUTS
merged those matrices on gene names that carry the sample's name, so no gene
was ever in two samples: on the first 2 Mb of CDC1551, M. bovis and H37Rv,
5379 rows in 6 columns, none of them in all six.

The genomes are told apart by file name, <name>.fasta: PATHOTYPR writes them
so, and with --skip_typing PANGENE_GENOME copies each genome of the sample
sheet to <sampleId>.fasta, as the sheet's names may clash (assembly.fasta from
two assemblies) or end in .fna.

These tests run both scripts from modules/local/pangene.nf the way Nextflow
renders them, every input staged as a link where the process stages it, on
tests/data/pangenome/: real genomes cut around RD4, which M. bovis lacks, and
around the RvD2 deletion of H37Rv, where CDC1551 and M. bovis carry three genes
no reference gene matches. miniprot, pangene, k8 and gfatools are stand-ins
that hand names on the way the real ones do, which is what the script relies
on; where the real ones are installed, as in the AMAP container, the last test
runs them.
"""

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'bin'))

import pangenome  # noqa: E402
from tests.test_sample_ids import process_script, render  # noqa: E402

MODULE = ROOT / 'modules' / 'local' / 'pangene.nf'
DATA = ROOT / 'tests' / 'data' / 'pangenome'
SAMPLES = ['A4_Mbovis', 'L4_CDC1551', 'L4_H37Rv']
GENOMES = [DATA / 'genomes' / f'{s}.fasta' for s in SAMPLES]
CORRECTED = sorted((DATA / 'corrected').iterdir())
REFERENCE = DATA / 'reference.faa'
RD4 = ['Rv1506c', 'Rv1507A', 'Rv1507c', 'Rv1508A', 'Rv1508c', 'Rv1509', 'Rv1510', 'Rv1511',
       'Rv1512', 'Rv1513', 'Rv1514c', 'Rv1515c', 'Rv1516c']


# ---------------------------------------------------------------------------
# The processes, as Nextflow stages and runs them
# ---------------------------------------------------------------------------

def _block(process):
    text = MODULE.read_text()
    return re.search(rf'^process {process} \{{\n(.*?)^\}}', text, re.S | re.M).group(1)


def _stage_patterns(process):
    """Each path input's stageAs pattern."""
    return dict(re.findall(r"path\(?\s*(\w+), stageAs: '([^']+)'", _block(process)))


def _expand(pattern):
    """A file pattern with one {a,b} group, the only one used here, as names."""
    m = re.search(r'\{([^{}]*)\}', pattern)
    if not m or m.group(0).startswith('{sampleId'):
        return [pattern]
    return [pattern[:m.start()] + alt + pattern[m.end():] for alt in m.group(1).split(',')]


def _outputs(process):
    """The file names a process's output block declares."""
    outputs = re.search(r'^    output:\n(.*?)\n\n', _block(process), re.S | re.M).group(1)
    return [name for pattern in re.findall(r'path\(?\s*"([^"]+)"', outputs) for name in _expand(pattern)]


def _link(work, staged, target):
    """A task input as Nextflow stages it: a link, where its stageAs puts it."""
    link = work / staged
    link.parent.mkdir(parents=True, exist_ok=True)
    link.symlink_to(Path(target).resolve())
    return staged


def _run(work, script, tools):
    """The script, with the shell nextflow.config gives every task."""
    (work / '.command.sh').write_text(script)
    env = dict(os.environ, PATH=f'{tools}{os.pathsep}{os.environ["PATH"]}',
               PYTHONDONTWRITEBYTECODE='1')
    done = subprocess.run(['bash', '-euo', 'pipefail', '.command.sh'], cwd=work, env=env,
                          capture_output=True, text=True)
    assert done.returncode == 0, f'{done.stdout[-2000:]}\n{done.stderr[-3000:]}'


def _copy(tmp_path, sample, genome, tools):
    """Run PANGENE_GENOME on one genome; the copy its output declares."""
    work = tmp_path / 'PANGENE_GENOME' / sample
    work.mkdir(parents=True)
    staged = _link(work, _stage_patterns('PANGENE_GENOME')['fasta'].replace('*', genome.name), genome)
    _run(work, render(process_script('pangene.nf', 'PANGENE_GENOME'),
                      {'sampleId': sample, 'fasta': staged}), tools)
    (output,) = _outputs('PANGENE_GENOME')
    return work / output.replace('${sampleId}', sample)


def _pangene(tmp_path, genomes, tools, cpus=4):
    """Run PANGENE on the genomes, the fixture's corrected tables and proteins
    and its reference proteome, staged as main.nf hands them over; the task
    directory."""
    stage = _stage_patterns('PANGENE')
    work = tmp_path / 'PANGENE'
    work.mkdir(parents=True)
    for genome in genomes:
        _link(work, stage['genomes'].replace('*', genome.name), genome)
    for f in CORRECTED:
        _link(work, stage['corrected'].replace('*', f.name), f)
    reference = _link(work, stage['ref_proteome'].replace('*', REFERENCE.name), REFERENCE)
    _run(work, render(process_script('pangene.nf', 'PANGENE'),
                      {'projectDir': ROOT, 'task.cpus': cpus, 'ref_proteome': reference}), tools)
    return work


def _matrix(work):
    """pangenome.Rtab: its genome columns, and each gene's row."""
    header, *rows = (work / 'pangenome.Rtab').read_text().splitlines()
    return header.split('\t')[1:], {r.split('\t')[0]: r.split('\t')[1:] for r in rows}


def _log(work, tool):
    """The arguments each call of a stand-in was given."""
    return [line.split('\t')[1:] for line in (work / 'tools.log').read_text().splitlines()
            if line.startswith(tool + '\t')]


def _genome(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


# What the real tools do with names, measured in the AMAP container
# (miniprot 0.18, pangene 1.1): pangene names a genome after its PAF without
# the folder, a gene after its protein's name up to the first colon, and
# leaves out a genome with no alignment; gfa2matrix names each genome
# <name>#0 and writes "Gene" and a tab when there is none.
MINIPROT = '''#!PYTHON
"""miniprot: a PAF line for every protein, on a genome with any sequence."""
import sys
genome, proteins = sys.argv[-2:]
with open('tools.log', 'a') as log:
    log.write('miniprot\\t%s\\t%s\\n' % (genome, proteins))
seq = ''.join(l.strip() for l in open(genome) if not l.startswith('>'))
if seq:
    for line in open(proteins):
        if line.startswith('>'):
            print('%s\\t100\\t0\\t100\\t+\\tctg\\t%d\\t0\\t300\\t300\\t300\\t0' % (line[1:].split()[0], len(seq)))
'''
PANGENE = '''#!PYTHON
"""pangene: a gene for each name up to the first colon, and a walk for each
PAF with a hit, named after the file; a PAF without one has no walk."""
import os, sys
with open('tools.log', 'a') as log:
    log.write('\\t'.join(['pangene'] + sys.argv[1:]) + '\\n')
genes, walks = [], []
for paf in [a for a in sys.argv[1:] if a.endswith('.paf')]:
    hits = [l.split('\\t')[0].split(':')[0] for l in open(paf) if l.strip()]
    genes += [g for g in dict.fromkeys(hits) if g not in genes]
    if hits:
        walks.append((os.path.basename(paf)[:-len('.paf')], hits))
for g in genes:
    print('S\\t%s\\t*\\tLN:i:1' % g)
for name, hits in walks:
    print('W\\t%s\\t0\\t%s\\t0\\t1\\t%s' % (name, name, ''.join('>' + g for g in hits)))
'''
K8 = '''#!PYTHON
"""k8 pangene.js: gfa2matrix, a gene by genome table with each genome named
<name>#<haplotype>; call, a header."""
import sys
_script, command, gfa = sys.argv[1:4]
genes, walks = [], []
for line in open(gfa):
    f = line.rstrip('\\n').split('\\t')
    if f[0] == 'S':
        genes.append(f[1])
    elif f[0] == 'W':
        walks.append(('%s#%s' % (f[1], f[2]), set(f[6].split('>')[1:])))
if command == 'gfa2matrix':
    print('Gene\\t' + '\\t'.join(name for name, _ in walks))
    for g in genes:
        print(g + '\\t' + '\\t'.join('1' if g in hits else '0' for _, hits in walks))
else:
    print('#stand-in')
'''
GFATOOLS = '''#!PYTHON
import sys
print('Number of segments: %d' % sum(1 for l in open(sys.argv[2]) if l.startswith('S')))
'''


def _tools(d, stand_ins):
    d.mkdir(parents=True, exist_ok=True)
    (d / 'python3').write_text(f'#!/bin/bash\nexec "{sys.executable}" "$@"\n')
    for name, text in stand_ins.items():
        (d / name).write_text(text.replace('PYTHON', sys.executable))
    for tool in d.iterdir():
        tool.chmod(0o755)
    return d


@pytest.fixture(scope='module')
def stand_ins(tmp_path_factory):
    return _tools(tmp_path_factory.mktemp('stand-ins'),
                  {'miniprot': MINIPROT, 'pangene': PANGENE, 'k8': K8, 'gfatools': GFATOOLS})


# ---------------------------------------------------------------------------
# PANGENE_GENOME
# ---------------------------------------------------------------------------

class TestPangeneGenome:

    def test_copy_is_named_after_the_sample(self, tmp_path, stand_ins):
        genome = _genome(tmp_path / 'data' / 'K12' / 'assembly.fna', '>NC_000913.3\nACGT\n')
        copy = _copy(tmp_path, 'K12', genome, stand_ins)
        assert copy.name == 'K12.fasta'
        assert copy.read_text() == '>NC_000913.3\nACGT\n'
        assert not copy.is_symlink()

    def test_genome_already_named_like_its_copy(self, tmp_path, stand_ins):
        # Staged at the task's top level, the copy would be its own input
        genome = _genome(tmp_path / 'data' / 'S1.fasta', '>S1\nACGT\n')
        copy = _copy(tmp_path, 'S1', genome, stand_ins)
        assert copy.read_text() == '>S1\nACGT\n'
        assert genome.read_text() == '>S1\nACGT\n'


# ---------------------------------------------------------------------------
# PANGENE
# ---------------------------------------------------------------------------

class TestPangene:

    def test_every_genome_is_searched_with_the_same_proteins(self, tmp_path, stand_ins):
        work = _pangene(tmp_path, GENOMES, stand_ins)
        aligned = _log(work, 'miniprot')
        # each genome once, the sample's own among them
        assert sorted(genome for genome, _ in aligned) == [f'genomes/{s}.fasta' for s in SAMPLES]
        assert {proteins for _, proteins in aligned} == {'proteins.faa'}
        # and one pangene run over all of them
        (run,) = _log(work, 'pangene')
        assert sorted(a for a in run if a.endswith('.paf')) == [f'paf/{s}.paf' for s in SAMPLES]

    def test_the_proteins_are_the_reference_and_the_genes_outside_it(self, tmp_path, stand_ins):
        work = _pangene(tmp_path, GENOMES, stand_ins)
        written = dict(pangenome.read_fasta(work / 'proteins.faa'))
        reference = [name for name, _ in pangenome.read_fasta(REFERENCE)]
        names = list(written)
        # every reference protein, named <gene>:<protein> as pangene reads it
        assert names[:len(reference)] == [f"{r.rsplit('_', 1)[0]}:{r}" for r in reference]
        # then the genes no reference gene matches, in order of sample
        outside = [n.split(':')[0] for n in names[len(reference):]]
        assert outside == ['A4_Mbovis_01832', 'A4_Mbovis_01833', 'A4_Mbovis_01834',
                           'L4_CDC1551_01839', 'L4_CDC1551_01841']
        # L4_CDC1551_01840 is A4_Mbovis_01833 to the residue, so it is not
        cdc = dict(pangenome.read_fasta(DATA / 'corrected' / 'L4_CDC1551.AMAP.faa'))
        bovis = dict(pangenome.read_fasta(DATA / 'corrected' / 'A4_Mbovis.AMAP.faa'))
        same = [n for n in cdc if n.startswith('L4_CDC1551_01840;')]
        assert len(same) == 1
        assert cdc[same[0]] == next(s for n, s in bovis.items() if n.startswith('A4_Mbovis_01833;'))
        # and none ends in the stop's asterisk
        assert not [n for n, s in written.items() if s.endswith('*')]

    def test_a_column_per_genome_named_after_it(self, tmp_path, stand_ins):
        # A genome none of the proteins is found in has no walk in the graph,
        # and still gets its column, of zeros
        empty = _genome(tmp_path / 'data' / 'NA_Empty.fasta', '>NA_Empty\n')
        work = _pangene(tmp_path, GENOMES + [empty], stand_ins)
        columns, rows = _matrix(work)
        assert columns == SAMPLES + ['NA_Empty']
        assert len(rows) == 27      # 22 reference genes and 5 outside the reference
        assert all(values == ['1', '1', '1', '0'] for values in rows.values())

    def test_untyped_genomes_under_their_sample_ids(self, tmp_path, stand_ins):
        # --skip_typing: two assemblies called assembly.fna and one already
        # called O157.fasta, each copied by PANGENE_GENOME first
        data = tmp_path / 'data'
        sheet = {'K12': _genome(data / 'K12' / 'assembly.fna', '>a\nACGT\n'),
                 'Sakai': _genome(data / 'Sakai' / 'assembly.fna', '>b\nACGT\n'),
                 'O157': _genome(data / 'O157.fasta', '>c\nACGT\n')}
        copies = [_copy(tmp_path, sample, genome, stand_ins) for sample, genome in sheet.items()]
        work = _pangene(tmp_path, copies, stand_ins)
        assert sorted(g for g, _ in _log(work, 'miniprot')) == [
            'genomes/K12.fasta', 'genomes/O157.fasta', 'genomes/Sakai.fasta']
        columns, _ = _matrix(work)
        assert columns == ['K12', 'O157', 'Sakai']

    def test_one_genome(self, tmp_path, stand_ins):
        work = _pangene(tmp_path, GENOMES[:1], stand_ins)
        columns, rows = _matrix(work)
        assert columns == ['A4_Mbovis']
        assert rows and all(values == ['1'] for values in rows.values())

    def test_as_many_genomes_as_jobs_at_once(self, tmp_path, stand_ins):
        # 4 CPUs run two miniprot jobs at once: after two genomes no job is
        # left to wait for, an empty array the bash macOS has as /bin/bash
        # calls unbound under set -u
        work = _pangene(tmp_path, GENOMES[:2], stand_ins, cpus=4)
        columns, _ = _matrix(work)
        assert columns == SAMPLES[:2]

    def test_what_it_publishes_is_what_it_declares(self, tmp_path, stand_ins):
        work = _pangene(tmp_path, GENOMES, stand_ins)
        declared = _outputs('PANGENE')
        assert sorted(declared) == ['pangenome.Rtab', 'pangenome.bubbles.txt',
                                    'pangenome.gfa', 'pangenome.stat']
        published = re.search(r'pattern: "([^"]+)"', _block('PANGENE')).group(1)
        assert sorted(_expand(published)) == sorted(declared)
        for name in declared:
            assert (work / name).is_file() and (work / name).stat().st_size > 0, name


# ---------------------------------------------------------------------------
# bin/pangenome.py on its own
# ---------------------------------------------------------------------------

class TestPangenomeScript:

    def _proteins(self, tmp_path, reference, tables):
        for name, text in tables.items():
            (tmp_path / name).write_text(text)
        (tmp_path / 'ref.faa').write_text(reference)
        out = tmp_path / 'proteins.faa'
        pangenome.write_proteins(tmp_path / 'ref.faa',
                                 [tmp_path / n for n in tables if n.endswith('.tsv')], out)
        return [name for name, _ in pangenome.read_fasta(out)]

    def test_no_two_genes_under_one_name(self, tmp_path):
        # A gene of a sample called like a reference gene, and a colon, which
        # pangene would cut the gene's name at
        names = self._proteins(tmp_path, '>K12_00001_1-30\r\nMKV\r\n>a:b_1-30\r\nMRL\r\n', {
            'K12.AMAP.tsv': 'id\tstrand\tstart\tend\trvs\nK12_00001\t+\t1\t30\t\n',
            'K12.AMAP.faa': '>K12_00001;+;0-30\nMTTT*\n',
        })
        assert names == ['K12_00001:K12_00001_1-30', 'a_b:a:b_1-30',
                         'K12_00001_2:K12_00001;+;0-30']

    def test_a_table_without_rvs_stops_it(self, tmp_path):
        with pytest.raises(SystemExit, match='no id and rvs columns'):
            self._proteins(tmp_path, '>Rv0001_1-30\nMKV\n', {
                'S1.AMAP.tsv': 'id\tstrand\tstart\tend\nS1_00001\t+\t1\t30\n',
                'S1.AMAP.faa': '>S1_00001;+;0-30\nMTTT\n'})

    def test_a_graph_with_no_genome(self, tmp_path):
        # gfa2matrix on an empty graph, as pangene leaves it when no protein
        # is found anywhere, measured in the container
        (tmp_path / 'raw.Rtab').write_text('Gene\t\n')
        assert pangenome.write_matrix(tmp_path / 'raw.Rtab', ['genomes/B.fasta', 'genomes/A.fasta'],
                                      tmp_path / 'out.Rtab') == (0, 2, 0)
        assert (tmp_path / 'out.Rtab').read_text() == 'Gene\tA\tB\n'

    def test_haplotypes_of_one_genome_are_one_column(self, tmp_path):
        (tmp_path / 'raw.Rtab').write_text('Gene\tA#0\tA#1\tB#0\nRv0001\t0\t1\t1\nRv0002\t0\t0\t1\n')
        pangenome.write_matrix(tmp_path / 'raw.Rtab', ['A.fasta', 'B.fasta'], tmp_path / 'out.Rtab')
        assert (tmp_path / 'out.Rtab').read_text() == 'Gene\tA\tB\nRv0001\t1\t1\nRv0002\t0\t1\n'


# ---------------------------------------------------------------------------
# The real tools
# ---------------------------------------------------------------------------

REAL_TOOLS = ['miniprot', 'pangene', 'k8', 'gfatools']


@pytest.mark.skipif(not all(shutil.which(t) for t in REAL_TOOLS)
                    or not Path('/usr/local/bin/pangene.js').is_file(),
                    reason='needs miniprot, pangene, k8 and gfatools (the AMAP container has them)')
def test_real_genomes(tmp_path):
    work = _pangene(tmp_path, GENOMES, _tools(tmp_path / 'tools', {}))
    columns, rows = _matrix(work)
    assert columns == SAMPLES
    # RD4, which M. bovis lacks
    assert {gene: rows[gene] for gene in RD4} == {gene: ['0', '1', '1'] for gene in RD4}
    # the genes around RD4 and RvD2, in every genome
    for gene in ['Rv1503c', 'Rv1504c', 'Rv1505c', 'Rv1517', 'Rv1754c', 'Rv1755c', 'Rv1758']:
        assert rows[gene] == ['1', '1', '1'], gene
    # the three genes no reference gene matches, in CDC1551 and M. bovis and
    # not in H37Rv: one row each, under whichever sample's copy pangene kept
    outside = {gene: values for gene, values in rows.items() if not gene.startswith('Rv')}
    assert len(outside) == 3
    assert all(values == ['1', '1', '0'] for values in outside.values()), outside
