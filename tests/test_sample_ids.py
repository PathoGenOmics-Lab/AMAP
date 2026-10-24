#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sample IDs with underscores, in the steps after the correction.

Lineage typing names a genome <lineage>_<sampleId>, L4_H37Rv_test, and under
--skip_typing a genome keeps its sample ID, H37Rv_test. Sample IDs are letters,
digits and underscores, and the steps after the correction told the lineage
from the sample by splitting names at underscores: L4_H37Rv_test counted as
genome H37Rv, the tutorial's Beijing_1 and a Beijing_2 as one genome Beijing,
an untyped H37Rv_test as lineage H37Rv, and each gene of an untyped sample with
no underscore as a genome of its own in the duplication network. PATHOTYPR now
says which lineage it named a genome after, and main.nf hands every step the
sample ID and the lineage.

The tests run the processes' own scripts the way Nextflow renders them, one
after the other as the pipeline chains them, for a typed run and an untyped
one, on files shaped like a real run's: Bakta's table, miniprot's PAF, FIX_CDS's
table and FASTA headers (L4_CDC1551d_00002;+;1616-2390), eggNOG's 21 columns
and Mycobrowser's 34. Diamond is a stand-in that reports identical proteins as
hits, unless the real one is installed. The report's own R code runs where R
is, as in the container.
"""

import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'bin'))

# (the name a genome's files carry, its sample ID, its lineage)
TYPED = [('L2_Beijing_1', 'Beijing_1', 'L2'),
         ('L2_Beijing_2', 'Beijing_2', 'L2'),
         ('L4_H37Rv_test', 'H37Rv_test', 'L4'),
         ('NA_K12_self', 'K12_self', 'NA')]          # typing found no lineage
UNTYPED = [('Beijing_1', 'Beijing_1', 'Untyped'),
           ('Beijing_2', 'Beijing_2', 'Untyped'),
           ('H37Rv', 'H37Rv', 'Untyped')]


def read_tsv(path, **kwargs):
    """A table as written, "NA" being a lineage and not a missing value."""
    return pd.read_csv(path, sep='\t', keep_default_na=False, **kwargs)


# ---------------------------------------------------------------------------
# A process's script, rendered and run the way Nextflow runs it
# ---------------------------------------------------------------------------

def process_script(nf_file, process):
    """A process's script block, with the indentation Nextflow strips."""
    text = (ROOT / 'modules' / 'local' / nf_file).read_text()
    body = re.search(rf'^process {process} \{{\n(.*?)^\}}', text, re.S | re.M).group(1)
    script = re.search(r'^[ \t]*script:\n.*?"""\n(.*?)\n[ \t]*"""', body, re.S | re.M).group(1)
    return textwrap.dedent(script) + '\n'


def render(script, values):
    """Groovy's triple-quoted string, for the little these scripts use: \\\\,
    \\$, ${...} and a backslash ending a line, which Groovy drops with the line
    break. Anything else fails rather than render unlike Nextflow."""
    out, i = [], 0
    while i < len(script):
        c = script[i]
        if c == '\\':
            nxt = script[i + 1]
            if nxt == '\n':
                i += 2
                continue
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
            out.append('' if expr.startswith('task.ext.code') else str(values[expr]))
            i = j
        else:
            out.append(c)
            i += 1
    return ''.join(out)


def run_task(work, nf_file, process, values, inputs, tools):
    """Stage the inputs as links in a task directory, as Nextflow does, and
    run the process's script there with the shell Nextflow gives it."""
    work.mkdir(parents=True)
    for f in inputs:
        (work / Path(f).name).symlink_to(Path(f).resolve())
    script = render(process_script(nf_file, process), dict(values, projectDir=ROOT))
    (work / '.command.sh').write_text(script)
    # Nextflow runs a script without a shebang with `bash -ue`
    shell = ['bash'] if script.startswith('#!') else ['bash', '-ue']
    env = dict(os.environ, PATH=f'{tools}{os.pathsep}{os.environ["PATH"]}',
               PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
    done = subprocess.run(shell + ['.command.sh'], cwd=work, env=env,
                          capture_output=True, text=True)
    assert done.returncode == 0, f'{process} failed:\n{done.stdout[-1500:]}\n{done.stderr[-3000:]}'
    return work


DIAMOND_STAND_IN = '''#!PYTHON
"""diamond: makedb makes nothing, and blastp reports every two identical
proteins of the query file as hits, in --outfmt 6 with AMAP's 15 fields."""
import sys

args = sys.argv[1:]
if args[0] == 'makedb':
    sys.exit(0)
seqs, name = {}, None
for line in open(args[args.index('-q') + 1]):
    line = line.strip()
    if line.startswith('>'):
        name = line[1:].split()[0]
        seqs[name] = ''
    elif name:
        seqs[name] += line
print('# DIAMOND stand-in')
print('# Fields: Query ID, Subject ID, Percentage of identical matches, ...')
for q, qs in seqs.items():
    for s, ss in seqs.items():
        if q != s and qs == ss:
            n = len(qs)
            print('\\t'.join(map(str, [q, s, 100.0, n, 0, 0, 1, n, 1, n, 1e-50, 250.0, n, n, n])))
'''


@pytest.fixture(scope='module')
def tools(tmp_path_factory):
    """python3 as the interpreter running the tests, and diamond as a
    stand-in unless the real one is installed."""
    d = tmp_path_factory.mktemp('tools')
    (d / 'python3').write_text(f'#!/bin/bash\nexec "{sys.executable}" "$@"\n')
    if not shutil.which('diamond'):
        (d / 'diamond').write_text(DIAMOND_STAND_IN.replace('PYTHON', sys.executable))
    for tool in d.iterdir():
        tool.chmod(0o755)
    return d


# ---------------------------------------------------------------------------
# A genome: three reference genes, Rv0003 missed by Bakta, and a gene of no
# reference in a few identical copies, which makes a duplication
# ---------------------------------------------------------------------------

_CODONS = [a + b + c for a in 'ACGT' for b in 'ACGT' for c in 'ACGT'
           if a + b + c not in ('TAA', 'TAG', 'TGA')]
_B, _AAS = 'TCAG', 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
_CODE = {a + b + c: _AAS[16 * i + 4 * j + k]
         for i, a in enumerate(_B) for j, b in enumerate(_B) for k, c in enumerate(_B)}
_rng = random.Random(7)
GENES = {key: 'ATG' + ''.join(_rng.choice(_CODONS) for _ in range(n - 2)) + 'TAA'
         for key, n in [('Rv0001', 100), ('Rv0002', 120), ('Rv0003', 90), ('novel', 80)]}
SPACERS = [''.join(_rng.choice('ACGT') for _ in range(150)) for _ in range(12)]
COG = {'Rv0001': 'L', 'Rv0002': 'S', 'Rv0003': 'K', 'novel': 'X'}
CATEGORY = {'Rv0001': 'information pathways', 'Rv0002': 'conserved hypotheticals',
            'Rv0003': 'PE/PPE'}


def _rc(nt):
    return nt[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def _protein(nt):
    return ''.join(_CODE[nt[i:i + 3]] for i in range(0, len(nt), 3))


def genome(copies):
    """The sequence, and its genes as (gene, start, end, strand), 1-based."""
    seq, genes = SPACERS[0], []
    order = [('Rv0001', '+'), ('Rv0002', '-'), ('Rv0003', '+')] + [('novel', '+')] * copies
    for n, (key, strand) in enumerate(order, 1):
        start = len(seq) + 1
        seq += GENES[key] if strand == '+' else _rc(GENES[key])
        genes.append((key, start, len(seq), strand))
        seq += SPACERS[n]
    return seq, genes


def write_sample(d, name, copies):
    """What TYPE_CASES and GET_OUTPUTS are handed for a sample, as the
    pipeline names it: its genome, a one-sequence genome being renamed after
    the sample, Bakta's table and miniprot's alignments, FIX_CDS's table and
    FASTA files, and eggNOG's annotation of them."""
    d.mkdir(parents=True)
    seq, genes = genome(copies)
    contig, length = name, len(seq)
    (d / f'{name}.fasta').write_text(f'>{contig}\n{seq}\n')

    # Bakta calls everything but Rv0003; locus tags are <name>_00001 on
    tags = {}
    with open(d / f'{name}.tsv', 'w') as f:
        f.write('# Annotated with Bakta\n# Software: v1.12.0\n'
                '#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n')
        for i, (key, s, e, st) in enumerate(g for g in genes if g[0] != 'Rv0003'):
            tags[s] = f'{name}_{i + 1:05d}'
            f.write(f'{contig}\tcds\t{s}\t{e}\t{st}\t{tags[s]}\t\thypothetical protein\t\n')

    # miniprot places the three reference proteins, NAME_START-END
    ref = [g for g in genes if g[0].startswith('Rv')]
    with open(d / 'ref.faa', 'w') as f:
        for key, s, e, st in ref:
            f.write(f'>{key}_{s}-{e}\n{_protein(GENES[key])}\n')
    with open(d / f'{name}.paf', 'w') as f:
        for key, s, e, st in ref:
            n, nt = len(GENES[key]) // 3, e - s + 1
            f.write(f'{key}_{s}-{e}\t{n}\t0\t{n}\t{st}\t{contig}\t{length}\t{s - 1}\t{e}\t{nt}\t{nt}\t0\t'
                    f'AS:i:500\tms:i:500\tnp:i:{n}\tfs:i:0\tst:i:1\tda:i:0\tdo:i:0\tcg:Z:{n}M\tcs:Z::{n}\n')
    (d / f'{name}-badbakta.paf').write_text('')
    # one contig each, so no contig end cuts a gene: MINIPROT's second pass places nothing
    (d / f'{name}-pieces.paf').write_text('')
    # Bakta's proteins of reference genes are all in the reference genome
    with open(d / f'{name}-ref.paf', 'w') as f:
        for key, s, e, st in ref:
            if s in tags:
                n, nt = (e - s + 1) // 3, e - s + 1
                f.write(f'{tags[s]}\t{n}\t0\t{n}\t{st}\tref\t{length}\t{s - 1}\t{e}\t{nt}\t{nt}\t0\tst:i:1\n')

    # FIX_CDS: Bakta's locus tags where Bakta called the gene, the reference
    # name where AMAP recovered it; the FASTA headers carry ID;strand;coords
    rows = [(tags.get(s, f'{name}_{key}'), st, s, e, key if key.startswith('Rv') else '', key)
            for key, s, e, st in genes]
    with open(d / f'{name}.AMAP.tsv', 'w') as f:
        f.write('id\tstrand\tstart\tend\trvs\tpseudogene\tcontig\tpartial\tnote\n')
        for gid, st, s, e, rvs, key in rows:
            f.write(f'{gid}\t{st}\t{s}\t{e}\t{rvs}\tno\t{contig}\t\t\n')
    with open(d / f'{name}.AMAP.faa', 'w') as faa, open(d / f'{name}.AMAP.fna', 'w') as fna:
        for gid, st, s, e, rvs, key in rows:
            faa.write(f'>{gid};{st};{s - 1}-{e}\n{_protein(GENES[key]).rstrip("*")}\n')
            fna.write(f'>{gid};{st};{s - 1}-{e}\n{GENES[key]}\n')
    with open(d / f'{name}.emapper.annotations', 'w') as f:
        f.write('## emapper-2.1.12\n##\n#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\t'
                'max_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\t'
                'KEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\t'
                'BiGG_Reaction\tPFAMs\n')
        for gid, st, s, e, rvs, key in rows:
            seed = f'83332.{key}' if rvs else '-'
            f.write('\t'.join([f'{gid};{st};{s - 1}-{e}', seed, '1.5e-60', '230.0', '-', '-',
                               COG[key], '-', '-'] + ['-'] * 12) + '\n')
    return d


def write_mycobrowser(path):
    """Mycobrowser's table, its 34 columns, for the three reference genes."""
    header = ['Refseq_ID', 'Mycobrowser_Version', 'Feature', 'Start', 'Stop', 'Score', 'Strand',
              'Frame', 'Locus', 'Name', 'Function', 'Product', 'Comments', 'UniProt_AC',
              'Is_Pseudogene', 'Functional_Category'] + [f'extra{i}' for i in range(18)]
    with open(path, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for key, category in CATEGORY.items():
            row = ['NC_000962.3', 'Mycobrowser_v4', 'CDS', '1', '300', '.', '+', '0', key, key,
                   'Unknown', 'Hypothetical protein', '', '', 'No', category] + [''] * 18
            f.write('\t'.join(row) + '\n')


# ---------------------------------------------------------------------------
# The steps after the correction, chained as main.nf chains them
# ---------------------------------------------------------------------------

def run_pipeline(base, samples, tools):
    """TYPE_CASES, GET_OUTPUTS, GET_COUNT and GET_DUPLICATIONS for every
    sample, then REPORT_INPUTS and GET_ALLDUP for them all. Returns the task
    directories."""
    myco = base / 'mycobrowser.txt'
    write_mycobrowser(myco)
    tasks = {'cases': {}, 'outputs': {}, 'count': {}, 'report_inputs': [], 'alldup': []}
    for k, (name, sample, lineage) in enumerate(samples):
        d = write_sample(base / 'inputs' / name, name, copies=2 + k)
        cases = run_task(base / 'TYPE_CASES' / name, 'cases.nf', 'TYPE_CASES', {
            'sampleId': name, 'sample': sample, 'lineage': lineage,
            'miniprot_paf': f'{name}.paf', 'badbakta_paf': f'{name}-badbakta.paf',
            'h37rv_paf': f'{name}-ref.paf', 'pieces_paf': f'{name}-pieces.paf',
            'fasta': f'{name}.fasta', 'bakta_tsv': f'{name}.tsv', 'ref_proteome': 'ref.faa',
        }, [d / f for f in (f'{name}.paf', f'{name}-badbakta.paf', f'{name}-ref.paf',
                            f'{name}-pieces.paf', f'{name}.fasta', f'{name}.tsv', 'ref.faa')], tools)

        downstream = [d / f for f in (f'{name}.AMAP.tsv', f'{name}.AMAP.faa', f'{name}.AMAP.fna',
                                      f'{name}.tsv', f'{name}.fasta', f'{name}.emapper.annotations')]
        files = {'sampleId': name, 'amap_tsv': f'{name}.AMAP.tsv', 'amap_faa': f'{name}.AMAP.faa',
                 'amap_fna': f'{name}.AMAP.fna', 'bakta_tsv': f'{name}.tsv',
                 'fasta': f'{name}.fasta', 'emapper': f'{name}.emapper.annotations'}
        outputs = run_task(base / 'GET_OUTPUTS' / name, 'get.nf', 'GET_OUTPUTS', dict(
            files, lineage=lineage, **{'myco_ref': myco,
                                       'params.skip_mycobrowser': 'false'}), downstream, tools)
        master = outputs / f'{name}-AMAPmastertable.tsv'
        count = run_task(base / 'GET_COUNT' / name, 'get.nf', 'GET_COUNT', {
            'sampleId': name, 'mastertable': master.name, 'sample': sample, 'lineage': lineage,
        }, [master], tools)
        dups = run_task(base / 'GET_DUPLICATIONS' / name, 'get.nf', 'GET_DUPLICATIONS',
                        dict(files, **{'task.cpus': 1}), downstream, tools)

        tasks['report_inputs'] += [p for p in cases.iterdir() if p.name.endswith(
            ('-summary.tsv', '-rvmissing.tsv'))]
        tasks['report_inputs'] += [p for p in count.iterdir() if p.suffix == '.tsv'
                                   and not p.is_symlink()]
        tasks['report_inputs'] += [master, outputs / f'{name}-AMAP.qc.tsv']
        tasks['alldup'] += [p for p in dups.iterdir() if p.suffix == '.dup']
        tasks['cases'][name], tasks['outputs'][name], tasks['count'][name] = cases, outputs, count
    # PANGENE's matrix, one for the run: a column per genome, named after it.
    # Every genome has the three reference genes and the same gene of no
    # reference, which goes by the copy of the first genome in name order.
    names = sorted(name for name, _sample, _lineage in samples)
    first = read_tsv(base / 'inputs' / names[0] / f'{names[0]}.AMAP.tsv')
    genes = ['Rv0001', 'Rv0002', 'Rv0003', first[first['rvs'] == '']['id'].iloc[0]]
    lines = ['\t'.join(['Gene'] + names)] + ['\t'.join([gene] + ['1'] * len(names)) for gene in genes]
    rtab = base / 'inputs' / 'pangenome.Rtab'
    rtab.write_text('\n'.join(lines) + '\n')
    tasks['report_inputs'].append(rtab)
    tasks['report'] = run_task(base / 'REPORT_INPUTS', 'freport.nf', 'REPORT_INPUTS', {},
                               tasks['report_inputs'], tools)
    tasks['dup'] = run_task(base / 'GET_ALLDUP', 'get.nf', 'GET_ALLDUP',
                            {'task.memory': '8 GB', 'task.cpus': 1}, tasks['alldup'], tools)
    return tasks


@pytest.fixture(scope='module')
def typed(tmp_path_factory, tools):
    return run_pipeline(tmp_path_factory.mktemp('typed'), TYPED, tools)


@pytest.fixture(scope='module')
def untyped(tmp_path_factory, tools):
    return run_pipeline(tmp_path_factory.mktemp('untyped'), UNTYPED, tools)


SCENARIOS = [('typed', TYPED), ('untyped', UNTYPED)]


@pytest.fixture(params=SCENARIOS, ids=[s for s, _ in SCENARIOS])
def run(request):
    """(the tasks of a run, its samples) for the typed run and the untyped one."""
    scenario, samples = request.param
    return request.getfixturevalue(scenario), samples


class TestEveryStepKnowsTheSample:

    def test_case_summary(self, run):
        tasks, samples = run
        for name, sample, lineage in samples:
            summary = read_tsv(tasks['cases'][name] / f'{name}-summary.tsv')
            assert set(summary['Sequence']) == {name}
            assert set(summary['Lineage']) == {lineage}
            assert set(summary['Sample']) == {sample}

    def test_master_table_lineage(self, run):
        tasks, samples = run
        for name, sample, lineage in samples:
            master = read_tsv(tasks['outputs'][name] / f'{name}-AMAPmastertable.tsv')
            assert set(master['lineage']) == {lineage}

    def test_category_counts(self, run):
        tasks, samples = run
        for k, (name, sample, lineage) in enumerate(samples):
            for table in ('COG', 'MycoBrowser'):
                counts = read_tsv(tasks['count'][name] / f'{name}-AMAP.{table}.tsv')
                assert set(zip(counts['lineage'], counts['genome'])) == {(lineage, sample)}, table
            cog = read_tsv(tasks['count'][name] / f'{name}-AMAP.COG.tsv')
            # every copy of the gene of no reference is X
            assert cog.loc[cog['COG_category'] == 'X', 'count'].tolist() == [2 + k]

    def test_genes_not_in_the_reference(self, run):
        tasks, samples = run
        for k, (name, sample, lineage) in enumerate(samples):
            noref = read_tsv(tasks['count'][name] / f'{name}-AMAP.noRef.tsv')
            assert noref.columns.tolist() == ['genome', 'lineage']
            assert noref.values.tolist() == [[sample, lineage]] * (2 + k)

    def test_one_genome_per_sample_in_the_duplication_network(self, run):
        tasks, samples = run
        communities = read_tsv(tasks['dup'] / 'homologous_communities.tsv')
        sizes = {}
        for genome, size in list(zip(communities['Genome1'], communities['Size1'])) + \
                list(zip(communities['Genome2'], communities['Size2'])):
            sizes.setdefault(genome, set()).add(size)
        # each genome's copies of the duplicated gene make one community of its own
        assert sizes == {name: {2 + k} for k, (name, _, _) in enumerate(samples)}
        assert len(communities) == len(samples) * (len(samples) - 1) // 2


# ---------------------------------------------------------------------------
# The report's own R code, on what REPORT_INPUTS and GET_ALLDUP made
# ---------------------------------------------------------------------------

REPORT = ROOT / 'bin' / 'report' / 'report_AMAP.Rmd'
R_PACKAGES = ['ggplot2', 'plotly', 'dplyr', 'tidyr', 'ggpubr', 'DT', 'stringr', 'knitr',
              'rstatix', 'ggrepel', 'scales', 'pheatmap', 'RColorBrewer']
R_CHUNKS = ['setup', 'sample_quality', 'cds_NoH37Rv', 'rvmissing_table', 'statsCOGs',
            'dup_table', 'peppe']
R_TABLES = {'data': 'data', 'quality_data': 'quality', 'datos_cog': 'cog', 'datos_myco': 'myco',
            'resultado': 'noref', 'rvmissing_raw': 'rvmissing', 'dup_display': 'dup',
            'peppe_count': 'peppe'}


def _report_chunks():
    """The R code of the report's chunks this needs, by label."""
    text = REPORT.read_text()
    chunks = dict(re.findall(r'^```\{r (\w+)[^}]*\}\n(.*?)^```', text, re.S | re.M))
    return [chunks[label] for label in R_CHUNKS]


@pytest.fixture(scope='module')
def rscript():
    exe = shutil.which('Rscript')
    if not exe:
        pytest.skip('R is not installed (the container has it)')
    libs = '; '.join(f'suppressMessages(library({p}))' for p in R_PACKAGES)
    if subprocess.run([exe, '-e', libs], capture_output=True).returncode:
        pytest.skip("the report's R packages are not installed (the container has them)")
    return exe


def stage_report_inputs(tasks, work, replace):
    """A run's report inputs, REPORT_INPUTS' tables and GET_ALLDUP's
    communities, in the folder `work`, those named in `replace` holding the
    text given there instead. Returns their paths."""
    staged = []
    for f in list(tasks['report'].glob('input_*.txt')) + [tasks['dup'] / 'homologous_communities.tsv']:
        if f.name in replace:
            (work / f.name).write_text(replace[f.name])
        else:
            (work / f.name).symlink_to(f)
        staged.append(work / f.name)
    return staged


def report_tables(tasks, rscript, empty=(), replace=None):
    """Run the report's chunks on a run's report inputs, the way FINAL_REPORT
    hands them over, those named in `empty` left empty and those in `replace`
    holding the text given there, and read back the tables they build."""
    work = Path(tempfile.mkdtemp(prefix='FINAL_REPORT-', dir=tasks['report'].parent))
    stage_report_inputs(tasks, work, {**dict.fromkeys(empty, ''), **(replace or {})})
    names = ['cog', 'myco', 'summary', 'noref', 'rvmissing', 'mastertables', 'pangene', 'qc']
    params = ', '.join([f'{"empty" if n == "noref" else "mastertable" if n == "mastertables" else n}'
                        f'_file="input_{n}.txt"' for n in names]
                       + ['dup_file="homologous_communities.tsv"', 'version="test"'])
    for i, code in enumerate(_report_chunks()):
        (work / f'chunk{i}.R').write_text(code)
    dumps = '\n'.join(f'write.table({var}, "out_{out}.tsv", sep="\\t", quote=FALSE, row.names=FALSE)'
                      for var, out in R_TABLES.items())
    (work / 'run.R').write_text(
        f'params <- list({params})\n'
        f'for (f in sprintf("chunk%d.R", 0:{len(R_CHUNKS) - 1})) source(f, print.eval=FALSE)\n'
        f'{dumps}\n')
    done = subprocess.run([rscript, 'run.R'], cwd=work, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr[-3000:]
    return {out: read_tsv(work / f'out_{out}.tsv') for out in R_TABLES.values()}


@pytest.fixture(scope='module')
def typed_report(typed, rscript):
    return report_tables(typed, rscript)


@pytest.fixture(scope='module')
def untyped_report(untyped, rscript):
    return report_tables(untyped, rscript)


@pytest.fixture(params=SCENARIOS, ids=[s for s, _ in SCENARIOS])
def report(request):
    """(the report's tables, the run's samples) for both runs."""
    scenario, samples = request.param
    return request.getfixturevalue(f'{scenario}_report'), samples


class TestReport:

    def test_cases_by_sample(self, report):
        tables, samples = report
        data = tables['data']
        assert set(zip(data['V3'], data['V4'], data['genome'])) == \
            {(name, lineage, sample) for name, sample, lineage in samples}

    def test_sample_quality_has_a_row_per_sample(self, report):
        tables, samples = report
        quality = tables['quality']
        assert sorted(zip(quality['Lineage'], quality['Genome'])) == \
            sorted((lineage, sample) for _, sample, lineage in samples)

    def test_category_percentages_are_per_sample(self, report):
        tables, samples = report
        for table in ('cog', 'myco'):
            counts = tables[table]
            assert set(zip(counts['Lineages'], counts['genome'])) == \
                {(lineage, sample) for _, sample, lineage in samples}, table
            # each sample's own categories add up to all of it, not half
            for genome, pct in counts.groupby('genome')['Percentage'].sum().items():
                assert pct == pytest.approx(1.0), (table, genome)

    def test_genes_not_in_the_reference_per_sample(self, report):
        tables, samples = report
        noref = tables['noref']
        assert sorted(zip(noref['genome'], noref['V2'], noref['n'])) == \
            sorted((sample, lineage, 2 + k) for k, (_, sample, lineage) in enumerate(samples))

    def test_every_gene_in_the_reference(self, typed, rscript):
        # As when H37Rv is annotated against itself: REPORT_INPUTS hands over
        # an empty table, whose columns R reads as logical, and cutting its
        # genomes at underscores stopped the report. Every sample has none.
        tables = report_tables(typed, rscript, empty=['input_noref.txt'])
        noref = tables['noref']
        assert sorted(zip(noref['genome'], noref['V2'], noref['n'])) == \
            sorted((sample, lineage, 0) for _, sample, lineage in TYPED)

    def test_a_sample_with_every_gene_in_the_reference_counts_zero(self, typed, rscript):
        # It stays in its lineage's box at 0 rather than drop out of it, which
        # moved the box up to the samples that have some
        _, none, _ = TYPED[0]
        lines = (typed['report'] / 'input_noref.txt').read_text().splitlines(keepends=True)
        others = ''.join(line for line in lines if line.split('\t')[0] != none)
        tables = report_tables(typed, rscript, replace={'input_noref.txt': others})
        noref = tables['noref']
        assert sorted(zip(noref['genome'], noref['V2'], noref['n'])) == \
            sorted((sample, lineage, 0 if sample == none else 2 + k)
                   for k, (_, sample, lineage) in enumerate(TYPED))

    def test_missing_genes_by_lineage(self, report):
        tables, samples = report
        rv = tables['rvmissing']
        assert set(zip(rv['Rv_gene'], rv['Sample'], rv['Lineage'])) == \
            {('Rv0003', name, lineage) for name, _, lineage in samples}

    def test_duplication_communities_by_genome(self, report):
        tables, samples = report
        dup = tables['dup']
        assert set(dup['Genome1']) | set(dup['Genome2']) == {name for name, _, _ in samples}

    def test_master_table_by_lineage(self, report):
        tables, samples = report
        assert set(tables['peppe']['lineage']) == {lineage for _, _, lineage in samples}


# ---------------------------------------------------------------------------
# The whole report, rendered by FINAL_REPORT's own script
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def pandoc(rscript):
    """rmarkdown renders the report with pandoc (the container has both)."""
    if subprocess.run([rscript, '-e', 'quit(status = !rmarkdown::pandoc_available())'],
                      capture_output=True).returncode:
        pytest.skip('rmarkdown or pandoc is not installed (the container has them)')


def final_report(tasks, tools, replace=None):
    """Run FINAL_REPORT's script on a run's report inputs, those named in
    `replace` holding the text given there, and return the report's HTML."""
    inputs = Path(tempfile.mkdtemp(prefix='inputs-', dir=tasks['report'].parent))
    staged = stage_report_inputs(tasks, inputs, replace or {})
    work = inputs.with_name(inputs.name.replace('inputs-', 'FINAL_REPORT-task-'))
    tables = {f'input_{n}': f'input_{n}.txt' for n in
              ['cog', 'myco', 'summary', 'noref', 'rvmissing', 'mastertables', 'pangene', 'qc']}
    run_task(work, 'freport.nf', 'FINAL_REPORT', dict(
        tables, dup_communities_tsv='homologous_communities.tsv',
        **{'workflow.manifest.version': 'test'}),
        staged + [tasks['dup'] / 'homologous_communities.html'], tools)
    return (work / 'final_report.html').read_text(encoding='utf-8')


def section(html, title):
    """A section of the report, from its heading to the next one."""
    start = html.index(f'>{title}</h2>')
    return html[start:html.index('<h2', start)]


class TestFinalReport:
    """Each renders the whole report, figures and all, so there are two."""

    def test_a_run(self, typed, tools, pandoc):
        assert '<img' in section(final_report(typed, tools), 'CDS not in Reference')

    def test_a_run_that_leaves_every_table_it_can_empty(self, typed, tools, pandoc):
        # Every gene in the reference (H37Rv against itself), no reference gene
        # missing, no Mycobrowser categories (--skip_mycobrowser), no pangenome
        # to merge, no community across genomes. The first of these stopped
        # FINAL_REPORT, and the run with it.
        communities = (typed['dup'] / 'homologous_communities.tsv').read_text()
        html = final_report(typed, tools, {
            'input_noref.txt': '', 'input_rvmissing.txt': '', 'input_myco.txt': '',
            'input_pangene.txt': '', 'homologous_communities.tsv': communities.splitlines(True)[0]})
        noref = section(html, 'CDS not in Reference')
        assert 'Every CDS of every sample matches a reference gene.' in noref
        assert '<img' not in noref


# ---------------------------------------------------------------------------
# The scripts on their own
# ---------------------------------------------------------------------------

class TestWriteSummary:

    def _summary(self, tmp_path, monkeypatch, name, **kwargs):
        from classify.io import write_summary
        monkeypatch.chdir(tmp_path)
        write_summary(name, [pd.DataFrame({'x': [1, 2]})], ['case2'], **kwargs)
        return read_tsv(tmp_path / f'{name}-summary.tsv')

    def test_typed(self, tmp_path, monkeypatch):
        summary = self._summary(tmp_path, monkeypatch, 'L4_H37Rv_test', sample='H37Rv_test', lineage='L4')
        assert summary.values.tolist() == [['case2', 2, 'L4_H37Rv_test', 'L4', 'H37Rv_test']]

    def test_untyped(self, tmp_path, monkeypatch):
        summary = self._summary(tmp_path, monkeypatch, 'H37Rv_test', sample='H37Rv_test', lineage='Untyped')
        assert summary.values.tolist() == [['case2', 2, 'H37Rv_test', 'Untyped', 'H37Rv_test']]

    def test_a_lineage_with_an_underscore(self, tmp_path, monkeypatch):
        # pathotypr's deeper calls hold underscores (M_bovis, M_caprae)
        summary = self._summary(tmp_path, monkeypatch, 'M_bovis_Mb_1', sample='Mb_1', lineage='M_bovis')
        assert summary[['Lineage', 'Sample']].values.tolist() == [['M_bovis', 'Mb_1']]

    def test_told_nothing_the_sample_is_its_name_and_untyped(self, tmp_path, monkeypatch):
        summary = self._summary(tmp_path, monkeypatch, 'L4_H37Rv_test')
        assert summary[['Sequence', 'Lineage', 'Sample']].values.tolist() == \
            [['L4_H37Rv_test', 'Untyped', 'L4_H37Rv_test']]


class TestNetworkGenomes:

    def _blast(self, pairs):
        rows = [[q, s, 100.0, 100, 0, 0, 1, 100, 1, 100, 1e-50, 200.0, 100, 100, 100]
                for q, s in pairs]
        df = pd.DataFrame(rows, columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gaps',
                                         'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
                                         'qlen', 'slen', 'nident'])
        return df

    def test_a_genome_named_like_the_start_of_another_keeps_its_own_genes(self):
        from network_genomes import build_intra_genome_networks, get_genome_ids
        genome_of = {'L2_Beijing_1_00003;X': 'L2_Beijing_1', 'L2_Beijing_1_00004;X': 'L2_Beijing_1',
                     'L2_Beijing_10_00003;X': 'L2_Beijing_10', 'L2_Beijing_10_00004;X': 'L2_Beijing_10'}
        df = self._blast([('L2_Beijing_1_00003;X', 'L2_Beijing_1_00004;X'),
                          ('L2_Beijing_10_00003;X', 'L2_Beijing_10_00004;X'),
                          ('L2_Beijing_1_00003;X', 'L2_Beijing_10_00003;X')])
        genomes = get_genome_ids(df, genome_of)
        assert genomes == {'L2_Beijing_1', 'L2_Beijing_10'}
        networks = build_intra_genome_networks(df, genomes, genome_of)
        assert set(networks['L2_Beijing_1'].nodes) == {'L2_Beijing_1_00003;X', 'L2_Beijing_1_00004;X'}
        assert set(networks['L2_Beijing_10'].nodes) == {'L2_Beijing_10_00003;X', 'L2_Beijing_10_00004;X'}

    def test_a_sequence_in_no_genome_stops_it(self):
        from network_genomes import get_genome_ids
        df = self._blast([('a;X', 'b;X')])
        with pytest.raises(ValueError, match='no genome'):
            get_genome_ids(df, {'a;X': 'S1'})

    def test_genome_table(self, tmp_path):
        from network_genomes import read_genomes
        table = tmp_path / 'genomes.tsv'
        table.write_text('L2_Beijing_1_00003;X\tL2_Beijing_1\nH37Rv_00004;+;10-250\tH37Rv\n')
        assert read_genomes(table) == {'L2_Beijing_1_00003;X': 'L2_Beijing_1',
                                       'H37Rv_00004;+;10-250': 'H37Rv'}


# ---------------------------------------------------------------------------
# PATHOTYPR says which lineage it named the genome after
# ---------------------------------------------------------------------------

PATHOTYPR_STAND_IN = """#!/bin/bash
# pathotypr: copies $FAKE_SUMMARY to <prefix>_summary.tsv, or fails without one
while [ $# -gt 0 ]; do
    case "$1" in
        --output-prefix) prefix="$2"; shift ;;
    esac
    shift
done
[ -n "$FAKE_SUMMARY" ] || exit 1
cp "$FAKE_SUMMARY" "${prefix}_summary.tsv"
"""
SUMMARY_HEADER = 'genome\tlineage:count\tmajor_lineage\n'


def _pathotypr(tmp_path, sample, summary=None):
    """(SAMPLE_NAME, LINEAGE) as PATHOTYPR exports them for a one-sequence
    genome whose header is its sample ID."""
    work = tmp_path / 'work'
    work.mkdir(parents=True)
    genome = tmp_path / f'{sample}.fasta'
    genome.write_text(f'>{sample}\nACGT\n')
    tools = tmp_path / 'tools'
    tools.mkdir()
    (tools / 'pathotypr').write_text(PATHOTYPR_STAND_IN)
    (tools / 'python3').write_text(f'#!/bin/bash\nexec "{sys.executable}" "$@"\n')
    for tool in tools.iterdir():
        tool.chmod(0o755)
    # where PATHOTYPR stages its input, wherever that is
    stage = re.search(r"path\(fasta, stageAs: '([^']+)'\)", (ROOT / 'modules/local/pathotypr.nf').read_text())
    staged = work / (stage.group(1).replace('*', genome.name) if stage else genome.name)
    staged.parent.mkdir(parents=True, exist_ok=True)
    staged.symlink_to(genome)
    script = render(process_script('pathotypr.nf', 'PATHOTYPR'), {
        'sampleId': sample, 'fasta': staged.relative_to(work).as_posix(), 'projectDir': ROOT,
        'markers': 'snps_lineages.tsv',
        'markers_ref': 'MTB_ancestor_reference.fasta'})
    # Nextflow reads env() outputs after the script, as this does
    (work / '.command.sh').write_text(
        script + 'printf %s "$SAMPLE_NAME" > .name\nprintf %s "${LINEAGE-}" > .lineage\n')
    env = dict(os.environ, PATH=f'{tools}{os.pathsep}{os.environ["PATH"]}')
    env.pop('FAKE_SUMMARY', None)
    if summary is not None:
        (tmp_path / 'summary.tsv').write_text(summary)
        env['FAKE_SUMMARY'] = str(tmp_path / 'summary.tsv')
    subprocess.run(['bash', '.command.sh'], cwd=work, env=env, check=True, capture_output=True)
    return (work / '.name').read_text(), (work / '.lineage').read_text()


class TestPathotyprLineage:

    def test_lineage_of_a_typed_genome(self, tmp_path):
        summary = SUMMARY_HEADER + 'H37Rv_test\tL4:79 L4;L4.9:12\tL4\n'
        assert _pathotypr(tmp_path, 'H37Rv_test', summary) == ('L4_H37Rv_test', 'L4')

    def test_a_lineage_with_an_underscore(self, tmp_path):
        # what the name puts first is not always cut off by the first underscore
        summary = SUMMARY_HEADER + 'Mb_1\tA4:40 A4;M_bovis:30\tM_bovis\n'
        assert _pathotypr(tmp_path, 'Mb_1', summary) == ('M_bovis_Mb_1', 'M_bovis')

    def test_no_lineage_is_NA(self, tmp_path):
        assert _pathotypr(tmp_path, 'K12_self') == ('NA_K12_self', 'NA')
        assert _pathotypr(tmp_path / 'b', 'K12_self', SUMMARY_HEADER)[1] == 'NA'
