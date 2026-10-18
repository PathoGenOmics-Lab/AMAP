#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The clean proteome leaves out the pseudogenes AMAP called.

export_clean_proteome.py took the pseudogenes to leave out of
`*-AMAP.clean.faa` from the first column the master table has of
Is_Pseudogene, pseudogene and is_pseudogene, and every master table has
Mycobrowser's Is_Pseudogene, empty when Mycobrowser is off. It then looked
each protein up by its FASTA ID, which FIX_CDS writes as ID;strand;start-end
where the master table has the ID alone, and never found one. Only proteins
with an internal stop codon were left out, so the pseudogenes without one
stayed in: 42 of the 123 AMAP calls in CDC1551 (a frameshift or lost stop
codon, no stop codon within reach, no stop codon as in the reference), and 27
of the 129 in E. coli Sakai against K-12. The clean proteome now leaves out the
genes whose `pseudogene` is yes, those the GFF3 marks pseudo=true and the
GenBank /pseudo, matched on the ID before the ';'.

The tests run the steps after the correction as Nextflow renders them, chained
as tests/test_sample_ids.py chains them, on its genomes with FIX_CDS calling
some of their genes pseudogenes, one with a stop codon in its protein and the
others without, and Mycobrowser marking a gene AMAP finds intact, with
Mycobrowser on and off.
"""

import sys
from pathlib import Path

import pytest
from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'bin'))

from tests import test_sample_ids as chain  # noqa: E402
from tests.test_sample_ids import tools  # noqa: E402,F401  (fixture)

SAMPLES = [('L2_Beijing_1', 'Beijing_1', 'L2'), ('L4_H37Rv_test', 'H37Rv_test', 'L4')]

# Each genome has Rv0001, Rv0002, Rv0003 and 2 + k copies of a gene of no
# reference, k its place. FIX_CDS calls Rv0002 a pseudogene for a frameshift,
# Rv0003, which Bakta missed, for a stop codon inside it, and every copy of the
# other gene but the first for having no stop codon within reach: the reasons
# of CDC1551's pseudogenes. The k-th sample has 3 + k, and only Rv0003, whose
# protein has a stop, was left out of the clean proteome before.
REASONS = {'Rv0002': 'frameshift or lost stop codon: its frame runs on 40 residues '
                     'past the end of the reference protein',
           'Rv0003': '1 internal stop codon(s)',
           'novel': 'no stop codon in its frame within 300 bp'}
STOP = 30   # the codon of Rv0003 that is a stop


def pseudogenes(k):
    """How many pseudogenes AMAP calls in the k-th sample."""
    return 3 + k


# tests/test_sample_ids.py's writers, kept before the runs patch them
_write_sample, _write_mycobrowser = chain.write_sample, chain.write_mycobrowser


def write_sample(d, name, copies):
    """tests/test_sample_ids.py's sample, with FIX_CDS's pseudogenes and their
    reasons in its table, and Rv0003's stop codon in its protein and its gene,
    as FIX_CDS writes them."""
    _write_sample(d, name, copies)
    table = d / f'{name}.AMAP.tsv'
    header, *rows = table.read_text().splitlines()
    novel = 0
    for i, row in enumerate(rows):
        cells = row.split('\t')
        key = cells[4] or 'novel'
        novel += key == 'novel'
        if key in REASONS and not (key == 'novel' and novel == 1):
            cells[5], cells[8] = 'yes', REASONS[key]
            rows[i] = '\t'.join(cells)
    table.write_text('\n'.join([header] + rows) + '\n')
    for fasta, stop in ((d / f'{name}.AMAP.faa', '*'), (d / f'{name}.AMAP.fna', 'TAA')):
        lines = fasta.read_text().splitlines()
        at = next(i for i, line in enumerate(lines) if line.startswith(f'>{name}_Rv0003;')) + 1
        n = len(stop)
        lines[at] = lines[at][:STOP * n] + stop + lines[at][(STOP + 1) * n:]
        fasta.write_text('\n'.join(lines) + '\n')
    return d


def write_mycobrowser(path):
    """Mycobrowser's table, which marks Rv0001 as a pseudogene of H37Rv, a gene
    AMAP finds intact in every sample here, as it finds 6 of the 13 in CDC1551."""
    _write_mycobrowser(path)
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines):
        cells = line.split('\t')
        if cells[8] == 'Rv0001':
            cells[14] = 'Yes'
            lines[i] = '\t'.join(cells)
    path.write_text('\n'.join(lines) + '\n')


def skipping_mycobrowser(run_task):
    """run_task, with GET_OUTPUTS run under --skip_mycobrowser, as a run with a
    reference from another genus is."""
    def run(work, nf_file, process, values, inputs, tools):
        if process == 'GET_OUTPUTS':
            values = dict(values, **{'params.skip_mycobrowser': 'true'})
        return run_task(work, nf_file, process, values, inputs, tools)
    return run


RUNS = ['mycobrowser', 'no_mycobrowser']


@pytest.fixture(scope='module', params=RUNS)
def run(request, tmp_path_factory, tools):
    """(the tasks of a run, whether Mycobrowser was on) for both runs."""
    mycobrowser = request.param == 'mycobrowser'
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(chain, 'write_sample', write_sample)
        mp.setattr(chain, 'write_mycobrowser', write_mycobrowser)
        if not mycobrowser:
            mp.setattr(chain, 'run_task', skipping_mycobrowser(chain.run_task))
        tasks = chain.run_pipeline(tmp_path_factory.mktemp(request.param), SAMPLES, tools)
    return tasks, mycobrowser


def records(path):
    """A FASTA file's (header, sequence) pairs, in its order."""
    out = []
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            out.append([line[1:], ''])
        elif out:
            out[-1][1] += line
    return [tuple(record) for record in out]


def ids(pairs):
    """The IDs of FASTA records, as FIX_CDS writes them before ';'."""
    return [header.split(';')[0] for header, _ in pairs]


def clean(tasks, name):
    return records(tasks['outputs'][name] / f'{name}-AMAP.clean.faa')


def called(tasks, name):
    """The IDs FIX_CDS calls pseudogenes in a sample's table, as GET_OUTPUTS
    was handed it."""
    table = tasks['outputs'][name].parent.parent / 'inputs' / name / f'{name}.AMAP.tsv'
    table = chain.read_tsv(table)
    return set(table.loc[table['pseudogene'] == 'yes', 'id'])


class TestGetOutputs:

    def test_the_master_table_has_both_columns(self, run):
        # as every master table has: Mycobrowser's, empty when it is off, and AMAP's
        tasks, mycobrowser = run
        for k, (name, _, _) in enumerate(SAMPLES):
            master = chain.read_tsv(tasks['outputs'][name] / f'{name}-AMAPmastertable.tsv')
            marked = master.loc[master['Is_Pseudogene'] != '', 'Locus']
            assert dict(zip(marked, master.loc[marked.index, 'Is_Pseudogene'])) == \
                ({'Rv0001': 'Yes', 'Rv0002': 'No', 'Rv0003': 'No'} if mycobrowser else {})
            assert set(master.loc[master['pseudogene'] == 'yes', 'simple_id']) == called(tasks, name)
            assert len(called(tasks, name)) == pseudogenes(k)

    def test_no_pseudogene_amap_called_is_in_the_clean_proteome(self, run):
        # all but Rv0003, whose protein has a stop, were in it before
        tasks, _ = run
        for name, _, _ in SAMPLES:
            assert set(ids(clean(tasks, name))) & called(tasks, name) == set()

    def test_every_other_protein_is_there_as_fix_cds_wrote_it(self, run):
        # Rv0001 among them, which Mycobrowser marks as a pseudogene of H37Rv
        tasks, _ = run
        for k, (name, _, _) in enumerate(SAMPLES):
            proteins = records(tasks['outputs'][name] / f'{name}-AMAP.faa')
            kept = [(h, s) for h, s in proteins if h.split(';')[0] not in called(tasks, name)]
            assert clean(tasks, name) == kept
            assert len(kept) == len(proteins) - pseudogenes(k)

    def test_the_genes_the_gff3_and_the_genbank_give_a_protein(self, run):
        # the CDS the GFF3 does not mark pseudo=true, and the GenBank translates
        tasks, _ = run
        for name, _, _ in SAMPLES:
            out = tasks['outputs'][name]
            features = [line.split('\t') for line in
                        (out / f'{name}-AMAP.gff3').read_text().split('##FASTA')[0].splitlines()
                        if line and not line.startswith('#')]
            gff3 = [dict(a.split('=', 1) for a in f[8].split(';'))['ID'] for f in features
                    if f[2] == 'CDS' and 'pseudo=true' not in f[8].split(';')]
            genbank = {f.qualifiers['locus_tag'][0]: f.qualifiers['translation'][0]
                       for record in SeqIO.parse(out / f'{name}-AMAP.gbk', 'genbank')
                       for f in record.features if f.type == 'CDS' and 'translation' in f.qualifiers}
            assert sorted(ids(clean(tasks, name))) == sorted(gff3)
            assert dict(zip(ids(clean(tasks, name)), (s for _, s in clean(tasks, name)))) == genbank
