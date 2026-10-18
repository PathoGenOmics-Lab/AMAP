#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A custom reference (--ref_annotation + --ref_genome) with a gene across the
origin of a circular sequence, and with more than one sequence.

PREPARE_REFERENCE extracts the reference's proteins, annotates the reference
with Bakta and lists what Bakta calls there that the annotation does not have:
those proteins go to miniprot as "bad Bakta", and a sample gene they cover is
removed. With E. coli O157:H7 Sakai (NCBI's GCF_000008865.2) as the reference:

- tagA on plasmid pO157 (92721 bp) is join(92527..92721,1..2502). From the
  GenBank its protein was named pO157p01_1-92721, the whole plasmid, and lost
  its stop (92721 bp is not 3 x 899 residues); from the GFF3, where NCBI
  writes it 92527..95223, it came out as its first 65 residues of 898.
- Bakta called stcE there as 92563..2502, start after end, which matched no
  gene at all, so its protein was listed as a false positive.
- The calls were compared with the genes of every sequence at once: a call on
  a plasmid counted as matched by a chromosome gene at the same numbers, and
  tagA's 1-92721 matched every call in the first 92721 bp of all three.

The reference here is a 2000 bp chromosome with gene A at 1501..1800 and the
3000 bp plasmid of tests/test_origin_genes.py, with gene W from 2851 across
the origin to 150 and gene G at 1001..1300.
"""

import os
import subprocess
import sys

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))

from classify import circular
from tests.test_origin_genes import N, PROTEIN, plasmid, revcomp, W
from tests.test_reference_files import _run, _write

BIN = os.path.join(os.path.dirname(__file__), '..', 'bin')
A = "ATG" + "GCT" * 98 + "TAA"                  # 300 bp at 1501..1800 of the chromosome
CHROMOSOME = "C" * 1500 + A + "C" * 200
G_PROTEIN = "M" + "K" * 98
FASTA = (f">chr chromosome, complete genome\n{CHROMOSOME}\n"
         f">p1 plasmid p1, complete sequence\n{plasmid()}\n")


@pytest.fixture(autouse=True)
def lengths():
    yield
    circular.set_lengths({})


def _genbank(path, plasmid_seq=None, w_location=None, plasmid_name='p1'):
    """The reference as NCBI writes it: W as join(2851..3000,1..150)."""
    w_location = w_location or CompoundLocation([SimpleLocation(2850, 3000, 1), SimpleLocation(0, 150, 1)])
    path.parent.mkdir(parents=True, exist_ok=True)
    records = []
    for name, seq, features in (
            ('chr', CHROMOSOME, [('gA', SimpleLocation(1500, 1800, 1))]),
            (plasmid_name, plasmid_seq or plasmid(), [('gG', SimpleLocation(1000, 1300, 1)), ('gW', w_location)])):
        record = SeqRecord(Seq(seq), id=name, name=name, description='',
                           annotations={'molecule_type': 'DNA', 'topology': 'circular'})
        record.features = [SeqFeature(loc, type='CDS', qualifiers={'locus_tag': [tag]}) for tag, loc in features]
        records.append(record)
    SeqIO.write(records, path, 'genbank')
    return path


def _gff3(path, w_lines, strand='+'):
    with open(path, 'w') as f:
        f.write("##gff-version 3\n"
                "chr\t.\tCDS\t1501\t1800\t.\t+\t0\tID=cds-A;locus_tag=gA\n"
                "p1\t.\tCDS\t1001\t1300\t.\t+\t0\tID=cds-G;locus_tag=gG\n")
        f.writelines(f"p1\t.\tCDS\t{s}\t{e}\t.\t{strand}\t0\tID=cds-W;locus_tag=gW\n" for s, e in w_lines)
    return path


def _extract(tmp_path, annotation, fasta=FASTA):
    genome = _write(tmp_path / 'ref.fasta', fasta)
    run = subprocess.run([sys.executable, os.path.join(BIN, 'extract_ref_proteins.py'),
                          '--ref-annotation', str(annotation), '--ref-genome', str(genome)],
                         cwd=tmp_path, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr[-2000:]
    proteins = {}
    for record in SeqIO.parse(tmp_path / 'custom_ref_ref_proteins.faa', 'fasta'):
        proteins[record.description] = str(record.seq)
    return proteins


WHOLE = {'gA_1501-1800 contig=chr': "M" + "A" * 98 + "*",
         'gG_1001-1300 contig=p1': G_PROTEIN + "*",
         'gW_2851-3150 contig=p1': PROTEIN + "*"}


# =============================================================================
# extract_ref_proteins.py
# =============================================================================

class TestExtract:

    def test_a_genbank_join_across_the_origin_ends_past_the_last_base(self, tmp_path):
        """It was gW_1-3000, the whole plasmid, and without its stop."""
        assert _extract(tmp_path, _genbank(tmp_path / 'ref.gbff')) == WHOLE

    def test_a_gff3_line_past_the_end_is_read_on_across_the_origin(self, tmp_path):
        """NCBI's way: the protein was cut at the plasmid's last base."""
        assert _extract(tmp_path, _gff3(tmp_path / 'ref.gff3', [(2851, 3150)])) == WHOLE

    def test_a_gff3_cds_cut_in_two_at_the_origin_is_joined_in_its_order(self, tmp_path):
        """Its lines were joined lowest first, 1..150 before 2851..3000."""
        assert _extract(tmp_path, _gff3(tmp_path / 'ref.gff3', [(1, 150), (2851, 3000)])) == WHOLE

    def test_the_same_on_the_minus_strand(self, tmp_path):
        rc = revcomp(W)
        minus = rc[150:] + "C" * (N - 300) + rc[:150]
        fasta = FASTA.replace(plasmid(), minus)
        # Biopython lists the parts of complement(join(2851..3000,1..150)) from the 5' end
        location = CompoundLocation([SimpleLocation(0, 150, -1), SimpleLocation(2850, 3000, -1)])
        want = {'gW_2851-3150 contig=p1': PROTEIN + "*"}
        got = _extract(tmp_path, _genbank(tmp_path / 'ref.gbff', minus, location), fasta)
        assert {k: v for k, v in got.items() if k.startswith('gW')} == want
        got = _extract(tmp_path, _gff3(tmp_path / 'ref.gff3', [(2851, 3150)], '-'), fasta)
        assert {k: v for k, v in got.items() if k.startswith('gW')} == want

    def test_a_join_that_does_not_cross_the_origin_keeps_its_span(self, tmp_path):
        """prfB's programmed frameshift on the minus strand, as Biopython
        lists its parts: from the 5' end."""
        location = CompoundLocation([SimpleLocation(1150, 1300, -1), SimpleLocation(1000, 1149, -1)])
        got = _extract(tmp_path, _genbank(tmp_path / 'ref.gbff', w_location=location))
        assert 'gW_1001-1300 contig=p1' in got

    def test_proteins_name_their_sequence_as_the_genome_does(self, tmp_path):
        """By the record's name when the FASTA has it, or else by the one
        sequence identical to it: Bakta's table names the FASTA's."""
        got = _extract(tmp_path, _genbank(tmp_path / 'ref.gbff', plasmid_name='pX'))
        assert got == WHOLE


# =============================================================================
# find_bakta_false_positives.py
# =============================================================================

CURATED = "".join(f">{h}\n{p}\n" for h, p in WHOLE.items())


def _false_positives(tmp_path, calls, curated=CURATED, fasta=FASTA):
    """The script on calls of (sequence, start, end, strand, locus tag);
    (exit code, false positives, output)."""
    tsv = _write(tmp_path / 'bakta_ref.tsv',
                 "# Annotated with Bakta\n#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
                 + "".join(f"{c}\tcds\t{s}\t{e}\t{st}\t{tag}\t\thypothetical protein\t\n"
                           for c, s, e, st, tag in calls))
    faa = _write(tmp_path / 'bakta_ref.faa', "".join(f">{c[4]} hypothetical protein\nMKKKKKKKKK\n" for c in calls))
    run = subprocess.run([sys.executable, os.path.join(BIN, 'find_bakta_false_positives.py'),
                          '--bakta-tsv', str(tsv), '--curated-faa', str(_write(tmp_path / 'curated.faa', curated)),
                          '--bakta-faa', str(faa), '--ref-genome', str(_write(tmp_path / 'ref.fasta', fasta))],
                         cwd=tmp_path, capture_output=True, text=True)
    bad = tmp_path / 'custom_ref_bad_bakta.faa'
    found = {r.id for r in SeqIO.parse(bad, 'fasta')} if bad.exists() else None
    return run.returncode, found, run.stdout + run.stderr


class TestFalsePositives:

    def test_a_call_across_the_origin_matches_its_gene(self, tmp_path):
        """stcE on pO157: Bakta starts it a few codons in, 2863..150 here."""
        code, bad, out = _false_positives(tmp_path, [('p1', 2863, 150, '+', 'b1')])
        assert code == 0, out
        assert bad == set()

    def test_a_piece_of_it_on_either_side_of_the_origin_matches_too(self, tmp_path):
        """Bakta's edge fragments of a gene it did not join (Rop on pOSAK1)."""
        code, bad, out = _false_positives(tmp_path, [('p1', 1, 120, '+', 'b1'), ('p1', 2900, 3000, '+', 'b2')])
        assert code == 0, out
        assert bad == set()

    def test_a_call_across_the_origin_where_no_gene_is_listed(self, tmp_path):
        code, bad, out = _false_positives(tmp_path, [('p1', 2863, 150, '-', 'b1'), ('chr', 1901, 99, '+', 'b2')],
                                          curated=">gG_1001-1300 contig=p1\nM\n>gA_1501-1800 contig=chr\nM\n")
        assert code == 0, out
        assert bad == {'b1', 'b2'}

    def test_a_call_is_matched_by_the_genes_of_its_own_sequence_only(self, tmp_path):
        """A at 1501..1800 is on the chromosome and G at 1001..1300 on the
        plasmid; calls at those numbers on the other sequence match nothing."""
        calls = [('p1', 1501, 1800, '+', 'b1'), ('chr', 1001, 1300, '+', 'b2'),
                 ('chr', 1501, 1800, '+', 'b3'), ('p1', 1001, 1300, '+', 'b4')]
        code, bad, out = _false_positives(tmp_path, calls)
        assert code == 0, out
        assert bad == {'b1', 'b2'}

    def test_calls_on_a_sequence_the_annotation_has_no_gene_on_are_not_judged(self, tmp_path):
        fasta = FASTA + f">p2 plasmid p2, complete sequence\n{'C' * 500}\n"
        code, bad, out = _false_positives(tmp_path, [('p2', 101, 400, '+', 'b1'), ('chr', 101, 400, '+', 'b2')],
                                          fasta=fasta)
        assert code == 0, out
        assert bad == {'b2'}
        assert 'no gene on p2' in out

    def test_sequences_named_differently_stop_the_run(self, tmp_path):
        """Bakta without --keep-contig-headers names them contig_1, contig_2."""
        code, bad, out = _false_positives(tmp_path, [('contig_1', 101, 400, '+', 'b1')])
        assert code != 0
        assert '--keep-contig-headers' in out

    def test_genes_whose_header_names_no_sequence_match_everywhere(self, tmp_path):
        """A proteome written before the sequences were named, as before."""
        code, bad, out = _false_positives(tmp_path, [('p1', 1501, 1800, '+', 'b1'), ('p1', 101, 400, '+', 'b2')],
                                          curated=">gA_1501-1800\nM\n")
        assert code == 0, out
        assert bad == {'b2'}


# =============================================================================
# PREPARE_REFERENCE, its script as Nextflow renders it
# =============================================================================

# Bakta, as far as the step needs it: it names each sequence contig_N unless
# told to keep the FASTA's names, writes a call across the origin start after
# end, and writes its proteins by locus tag.
FAKE_BAKTA = '''
    import sys
    from pathlib import Path
    args = sys.argv[1:]
    def option(name):
        return args[args.index(name) + 1]
    CALLS = {'chr': [(1501, 1800, '+'), (101, 400, '+')],
             'p1': [(2863, 150, '+'), (1001, 1300, '+'), (1501, 1800, '-')]}
    names = [line[1:].split()[0] for line in open(args[-1]) if line.startswith('>')]
    out = Path(option('--output'))
    out.mkdir(exist_ok=True)
    rows, proteins, n = [], [], 0
    for i, name in enumerate(names, 1):
        for start, end, strand in CALLS.get(name, []):
            n += 1
            tag = f"{option('--locus-tag')}_{n:05d}"
            seq_id = name if '--keep-contig-headers' in args else f'contig_{i}'
            rows.append(f"{seq_id}\\tcds\\t{start}\\t{end}\\t{strand}\\t{tag}\\t\\thypothetical protein\\t\\n")
            proteins.append(f">{tag} hypothetical protein\\nMKKKKKKKKK\\n")
    prefix = option('--prefix')
    (out / f'{prefix}.tsv').write_text('# Annotated with Bakta\\n#Sequence Id\\tType\\tStart\\tStop\\tStrand'
                                       '\\tLocus Tag\\tGene\\tProduct\\tDbXrefs\\n' + ''.join(rows))
    (out / f'{prefix}.faa').write_text(''.join(proteins))
'''

PYTHON = '''
    import os, sys
    os.execv(sys.executable, [sys.executable] + sys.argv[1:])
'''


class TestPrepareReference:

    def test_only_what_bakta_invents_is_listed(self, tmp_path):
        """W across the origin, G and A are the annotation's; the call on the
        chromosome at 101..400 and the one on the plasmid at A's numbers
        are not. W was listed and the plasmid's invention was not."""
        inputs = {'ref_annotation': _genbank(tmp_path / 'store' / 'ref.gbff'),
                  'ref_genome': _write(tmp_path / 'store' / 'ref.fasta', FASTA),
                  'bakta_db': []}
        params = {'bakta_db': '/data/bakta_db/db-light', 'genus': 'Escherichia', 'species': 'Escherichia coli'}
        code, err, work = _run(tmp_path, 'modules/local/reference.nf', 'PREPARE_REFERENCE', inputs, params,
                               {'bakta': FAKE_BAKTA, 'python3': PYTHON})
        assert code == 0, err
        bad = {r.id for r in SeqIO.parse(work / 'custom_ref_bad_bakta.faa', 'fasta')}
        assert bad == {'bakta_ref_00002', 'bakta_ref_00005'}
        proteins = {r.description for r in SeqIO.parse(work / 'custom_ref_ref_proteins.faa', 'fasta')}
        assert proteins == set(WHOLE)
