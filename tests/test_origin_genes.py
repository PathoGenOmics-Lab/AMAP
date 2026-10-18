#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Genes that run across the origin of a circular sequence.

Bakta joins the ends of a sequence it takes for a complete replicon and calls
a gene across its origin from its start, near the end, to its stop, near the
beginning. On E. coli O157:H7 Sakai with NCBI's headers it called stcE on
plasmid pO157 as 92563..2502, and TYPE_CASES stopped on it ("Null Interval
objects not allowed in IntervalTree"); the sample was dropped and the run
ended with exit 0. Everything after it read genes as start..end too: check_cds
dropped such a gene as invalid coordinates, deduplicate_rvs measured it as a
negative length, the GenBank writer stopped on it, and the GFF3 and Panaroo
files would have written it backwards or over the whole sequence.

The genome here is a 3000 bp plasmid described as complete, with gene W
(ATG, 98 codons that repeat nothing, TAA) from 2851 across the origin to 150.
"""

import os
import subprocess
import sys
import tempfile

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))

from classify import circular

BIN = os.path.join(os.path.dirname(__file__), '..', 'bin')
N = 3000
# One codon for each amino acid but M, in an order that does not repeat: a
# protein that lines up with its reference from any residue would hide a
# missing N-terminus.
CODONS = dict(zip("AKDWRGHFIEPSTYVLNQC", ["GCT", "AAA", "GAT", "TGG", "CGT", "GGC", "CAT", "TTC", "ATC",
                                        "GAG", "CCA", "AGC", "ACC", "TAT", "GTG", "CTG", "AAC", "CAG",
                                        "TGC"]))
BODY = "".join("AKDWRGHFIEPSTYVLNQC"[(1103515245 * i * i + 12345 * i + 7) % 2 ** 31 % 19] for i in range(98))
W = "ATG" + "".join(CODONS[a] for a in BODY) + "TAA"     # 300 bp
PROTEIN = "M" + BODY
G = "ATG" + "AAA" * 98 + "TAG"                 # 300 bp at 1001..1300
FILL = "C"


def plasmid():
    """W from 2851 across the origin to 150, G at 1001..1300, C elsewhere."""
    seq = W[150:] + FILL * 850 + G + FILL * (N - 150 - 850 - 300 - 150) + W[:150]
    assert len(seq) == N and seq[2850:] + seq[:150] == W and seq[1000:1300] == G
    return seq


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


@pytest.fixture(autouse=True)
def lengths():
    """Sequence lengths, as classify_cases.py keeps them from the genome."""
    circular.set_lengths({'c1': N})
    yield
    circular.set_lengths({})


# =============================================================================
# classify/circular.py
# =============================================================================

class TestCircular:

    def test_a_gene_across_the_origin_has_its_start_after_its_end(self):
        assert circular.wraps(2851, 150)
        assert not circular.wraps(150, 2851)
        # miniprot's unplaced proteins read as 1..0: no place at all
        assert not circular.wraps(1, 0)

    def test_pieces_copies_and_length(self):
        assert circular.pieces(2851, 150, 'c1') == [(2851, 3000), (1, 150)]
        assert circular.copies(2851, 150, 'c1') == [(2851, 3150), (-149, 150)]
        assert circular.length(2851, 150, 'c1') == 300
        assert circular.length(100, 399, 'c1') == 300
        assert circular.copies(100, 399, 'c1') == [(100, 399)]

    def test_shared_bases_on_either_side_of_the_origin(self):
        assert circular.shared((2851, 150), (1, 100), 'c1') == 100
        assert circular.shared((2851, 150), (2901, 3000), 'c1') == 100
        assert circular.shared((2851, 150), (2801, 50), 'c1') == 150 + 50
        assert circular.shared((2851, 150), (200, 2800), 'c1') == 0
        # a linear feature reaching both pieces
        assert circular.shared((2851, 150), (100, 2900), 'c1') == 51 + 50

    def test_frames_compare_on_the_copy_over_the_other(self):
        """3000 is a multiple of three, so try a length that is not."""
        circular.set_lengths({'c1': 3001})
        assert circular.facing((2851, 150), (1, 90), 'c1') == (2851 - 3001, 150)
        assert circular.facing((2851, 150), (2900, 2990), 'c1') == (2851, 150 + 3001)

    def test_an_unknown_sequence_length_is_an_error_not_a_guess(self):
        with pytest.raises(ValueError, match='c9'):
            circular.length(2851, 150, 'c9')

    def test_what_bakta_takes_for_a_circular_sequence(self):
        """Bakta 1.12 reads the words after the first space."""
        assert circular.bakta_circular('NC_002128.1 Escherichia coli O157:H7 str. Sakai plasmid pO157, '
                                       'complete sequence')
        assert circular.bakta_circular('1 length=4411532 depth=1.00x circular=true')
        assert not circular.bakta_circular('contig_1 length=4411532 circular=false')
        assert not circular.bakta_circular('c1 [completeness=complete=false]')
        assert not circular.bakta_circular('L4_H37Rv')
        assert not circular.bakta_circular('complete_genome')
        # Bakta's test is a substring: "incomplete" reads as complete too, and
        # what matters is what Bakta did
        assert circular.bakta_circular('c1 incomplete')

    def test_read_lengths(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'g.fasta')
            with open(path, 'w') as f:
                f.write(">c1 plasmid, complete sequence\nACGT\nAC\n>c2\nA\n")
            assert circular.read_lengths(path) == {'c1': 6, 'c2': 1}


# =============================================================================
# The classifier
# =============================================================================

class TestClassifier:

    @staticmethod
    def _bakta(rows):
        df = pd.DataFrame(rows, columns=['Locus Tag', 'Start', 'End', 'Strand', 'Contig'])
        df['Length_nt'] = [circular.length(s, e, c) for s, e, c in zip(df['Start'], df['End'], df['Contig'])]
        return df[['Locus Tag', 'Start', 'End', 'Strand', 'Length_nt', 'Contig']]

    @staticmethod
    def _mp(rows):
        return pd.DataFrame(rows, columns=['Prot Id', 'Prot len', 'Strand', 'Start', 'End',
                                           'Match Nt', 'Num nt', 'Contig'])

    def test_a_call_across_the_origin_is_measured_whole(self):
        from classify.io import read_bakta
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'S1.tsv')
            with open(path, 'w') as f:
                f.write("#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
                        "c1\tcds\t2851\t150\t+\tS1_00001\tstcE\tmetalloprotease StcE\t\n")
            assert list(read_bakta(path)['Length_nt']) == [300]

    def test_an_alignment_past_the_origin_overlaps_the_call(self):
        """The call used to overlap nothing: case 1, and the alignment went in
        as a missing gene beside it."""
        from classify.cases import classify_case1_case4, find_rv_missing
        bk = self._bakta([['S1_00001', 2851, 150, '+', 'c1']])
        mp = self._mp([['Rv0001', 50, '+', 1, 150, 150, 150, 'c1']])
        case4, case1 = classify_case1_case4(bk, mp)
        assert list(case4['Locus Tag']) == ['S1_00001'] and case1.empty
        assert find_rv_missing(mp, bk).empty

    def test_a_call_across_the_origin_that_bad_calls_cover_is_removed(self):
        from classify.cases import filter_case1
        case1 = self._bakta([['S1_00001', 2851, 150, '+', 'c1']])
        bad = pd.DataFrame([['bad1', 100, '+', 2851, 3000, 150, 150, 'c1'],
                            ['bad2', 50, '+', 1, 120, 120, 120, 'c1']],
                           columns=['Prot Id', 'Prot len', 'Strand', 'Start', 'End', 'Match Nt',
                                    'Num nt', 'Contig'])
        with tempfile.NamedTemporaryFile('w', suffix='.paf', delete=False) as f:
            pass
        try:
            keep, remove = filter_case1(case1, bad, f.name)
        finally:
            os.unlink(f.name)
        assert list(remove['Locus Tag']) == ['S1_00001'] and keep.empty
        assert remove.iloc[0]['Overlap_Percentage'] == 90.0

    def test_a_much_longer_call_across_the_origin_wins_over_a_fragment(self):
        """Its length came out negative, so any alignment beat it."""
        from classify.overlaps import resolve_single_overlaps
        mp = self._mp([['Rv0001', 300, '+', 1, 150, 150, 150, 'c1']]).assign(Complete=False)
        circular.set_lengths({'c1': 5000})
        groups = {0: {'case4': [(3301, 150)], 'case4_strands': ['+'], 'case4_contigs': ['c1'],
                      'miniprot_df': [(1, 150)], 'contig': 'c1', 'Prot Id': 'Rv0001',
                      'Locus Tags': ['S1_00001']}}
        combined = {}
        resolve_single_overlaps(groups, mp, combined, {})
        assert combined == {(3301, 150, '+', 'c1'): 'S1_00001'}

    def test_an_alignment_over_most_of_a_call_across_the_origin_takes_its_place(self):
        """As over any call of about its length. They were compared as if
        they did not overlap, and both stayed."""
        from classify.overlaps import _insert_with_dedup
        combined = {(2981, 280, '+', 'c1'): 'S1_00001'}
        _insert_with_dedup(combined, (1, 280, '+', 'c1'), 'Rv0001', reference_ids={'Rv0001'})
        assert combined == {(1, 280, '+', 'c1'): 'Rv0001'}

    def test_recovery_reads_a_gene_across_the_origin_on_both_sides(self):
        """It stopped the run: the trees refuse an interval that ends before
        it starts. Now it accounts for a placement in its frame past the
        origin, whose name it takes as it has none of its own, and not for
        one in another frame, which goes in beside it."""
        from classify.reconciliation import recover_missing_genes
        # 3001 bases: the gene is read in frame 1 (2851 % 3) before the origin
        # and in frame 0 ((2851 - 3001) % 3) after it
        circular.set_lengths({'c1': 3001})
        mp = pd.DataFrame([['Rv0001', 33, 149, '+', 117, 117, 'c1'],
                           ['Rv0002', 34, 141, '+', 108, 108, 'c1']],
                          columns=['Prot Id', 'Start', 'End', 'Strand', 'Match Nt', 'Num nt', 'Contig'])
        combined = {(2851, 150, '+', 'c1'): 'S1_00001'}
        recover_missing_genes(combined, mp)
        assert combined == {(2851, 150, '+', 'c1'): 'Rv0001', (34, 141, '+', 'c1'): 'Rv0002'}

    def test_a_copy_inside_another_gene_across_the_origin_goes(self):
        """As anywhere else: a placement of a gene placed twice that lies
        inside another gene is dropped."""
        from classify.reconciliation import deduplicate
        combined = {(2851, 150, '+', 'c1'): 'Rv0001', (1001, 1300, '+', 'c1'): 'Rv0001',
                    (2801, 200, '+', 'c1'): 'S1_00009'}
        deduplicate(combined)
        assert combined == {(1001, 1300, '+', 'c1'): 'Rv0001', (2801, 200, '+', 'c1'): 'S1_00009'}


# =============================================================================
# classify_cases.py and check_cds.py on the plasmid
# =============================================================================

def _paf(prot, plen, strand, contig, start0, end, matched):
    return (f"{prot}\t{plen}\t0\t{plen}\t{strand}\t{contig}\t{N}\t{start0}\t{end}\t{matched}\t{matched}\t0"
            f"\tAS:i:500\tms:i:500\tnp:i:{plen}\tst:i:1\tcg:Z:{plen}M\n")


def classify(d, genome, bakta_rows, paf_rows, proteome, description='plasmid, complete sequence'):
    """classify_cases.py as TYPE_CASES runs it; the corrected table."""
    with open(os.path.join(d, 'g.fasta'), 'w') as f:
        f.write(f">c1 {description}\n{genome}\n")
    with open(os.path.join(d, 'S1.tsv'), 'w') as f:
        f.write("# Annotated with Bakta\n#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n")
        f.writelines(f"c1\tcds\t{s}\t{e}\t{st}\t{tag}\t\tproduct\t\n" for tag, s, e, st in bakta_rows)
    with open(os.path.join(d, 'S1.paf'), 'w') as f:
        f.writelines(paf_rows)
    for name in ('S1-badbakta.paf', 'S1-ref.paf'):
        open(os.path.join(d, name), 'w').close()
    with open(os.path.join(d, 'ref.faa'), 'w') as f:
        f.writelines(f">{h}\n{seq}\n" for h, seq in proteome.items())
    run = subprocess.run([sys.executable, os.path.join(BIN, 'classify_cases.py'), '-b', 'S1.tsv',
                          '-p', 'S1.paf', '-bb', 'S1-badbakta.paf', '-c', 'S1-ref.paf',
                          '--genome', 'g.fasta', '--ref-proteome', 'ref.faa'],
                         cwd=d, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr[-2000:]
    return pd.read_csv(os.path.join(d, 'S1_corrected.tsv'), sep='\t', header=None,
                       names=['id', 'strand', 'start', 'end', 'rvs', 'contig'], dtype={'rvs': str})


def check(d, rows, proteome=None, description='plasmid, complete sequence', genome=None, dedup=False):
    """check_cds.py on rows of (id, strand, start, end, rvs), and with dedup
    deduplicate_rvs.py after it, as FIX_CDS runs them; the table and the
    proteins."""
    with open(os.path.join(d, 'g.fasta'), 'w') as f:
        f.write(f">c1 {description}\n{genome or plasmid()}\n")
    with open(os.path.join(d, 'S1_corrected.tsv'), 'w') as f:
        f.writelines("\t".join(map(str, r)) + "\tc1\n" for r in rows)
    with open(os.path.join(d, 'ref.faa'), 'w') as f:
        f.writelines(f">{h}\n{seq}\n" for h, seq in (proteome or {'gW_1-300': PROTEIN + '*'}).items())
    steps = [[os.path.join(BIN, 'check_cds.py'), 'g.fasta', 'S1_corrected.tsv', '--ref-proteome', 'ref.faa']]
    if dedup:
        steps.append([os.path.join(BIN, 'deduplicate_rvs.py'), 'S1.AMAP.tsv', '--faa', 'S1.AMAP.faa',
                      '--fna', 'S1.AMAP.fna', '--ref-proteome', 'ref.faa', '--genome', 'g.fasta'])
    for cmd in steps:
        run = subprocess.run([sys.executable] + cmd, cwd=d, capture_output=True, text=True)
        assert run.returncode == 0, run.stderr[-2000:]
    table = pd.read_csv(os.path.join(d, 'S1.AMAP.tsv'), sep='\t', dtype=str, keep_default_na=False)
    proteins = {}
    with open(os.path.join(d, 'S1.AMAP.faa')) as fh:
        for line in fh:
            if line.startswith('>'):
                key = line[1:].strip()
                proteins[key] = ''
            else:
                proteins[key] += line.strip()
    return table.set_index('id'), proteins


class TestPipelineSteps:

    PROTEOME = {'gG_1001-1300': 'M' + 'K' * 98 + '*'}

    def test_type_cases_keeps_a_call_across_the_origin(self):
        """stcE on pO157 of E. coli Sakai: the call nothing in the reference
        matches is kept as Bakta made it, start after end."""
        with tempfile.TemporaryDirectory() as d:
            out = classify(d, plasmid(), [('S1_00001', 2851, 150, '+'), ('S1_00002', 1001, 1300, '+')],
                           [_paf('gG_1001-1300', 100, '+', 'c1', 1000, 1300, 300)], self.PROTEOME)
        got = {r.id: (r.start, r.end, r.strand) for r in out.itertuples()}
        assert got == {'S1_00001': (2851, 150, '+'), 'S1_00002': (1001, 1300, '+')}

    def test_fix_cds_translates_it_whole(self):
        """It was dropped as invalid coordinates."""
        with tempfile.TemporaryDirectory() as d:
            table, proteins = check(d, [('S1_00001', '+', 2851, 150, '')])
        assert (table.at['S1_00001', 'start'], table.at['S1_00001', 'end']) == ('2851', '150')
        assert (table.at['S1_00001', 'pseudogene'], table.at['S1_00001', 'partial']) == ('no', '')
        assert proteins == {'S1_00001;+;2850-150': PROTEIN}

    def test_the_same_on_the_minus_strand(self):
        genome = revcomp(plasmid())
        with tempfile.TemporaryDirectory() as d:
            table, proteins = check(d, [('S1_00001', '-', N - 149, 150, '')], genome=genome)
        assert (table.at['S1_00001', 'start'], table.at['S1_00001', 'end']) == (str(N - 149), '150')
        assert list(proteins.values()) == [PROTEIN]

    def test_a_gene_that_reaches_the_origin_runs_on_to_its_stop(self):
        """Bakta called the first 141 bp of Rop on pOSAK1 (3306 bp); its stop
        lies past the origin. Read across it, the gene is whole."""
        with tempfile.TemporaryDirectory() as d:
            table, proteins = check(d, [('S1_00001', '+', 2851, 3000, '')])
        assert (table.at['S1_00001', 'start'], table.at['S1_00001', 'end']) == ('2851', '150')
        assert table.at['S1_00001', 'partial'] == ''
        assert list(proteins.values()) == [PROTEIN]

    def test_on_a_linear_sequence_the_end_still_cuts_it(self):
        """No "complete" in the description: Bakta read the sequence as linear,
        and so does AMAP."""
        with tempfile.TemporaryDirectory() as d:
            table, _ = check(d, [('S1_00001', '+', 2851, 3000, '')], description='plasmid')
        assert (table.at['S1_00001', 'start'], table.at['S1_00001', 'end']) == ('2851', '3000')
        assert table.at['S1_00001', 'partial'] == "3'"

    def test_a_missing_n_terminus_is_found_past_the_origin(self):
        """miniprot cannot align across the origin; the piece it placed after
        it is completed from the reference protein."""
        with tempfile.TemporaryDirectory() as d:
            table, proteins = check(d, [('S1_gW', '+', 1, 150, 'gW')])
        assert (table.at['S1_gW', 'start'], table.at['S1_gW', 'end']) == ('2851', '150')
        assert list(proteins.values()) == [PROTEIN]

    def test_genes_away_from_the_origin_read_as_before(self):
        rows = [('S1_00002', '+', 1001, 1290, '')]          # no stop: extended to 1300
        with tempfile.TemporaryDirectory() as d:
            circ, _ = check(d, rows)
        with tempfile.TemporaryDirectory() as d:
            lin, _ = check(d, rows, description='plasmid')
        assert circ.to_dict() == lin.to_dict()
        assert (circ.at['S1_00002', 'start'], circ.at['S1_00002', 'end']) == ('1001', '1300')

    def test_fix_cds_keeps_the_whole_gene_and_drops_its_piece(self):
        """miniprot's piece past the origin, completed, is the gene Bakta
        called; FIX_CDS keeps one of them."""
        rows = [('S1_00001', '+', 2851, 150, ''), ('S1_gW', '+', 1, 150, 'gW')]
        with tempfile.TemporaryDirectory() as d:
            table, proteins = check(d, rows, dedup=True)
        assert len(table) == 1 and (table.iloc[0]['start'], table.iloc[0]['end']) == ('2851', '150')
        assert list(proteins.values()) == [PROTEIN]

    def test_dedup_keeps_the_whole_gene_over_a_piece_of_it(self):
        """Measured as a negative length, the whole gene lost to a piece."""
        from deduplicate_rvs import deduplicate
        genes = [{'id': 'S1_gW', 'strand': '+', 'start': '1', 'end': '150', 'rvs': 'gW',
                  'pseudogene': 'no', 'contig': 'c1'},
                 {'id': 'S1_gW', 'strand': '+', 'start': '2851', 'end': '150', 'rvs': 'gW',
                  'pseudogene': 'no', 'contig': 'c1'}]
        kept, _ = deduplicate(genes)
        assert [(g['start'], g['end']) for g in kept] == [('2851', '150')]

    def test_dedup_overlap_counts_across_the_origin(self):
        from deduplicate_rvs import deduplicate
        genes = [{'id': 'S1_gW', 'strand': '+', 'start': '2851', 'end': '150', 'rvs': 'gW',
                  'pseudogene': 'no', 'contig': 'c1'},
                 {'id': 'S1_00001', 'strand': '+', 'start': '2861', 'end': '150', 'rvs': '',
                  'pseudogene': 'no', 'contig': 'c1'}]
        kept, _ = deduplicate(genes)
        assert [g['id'] for g in kept] == ['S1_gW']


# =============================================================================
# GET_OUTPUTS: the GFF3, the GenBank, Panaroo's GFF3 and the QC table
# =============================================================================

class TestOutputs:

    def _write(self, d, description='plasmid, complete sequence'):
        from generate_tsv_gff_gbk import generate_tsv, generate_gff, tsv_to_gbk
        merged = pd.DataFrame({
            'id': ['S1_00001', 'S1_00002'], 'strand': ['+', '+'], 'start': [2851, 1001],
            'end': [150, 1300], 'rvs': ['', 'gG'], 'Preferred_name': ['-', '-'], 'Gene': ['w', 'g'],
            'Product_x': ['p', 'q'], 'Description': ['-', '-'], 'contig': ['c1', 'c1'],
            'pseudogene': ['no', 'no'], 'partial': ['', ''], 'note': ['', '']})
        rest = pd.DataFrame(columns=["Sequence Id", "Type", "Start", "Stop", "Strand",
                                     "Locus Tag", "Gene", "Product", "DbXrefs"])
        fasta = os.path.join(d, 'g.fasta')
        with open(fasta, 'w') as f:
            f.write(f">c1 {description}\n{plasmid()}\n")
        tsv = os.path.join(d, 'S1-AMAP.tsv')
        out = generate_tsv(merged, rest, fasta, tsv)
        generate_gff(None, out, os.path.join(d, 'S1-AMAP.gff3'), fasta)
        tsv_to_gbk(tsv, fasta, os.path.join(d, 'S1-AMAP.gbk'))
        return out

    def test_gff3_ends_past_the_last_base_on_a_circular_region(self):
        """GFF3 writes a feature across the origin of a circular landmark to
        its end plus the landmark's length, as Bakta's does."""
        with tempfile.TemporaryDirectory() as d:
            self._write(d)
            with open(os.path.join(d, 'S1-AMAP.gff3')) as f:
                lines = [l.rstrip('\n').split('\t') for l in f if '\t' in l and not l.startswith('>')]
        region = [l for l in lines if l[2] == 'region']
        cds = {l[8].split(';')[0][3:]: (int(l[3]), int(l[4])) for l in lines if l[2] == 'CDS'}
        assert region == [['c1', 'AMAP', 'region', '1', str(N), '.', '+', '.', 'ID=c1;Name=c1;Is_circular=true']]
        assert cds == {'S1_00001': (2851, N + 150), 'S1_00002': (1001, 1300)}

    def test_a_linear_sequence_gets_no_region(self):
        with tempfile.TemporaryDirectory() as d:
            out = self._write(d, description='plasmid')
            out = out[out['Locus Tag'] == 'S1_00002']
            from generate_tsv_gff_gbk import generate_gff
            generate_gff(None, out, os.path.join(d, 'only.gff3'), os.path.join(d, 'g.fasta'))
            with open(os.path.join(d, 'only.gff3')) as f:
                assert '\tregion\t' not in f.read()

    def test_genbank_joins_it_across_the_origin(self):
        from Bio import SeqIO
        with tempfile.TemporaryDirectory() as d:
            self._write(d)
            record = next(SeqIO.parse(os.path.join(d, 'S1-AMAP.gbk'), 'genbank'))
            with open(os.path.join(d, 'S1-AMAP.gbk')) as f:
                text = f.read()
        assert record.annotations['topology'] == 'circular'
        w = [f for f in record.features if f.qualifiers['locus_tag'] == ['S1_00001']][0]
        assert 'join(2851..3000,1..150)' in text
        assert w.qualifiers['translation'] == [PROTEIN]

    def test_panaroo_gff3_does_not_span_the_whole_sequence(self):
        """Biopython spans a join across the origin as the whole sequence,
        1..3000."""
        with tempfile.TemporaryDirectory() as d:
            self._write(d)
            subprocess.run([sys.executable, os.path.join(BIN, 'to_panaroo.py'), '--id', 'S1'],
                           cwd=d, check=True, capture_output=True)
            with open(os.path.join(d, 'S1.panaroo.gff')) as f:
                cds = [l.split('\t') for l in f if '\tCDS\t' in l]
        assert {(c[3], c[4]) for c in cds} == {('2851', str(N + 150)), ('1001', '1300')}

    def test_qc_measures_it_whole(self):
        from annotation_qc import compute_qc
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'g.fasta'), 'w') as f:
                f.write(f">c1 plasmid, complete sequence\n{plasmid()}\n")
            pd.DataFrame({'id': ['S1_00001', 'S1_00002'], 'start': [2851, 1001], 'end': [150, 1300],
                          'contig': ['c1', 'c1']}).to_csv(os.path.join(d, 'm.tsv'), sep='\t', index=False)
            compute_qc(os.path.join(d, 'm.tsv'), os.path.join(d, 'g.fasta'), os.path.join(d, 'qc.tsv'))
            qc = pd.read_csv(os.path.join(d, 'qc.tsv'), sep='\t').iloc[0]
        assert qc['mean_gene_length'] == 300 and qc['coding_density_pct'] == 20.0
