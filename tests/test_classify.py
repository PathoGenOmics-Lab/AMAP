#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the AMAP classification module.

Tests cover:
- ID cleaning and formatting (ids.py)
- Overlap calculations (overlaps.py)
- Recovery and deduplication (reconciliation.py)
- Case classification (cases.py)
- I/O functions (io.py)
"""

import os
import sys
import tempfile
import pytest
import pandas as pd
import numpy as np

# Add bin/ to path so we can import classify
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))

from classify.ids import clean_id, format_cluster_id
from classify.overlaps import (
    get_overlap_percentage, get_subset_overlap_percentage,
    _range_length, _insert_with_dedup, _is_reference_gene
)


# =============================================================================
# IDS MODULE
# =============================================================================

class TestCleanId:
    """Tests for clean_id — strips embedded coordinates from gene IDs."""
    
    def test_simple_rv_with_coords(self):
        assert clean_id('Rv0609A_703830-704057') == 'Rv0609A'
    
    def test_prefixed_rv_with_coords(self):
        assert clean_id('L8_G26689_Rv0609A_703830-704057') == 'L8_G26689_Rv0609A'
    
    def test_pe_pgrs_with_coords(self):
        """PE_PGRS genes have 5+ underscore parts — old 4-part rule failed here."""
        assert clean_id('PE_PGRS1_703830-704057') == 'PE_PGRS1'
    
    def test_prefixed_pe_pgrs_with_coords(self):
        assert clean_id('L8_G26689_PE_PGRS1_703830-704057') == 'L8_G26689_PE_PGRS1'
    
    def test_no_coords_unchanged(self):
        assert clean_id('L8_G26689_00042') == 'L8_G26689_00042'
    
    def test_simple_id_unchanged(self):
        assert clean_id('Rv0001') == 'Rv0001'
    
    def test_empty_string(self):
        assert clean_id('') == ''
    
    def test_only_coords(self):
        """Edge case: just coordinates, no gene name."""
        assert clean_id('_100-200') == ''
    
    def test_large_coords(self):
        assert clean_id('dnaA_1-1524') == 'dnaA'
    
    def test_non_mtbc_gene(self):
        """Organism-agnostic: non-Rv gene IDs should also work."""
        assert clean_id('thrA_337-2799') == 'thrA'
        assert clean_id('ECOLI_thrA_337-2799') == 'ECOLI_thrA'


class TestFormatClusterId:
    """Tests for format_cluster_id — adds sample prefix to gene IDs."""
    
    def test_simple_gene(self):
        assert format_cluster_id('Rv0001', 'L8_G26689') == 'L8_G26689_Rv0001'
    
    def test_gene_with_coords(self):
        assert format_cluster_id('Rv0001_1-1524', 'L8_G26689') == 'L8_G26689_Rv0001'
    
    def test_already_prefixed(self):
        assert format_cluster_id('L8_G26689_Rv0001', 'L8_G26689') == 'L8_G26689_Rv0001'
    
    def test_different_prefix(self):
        result = format_cluster_id('Rv0001', 'L4_G12345')
        assert result == 'L4_G12345_Rv0001'


# =============================================================================
# OVERLAPS MODULE
# =============================================================================

class TestRangeLength:
    """Tests for _range_length — 1-based inclusive length."""
    
    def test_single_base(self):
        assert _range_length((100, 100)) == 1
    
    def test_normal_range(self):
        assert _range_length((100, 200)) == 101
    
    def test_gene_length(self):
        assert _range_length((1, 1524)) == 1524


class TestOverlapPercentage:
    """Tests for get_overlap_percentage."""
    
    def test_identical_ranges(self):
        assert get_overlap_percentage((100, 200), (100, 200)) == 100.0
    
    def test_no_overlap(self):
        assert get_overlap_percentage((100, 200), (300, 400)) == 0.0
    
    def test_partial_overlap(self):
        # Overlap region: 150-200 = 51 bases, max length = 101
        pct = get_overlap_percentage((100, 200), (150, 250))
        assert 40 < pct < 60  # ~50.5%
    
    def test_contained(self):
        # (150, 170) fully inside (100, 200)
        pct = get_overlap_percentage((150, 170), (100, 200))
        assert pct > 0  # overlap exists
    
    def test_adjacent_no_overlap(self):
        # Adjacent ranges don't overlap (200 end, 201 start)
        assert get_overlap_percentage((100, 200), (201, 300)) == 0.0


class TestSubsetOverlap:
    """Tests for get_subset_overlap_percentage."""
    
    def test_fully_contained(self):
        assert get_subset_overlap_percentage((150, 170), (100, 200)) == 100.0
    
    def test_no_overlap(self):
        assert get_subset_overlap_percentage((100, 200), (300, 400)) == 0.0
    
    def test_partial(self):
        pct = get_subset_overlap_percentage((100, 200), (150, 250))
        assert 0 < pct < 100


class TestIsReferenceGene:
    """A gene in the dict is a reference gene when its ID is a reference protein's."""

    def test_reference_protein(self):
        assert _is_reference_gene('Rv0001', {'Rv0001', 'Rv0002'}) is True

    def test_locus_tag(self):
        assert _is_reference_gene('ABCDE_00042', {'Rv0001'}) is False

    def test_pe_pgrs(self):
        assert _is_reference_gene('PE_PGRS1', {'PE_PGRS1'}) is True

    def test_nothing_is_one_without_the_reference_ids(self):
        """Coordinates in an ID are no sign: clean_id has removed them by then."""
        assert _is_reference_gene('Rv0001_1-1524') is False


class TestInsertWithDedup:
    """Tests for _insert_with_dedup — strand-aware overlap resolution."""
    
    def test_no_conflict(self):
        d = {}
        _insert_with_dedup(d, (100, 200, '+'), 'gene1')
        assert d == {(100, 200, '+'): 'gene1'}
    
    def test_opposite_strands_no_conflict(self):
        """Genes on opposite strands should NOT conflict."""
        d = {(100, 200, '+'): 'gene1'}
        _insert_with_dedup(d, (100, 200, '-'), 'gene2')
        assert len(d) == 2
        assert d[(100, 200, '+')] == 'gene1'
        assert d[(100, 200, '-')] == 'gene2'
    
    def test_same_strand_overlap_larger_wins(self):
        """When two genes on same strand overlap >70%, larger wins."""
        d = {(100, 200, '+'): 'small_gene'}
        _insert_with_dedup(d, (50, 800, '+'), 'big_gene', threshold=0)
        # big_gene should replace small_gene
        assert (50, 800, '+') in d
    
    def test_no_overlap_different_regions(self):
        d = {(100, 200, '+'): 'gene1'}
        _insert_with_dedup(d, (500, 600, '+'), 'gene2')
        assert len(d) == 2

    def test_rejected_insert_does_not_delete_existing(self):
        """A discarded candidate must not take other genes down with it.

        The candidate overlaps two entries: the first is queued for removal,
        the second (a reference gene) wins and aborts the insert. The queued
        removal was only valid as a replacement, so nothing may be deleted.
        """
        d = {
            (100, 900, '+'): 'baktaGene_short',
            (150, 950, '+'): 'Rv0001',
        }
        _insert_with_dedup(d, (120, 920, '+'), 'newGene', reference_ids={'Rv0001'})
        assert (120, 920, '+') not in d          # candidate rejected
        assert set(d.values()) == {'baktaGene_short', 'Rv0001'}

    def test_similar_call_does_not_displace_a_reference_gene(self):
        """Bakta's version of Rv1614, one codon longer, replaced the reference
        gene: telling reference genes apart relied on coordinates in the ID,
        which clean_id had already removed."""
        d = {(1000, 2500, '+', 'c1'): 'Rv1614'}
        _insert_with_dedup(d, (997, 2500, '+', 'c1'), 'S1_01690', reference_ids={'Rv1614'})
        assert d == {(1000, 2500, '+', 'c1'): 'Rv1614'}

    def test_a_much_longer_call_still_displaces_a_reference_gene(self):
        d = {(1000, 2400, '+', 'c1'): 'Rv0001'}
        _insert_with_dedup(d, (449, 2400, '+', 'c1'), 'S1_00001', reference_ids={'Rv0001'})
        assert d == {(449, 2400, '+', 'c1'): 'S1_00001'}


# =============================================================================
# RECONCILIATION MODULE
# =============================================================================

class TestDeduplicate:
    """deduplicate must never drop a gene entirely."""

    def test_keeps_one_entry_per_gene(self):
        from classify.reconciliation import deduplicate
        combined = {
            (100, 200, '+'): 'geneA',
            (100, 200, '-'): 'geneA',       # same gene, both fully inside geneB
            (50, 500, '+'): 'geneB',
        }
        deduplicate(combined)
        assert 'geneA' in combined.values()
        assert 'geneB' in combined.values()

    def test_no_duplicates_is_a_noop(self):
        from classify.reconciliation import deduplicate
        combined = {(100, 200, '+'): 'geneA', (300, 400, '+'): 'geneB'}
        before = dict(combined)
        deduplicate(combined)
        assert combined == before


class TestRecoverMissingGenes:
    """A reference protein goes back in only when no gene accounts for it."""

    @staticmethod
    def _mp(rows):
        return pd.DataFrame(rows, columns=['Prot Id', 'Start', 'End', 'Strand', 'Match Nt',
                                           'Num nt', 'Contig'])

    def test_placement_a_codon_off_a_call_is_not_added_again(self):
        """Genes matched in case 2 and 3 are stored under Bakta locus tags, so
        every protein looked missing by name, and only an exact-coordinate
        check kept it out: a placement a codon away went in as a second gene."""
        from classify.reconciliation import recover_missing_genes
        combined = {(1000, 1900, '+', 'c1'): 'S1_00001'}
        recover_missing_genes(combined, self._mp([['Rv0001', 1003, 1900, '+', 897, 897, 'c1']]))
        assert combined == {(1000, 1900, '+', 'c1'): 'S1_00001'}

    def test_a_gene_in_another_frame_accounts_for_nothing(self):
        """Genes that really overlap sit in different frames."""
        from classify.reconciliation import recover_missing_genes
        combined = {(1000, 1900, '+', 'c1'): 'S1_00001'}
        recover_missing_genes(combined, self._mp([['Rv0002', 1001, 1600, '+', 600, 600, 'c1']]))
        assert combined[(1001, 1600, '+', 'c1')] == 'Rv0002'

    def test_another_reference_gene_in_the_same_frame_accounts_for_nothing(self):
        """E. coli nests hokC (b4412) inside mokC (b0018), in one frame."""
        from classify.reconciliation import recover_missing_genes
        mp = self._mp([['b0018', 16751, 16960, '-', 210, 210, 'c1'],
                       ['b4412', 16751, 16903, '-', 153, 153, 'c1']])
        combined = {(16751, 16960, '-', 'c1'): 'b0018'}
        recover_missing_genes(combined, mp)
        assert combined[(16751, 16903, '-', 'c1')] == 'b4412'

        # The same with mokC under the locus tag of the Bakta call that matched it.
        combined = {(16751, 16960, '-', 'c1'): 'S1_00018'}
        recover_missing_genes(combined, mp)
        assert combined[(16751, 16903, '-', 'c1')] == 'b4412'

    def test_a_missing_protein_goes_in_once_where_it_fits_best(self):
        from classify.reconciliation import recover_missing_genes
        combined = {}
        recover_missing_genes(combined, self._mp([['Rv0003', 5000, 5900, '+', 850, 900, 'c1'],
                                                  ['Rv0003', 9000, 9900, '+', 900, 900, 'c1']]))
        assert combined == {(9000, 9900, '+', 'c1'): 'Rv0003'}


# =============================================================================
# INTEGRATION-STYLE TESTS
# =============================================================================

class TestCoordsDict:
    """Tests for the combined_dict key structure."""
    
    def test_strand_aware_keys(self):
        """Verify that (Start, End, Strand) keys prevent cross-strand collisions."""
        from classify_cases import create_coords_dict
        df = pd.DataFrame({
            'Start': [100, 100],
            'End': [200, 200],
            'Strand': ['+', '-'],
            'Locus Tag': ['gene_fwd', 'gene_rev']
        })
        result = create_coords_dict(df, 'Locus Tag')
        assert len(result) == 2  # Both genes preserved
        assert result[(100, 200, '+', '')] == 'gene_fwd'
        assert result[(100, 200, '-', '')] == 'gene_rev'

    def test_contig_aware_keys(self):
        """Every contig of a draft starts at 1: the same coordinates on two
        contigs are two genes."""
        from classify_cases import create_coords_dict
        df = pd.DataFrame({'Start': [100, 100], 'End': [200, 200], 'Strand': ['+', '+'],
                           'Contig': ['c1', 'c2'], 'Locus Tag': ['g1', 'g2']})
        result = create_coords_dict(df, 'Locus Tag')
        assert result == {(100, 200, '+', 'c1'): 'g1', (100, 200, '+', 'c2'): 'g2'}
    
    def test_duplicate_coords_first_wins(self):
        from classify_cases import create_coords_dict
        df = pd.DataFrame({
            'Start': [100, 100],
            'End': [200, 200],
            'Strand': ['+', '+'],
            'Locus Tag': ['first', 'second']
        })
        result = create_coords_dict(df, 'Locus Tag')
        assert len(result) == 1
        assert result[(100, 200, '+', '')] == 'first'


# =============================================================================
# VALIDATION TESTS
# =============================================================================

class TestValidateInputs:
    """Tests for validate_inputs.py."""
    
    def test_parse_fasta(self):
        from validate_inputs import parse_fasta
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">seq1\nATGCATGC\n>seq2\nGGGGCCCC\n")
            f.flush()
            records = list(parse_fasta(f.name))
        os.unlink(f.name)
        assert len(records) == 2
        assert records[0] == ('seq1', 'ATGCATGC')
        assert records[1] == ('seq2', 'GGGGCCCC')
    
    def test_calc_n50(self):
        from validate_inputs import calc_n50
        assert calc_n50([1000, 500, 200, 100]) == 1000
        assert calc_n50([500, 500]) == 500
        assert calc_n50([]) == 0
    
    def test_validate_fasta_good(self):
        from validate_inputs import validate_fasta
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">chr1\n" + "ATGC" * 250000 + "\n")  # 1 Mb genome
            f.flush()
            errors, warnings = validate_fasta(f.name, 'test_sample')
        os.unlink(f.name)
        assert len(errors) == 0
    
    def test_validate_fasta_too_small(self):
        from validate_inputs import validate_fasta
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">chr1\nATGC\n")  # 4 bp genome
            f.flush()
            errors, warnings = validate_fasta(f.name, 'test_sample')
        os.unlink(f.name)
        assert any('too small' in e for e in errors)
    
    def test_validate_fasta_high_n_content(self):
        from validate_inputs import validate_fasta
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            seq = "N" * 150000 + "ATGC" * 100000  # >25% N
            f.write(f">chr1\n{seq}\n")
            f.flush()
            errors, warnings = validate_fasta(f.name, 'test_sample')
        os.unlink(f.name)
        assert any('N-content' in e for e in errors)
    
    def test_validate_fasta_empty(self):
        from validate_inputs import validate_fasta
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write("")
            f.flush()
            errors, warnings = validate_fasta(f.name, 'test_sample')
        os.unlink(f.name)
        assert any('empty' in e.lower() for e in errors)


# =============================================================================
# EXPORT CLEAN PROTEOME TESTS
# =============================================================================

class TestCleanProteome:
    """Tests for export_clean_proteome.py."""
    
    def test_strips_internal_stops(self):
        from export_clean_proteome import clean_proteome
        with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as faa:
            faa.write(">good_gene\nMKTLLFGAA\n>bad_gene\nMKT*LLFGAA\n")
            faa.flush()
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as mt:
            mt.write("simple_id\tIs_Pseudogene\ngood_gene\tFalse\nbad_gene\tFalse\n")
            mt.flush()
        outf = tempfile.mktemp(suffix='.faa')
        clean_proteome(faa.name, mt.name, outf)
        with open(outf) as f:
            content = f.read()
        os.unlink(faa.name)
        os.unlink(mt.name)
        os.unlink(outf)
        assert '>good_gene' in content
        assert '>bad_gene' not in content
    
    def test_strips_pseudogenes(self):
        from export_clean_proteome import clean_proteome
        with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as faa:
            faa.write(">real_gene\nMKTLLFGAA\n>pseudo_gene\nMKTLLFGAA\n")
            faa.flush()
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as mt:
            mt.write("simple_id\tIs_Pseudogene\nreal_gene\tFalse\npseudo_gene\tTrue\n")
            mt.flush()
        outf = tempfile.mktemp(suffix='.faa')
        clean_proteome(faa.name, mt.name, outf)
        with open(outf) as f:
            content = f.read()
        os.unlink(faa.name)
        os.unlink(mt.name)
        os.unlink(outf)
        assert '>real_gene' in content
        assert '>pseudo_gene' not in content
    
    def test_strips_trailing_stop(self):
        from export_clean_proteome import clean_proteome
        with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as faa:
            faa.write(">gene1\nMKTLLFGAA*\n")
            faa.flush()
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as mt:
            mt.write("simple_id\tIs_Pseudogene\ngene1\tFalse\n")
            mt.flush()
        outf = tempfile.mktemp(suffix='.faa')
        clean_proteome(faa.name, mt.name, outf)
        with open(outf) as f:
            content = f.read()
        os.unlink(faa.name)
        os.unlink(mt.name)
        os.unlink(outf)
        assert 'MKTLLFGAA' in content
        assert '*' not in content


# =============================================================================
# ANNOTATION QC TESTS
# =============================================================================

class TestAnnotationQC:
    """Tests for annotation_qc.py."""
    
    def test_parse_fasta_stats(self):
        from annotation_qc import parse_fasta_stats
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            # ATGCGCGC = 8 chars, 50% GC (A,T,G,C each 2x per repeat)
            seq = "ATGCGCGC" * 125000  # 1,000,000 bp
            f.write(f">chr1\n{seq}\n")
            f.flush()
            stats = parse_fasta_stats(f.name)
        os.unlink(f.name)
        assert stats['genome_size'] == 1000000
        assert stats['n_contigs'] == 1
        assert stats['n50'] == 1000000
        # ATGCGCGC: G=2, C=2 out of 8 = 50% GC... but wait:
        # A=1, T=1, G=2, C=2 → GC = 4/6? No: ATGCGCGC → A=1,T=1,G=3,C=3 per 8
        # Actually: A-T-G-C-G-C-G-C → A=1, T=1, G=3, C=3 → GC = 6/8 = 75%
        assert 74 < stats['gc_pct'] < 76

    def test_annotation_rate_ignores_emapper_dashes(self):
        """emapper writes '-' for fields it could not fill.

        Those survive into the master table, so counting with notna() alone
        scored every unannotated gene as annotated and the rate sat at 100%.
        """
        from annotation_qc import compute_qc
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as fa:
            fa.write(">chr1\n" + "ATGC" * 250000 + "\n")
            fa.flush()
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as mt:
            mt.write("simple_id\tCOG_category\tIs_Pseudogene\tstart\tend\n")
            for i in range(10):
                cog = 'S' if i < 4 else '-'      # 4 of 10 genuinely annotated
                mt.write(f"g{i}\t{cog}\tFalse\t{i * 1000 + 1}\t{i * 1000 + 900}\n")
            mt.flush()
        outf = tempfile.mktemp(suffix='.qc.tsv')
        compute_qc(mt.name, fa.name, outf)
        qc = pd.read_csv(outf, sep='\t')
        os.unlink(fa.name)
        os.unlink(mt.name)
        os.unlink(outf)
        assert qc['eggnog_annotation_pct'].iloc[0] == 40.0
        assert qc['total_cds'].iloc[0] == 10

    def test_annotation_rate_counts_real_annotations(self):
        """With nothing missing the rate must still reach 100%."""
        from annotation_qc import compute_qc
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as fa:
            fa.write(">chr1\n" + "ATGC" * 250000 + "\n")
            fa.flush()
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as mt:
            mt.write("simple_id\tCOG_category\tIs_Pseudogene\tstart\tend\n")
            for i in range(5):
                mt.write(f"g{i}\tS\tFalse\t{i * 1000 + 1}\t{i * 1000 + 900}\n")
            mt.flush()
        outf = tempfile.mktemp(suffix='.qc.tsv')
        compute_qc(mt.name, fa.name, outf)
        qc = pd.read_csv(outf, sep='\t')
        os.unlink(fa.name)
        os.unlink(mt.name)
        os.unlink(outf)
        assert qc['eggnog_annotation_pct'].iloc[0] == 100.0


# =============================================================================
# REFERENCE PROTEIN EXTRACTION TESTS
# =============================================================================

class TestAppendStopCodon:
    """The reference proteome has to carry the stop codon.

    The bundled H37Rv file does: 3 * len(protein) equals the coordinate span
    and the sequence ends in '*'. A protein written without it aligns three
    bases short of the CDS Bakta called, so no gene can ever reach case 2.
    """

    def test_appends_when_span_is_one_codon_longer(self):
        from extract_ref_proteins import append_stop_codon
        # 1-36 is 36 nt, 12 codons; an 11 residue protein leaves one spare
        assert append_stop_codon('M' * 11, 1, 36) == 'M' * 11 + '*'

    def test_leaves_matching_lengths_alone(self):
        from extract_ref_proteins import append_stop_codon
        # 12 residues already account for all 36 nt
        assert append_stop_codon('M' * 12, 1, 36) == 'M' * 12

    def test_does_not_double_up_an_existing_stop(self):
        from extract_ref_proteins import append_stop_codon
        assert append_stop_codon('M' * 11 + '*', 1, 36) == 'M' * 11 + '*'

    def test_leaves_other_mismatches_alone(self):
        """Joined or frameshifted features span more than their translation."""
        from extract_ref_proteins import append_stop_codon
        assert append_stop_codon('M' * 11, 1, 300) == 'M' * 11

    def test_handles_empty(self):
        from extract_ref_proteins import append_stop_codon
        assert append_stop_codon('', 1, 36) == ''


# =============================================================================
# CHECK CDS TESTS
# =============================================================================

class TestCheckCds:
    """Tests for check_cds.py."""

    GENOME = "ATG" + "GCT" * 11 + "TAA" + "GCT" * 20

    @staticmethod
    def _run(tmpdir, gene_id, rvs, start=1, end=36, genome=None, proteome=None):
        """Run check_cds.py on one gene and return (start, end, pseudogene).

        The default genome holds a 36 bp gene with no stop codon and a TAA just
        past it, so whether the gene is left alone or extended shows up as
        keeping 1-36 or ending at 39. rvs carries plain reference names, as
        classify_cases.py writes them; whatever check_cds knows about those
        references comes from the proteome, when one is given.
        """
        import subprocess
        fasta = os.path.join(tmpdir, 'g.fasta')
        with open(fasta, 'w') as f:
            f.write(">chr1\n" + (genome or TestCheckCds.GENOME) + "\n")
        tsv = os.path.join(tmpdir, 'sample_corrected.tsv')
        with open(tsv, 'w') as f:
            f.write(f"{gene_id}\t+\t{start}\t{end}\t{rvs}\n")
        cmd = [sys.executable, os.path.join(os.path.dirname(__file__), '..', 'bin', 'check_cds.py'),
               fasta, tsv]
        if proteome:
            faa = os.path.join(tmpdir, 'ref.faa')
            with open(faa, 'w') as f:
                f.writelines(f">{h}\n{seq}\n" for h, seq in proteome.items())
            cmd += ['--ref-proteome', faa]
        subprocess.run(cmd, cwd=tmpdir, check=True, capture_output=True)
        out = pd.read_csv(os.path.join(tmpdir, 'sample.AMAP.tsv'), sep='\t')
        return int(out.iloc[0]['start']), int(out.iloc[0]['end']), out.iloc[0]['pseudogene']

    # 12 residues: the 36 bp gene translated, no stop.
    STOPLESS = "MAAAAAAAAAAA"

    def test_reference_pseudogene_kept_for_mtbc(self):
        """A reference protein without a stop marks a reference pseudogene."""
        import tempfile as tf
        with tf.TemporaryDirectory() as d:
            got = self._run(d, 'S1_Rv0001', 'Rv0001', proteome={'Rv0001_1-36': self.STOPLESS})
            assert got == (1, 36, 'yes')

    def test_reference_pseudogene_kept_for_other_organisms(self):
        """Nothing depends on the gene being called RvNNNN."""
        import tempfile as tf
        with tf.TemporaryDirectory() as d:
            got = self._run(d, 'S1_thrA', 'thrA', proteome={'thrA_1-36': self.STOPLESS})
            assert got == (1, 36, 'yes')

    def test_reference_pseudogene_kept_with_several_references(self):
        """rvs is comma separated when a gene matches several references."""
        import tempfile as tf
        with tf.TemporaryDirectory() as d:
            got = self._run(d, 'S1_multi', 'Rv0001, Rv0002',
                            proteome={'Rv0001_1-36': self.STOPLESS, 'Rv0002_500-600': self.STOPLESS})
            assert got == (1, 36, 'yes')

    def test_fragment_of_a_frameshifted_reference_runs_to_its_stop(self):
        """E. coli b0240 is a frameshifted pseudogene: translated straight
        through, its reference protein has no final stop, because the real one
        is in another frame. The C-terminal fragment a sample carries in that
        frame still ends at a stop, and is extended to it as before."""
        import tempfile as tf
        with tf.TemporaryDirectory() as d:
            # The reference protein ends in residues the fragment never reaches.
            got = self._run(d, 'S1_b0240', 'b0240',
                            proteome={'b0240_1-200': 'MAAAAAAAAAAA' + 'WRWRWRWRWR'})
            assert got[:2] == (1, 39)

    def test_join_cds_without_a_final_star_still_runs_to_its_stop(self):
        """For a frameshifted (join) CDS the extracted protein has no '*' even
        though the gene has a stop: b0240 is 133 residues over 1178 bp. Its
        last residues match the gene's, so only the span says it is not a
        reference annotated without a stop."""
        import tempfile as tf
        with tf.TemporaryDirectory() as d:
            got = self._run(d, 'S1_b0240', 'b0240', proteome={'b0240_1-1178': self.STOPLESS})
            assert got[:2] == (1, 39)

    def test_reference_with_a_stop_means_the_sample_lost_it(self):
        """If the reference protein has its stop, the sample's gene ran past it
        and is extended to the next one, as it always was."""
        import tempfile as tf
        with tf.TemporaryDirectory() as d:
            got = self._run(d, 'S1_Rv0001', 'Rv0001', proteome={'Rv0001_1-39': self.STOPLESS + '*'})
            assert got[:2] == (1, 39)

    def test_rvs_never_carries_coordinates(self):
        """The pipeline writes plain names into rvs, so without the proteome
        there is nothing to recognise a reference pseudogene by.

        The previous check looked for coordinates in rvs and was tested with
        hand-written ones; classify_cases.py had stopped writing them in March
        2026, so in a real run it never fired.
        """
        import tempfile as tf
        with tf.TemporaryDirectory() as d:
            assert self._run(d, 'S1_Rv0001', 'Rv0001')[:2] == (1, 39)

    def test_gene_without_reference_coordinates_is_extended(self):
        """No reference means nothing vouches for the boundary."""
        import tempfile as tf
        with tf.TemporaryDirectory() as d:
            start, end, _ = self._run(d, 'S1_00042', '', proteome={'Rv0001_1-36': self.STOPLESS})
            assert (start, end) == (1, 39)

    def test_partial_final_codon_of_a_reference_pseudogene_is_kept(self):
        """oxyR' (Rv2427A) spans 517 bp in the reference; its alignment covers
        the 516 that make whole codons. The trailing base goes back on."""
        import tempfile as tf
        with tf.TemporaryDirectory() as d:
            # Reference span 1-37 carries one base past the twelve codons.
            got = self._run(d, 'S1_Rv2427A', 'Rv2427A',
                            proteome={'Rv2427A_1-37': self.STOPLESS + 'A*'})
            assert got == (1, 37, 'yes')

    def test_clipped_n_terminus_is_restored(self):
        """Miniprot aligned Rv3611 from residue 30 on CP110619."""
        import tempfile as tf
        genome = "ATG" + "AAA" * 5 + "GCT" * 10 + "TAA" + "GCT" * 10
        with tf.TemporaryDirectory() as d:
            got = self._run(d, 'S1_Rv3611', 'Rv3611', start=19, end=51, genome=genome,
                            proteome={'Rv3611_1-51': 'MKKKKK' + 'A' * 10 + '*'})
            assert got == (1, 51, 'no')

    def test_n_terminus_is_not_restored_across_a_stop(self):
        """If the sample has a stop where the reference has residues, its
        protein really is shorter."""
        import tempfile as tf
        genome = "ATG" + "AAA" * 2 + "TAA" + "AAA" * 2 + "GCT" * 10 + "TAA" + "GCT" * 10
        with tf.TemporaryDirectory() as d:
            got = self._run(d, 'S1_Rv3611', 'Rv3611', start=19, end=51, genome=genome,
                            proteome={'Rv3611_1-51': 'MKKKKK' + 'A' * 10 + '*'})
            assert got[:2] == (19, 51)

    def test_n_terminus_is_not_restored_when_upstream_differs(self):
        """The gene is clipped, but the sample's upstream codons translate to
        something else than the reference's first residues: its N-terminus is
        not the reference's, so nothing is borrowed from it."""
        import tempfile as tf
        genome = "ATG" + "GAA" * 5 + "GCT" * 10 + "TAA" + "GCT" * 10
        with tf.TemporaryDirectory() as d:
            got = self._run(d, 'S1_Rv3611', 'Rv3611', start=19, end=51, genome=genome,
                            proteome={'Rv3611_1-51': 'MKKKKK' + 'A' * 10 + '*'})
            assert got[:2] == (19, 51)

    def test_internal_deletion_is_not_mistaken_for_a_clipped_start(self):
        """M. bovis Rv3159c (PPE53) is one residue short through an in-frame
        deletion and starts at the same ATG as H37Rv. With a CTG in front, an
        earlier version of the completion put it on, giving LMNYSV..."""
        import tempfile as tf
        genome = "CTG" + "ATG" + "GCT" * 20 + "TAA" + "GCT" * 10
        with tf.TemporaryDirectory() as d:
            got = self._run(d, 'S1_Rv3159c', 'Rv3159c', start=4, end=69, genome=genome,
                            proteome={'Rv3159c_1-72': 'MD' + 'A' * 20 + '*'})
            assert got[:2] == (4, 69)

    def test_internal_deletion_with_a_start_codon_upstream(self):
        """Same, with an ATG in front, which translates to the reference's
        first residue: only the gene already lining up with the reference's
        start says there is nothing to restore."""
        import tempfile as tf
        genome = "ATG" + "ATG" + "GCT" * 20 + "TAA" + "GCT" * 10
        with tf.TemporaryDirectory() as d:
            got = self._run(d, 'S1_Rv3159c', 'Rv3159c', start=4, end=69, genome=genome,
                            proteome={'Rv3159c_1-72': 'MD' + 'A' * 20 + '*'})
            assert got[:2] == (4, 69)

    def test_far_stop_is_found_when_the_reference_says_how_far(self):
        """The last 128 residues of a PE_PGRS went unaligned, 387 bp from the
        stop; the default search stops at 300 bp."""
        import tempfile as tf
        genome = "ATG" + "GCT" * 150 + "TAA" + "GCT" * 20
        protein = 'M' + 'A' * 150 + '*'
        with tf.TemporaryDirectory() as d:
            assert self._run(d, 'S1_g', 'gX', end=90, genome=genome,
                             proteome={'gX_1-456': protein})[:2] == (1, 456)
        with tf.TemporaryDirectory() as d:
            assert self._run(d, 'S1_g', 'gX', end=90, genome=genome)[:2] == (1, 90)

    def test_fasta_matches_the_table_after_deduplication(self):
        """Two genes that extend onto the same coordinates leave one row.

        The sequence dictionaries are filled during the loop, before that row
        is dropped, so the .faa used to keep a protein the .tsv no longer
        listed.
        """
        import subprocess
        import tempfile as tf
        with tf.TemporaryDirectory() as d:
            genome = os.path.join(d, 'g.fasta')
            with open(genome, 'w') as f:
                # 1-33 and 1-36 both run on to the stop at 37-39
                f.write(">chr1\n" + "ATG" + "GCT" * 11 + "TAA" + "GCT" * 20 + "\n")
            tsv = os.path.join(d, 'sample_corrected.tsv')
            with open(tsv, 'w') as f:
                f.write("S1_geneA\t+\t1\t33\t\n")
                f.write("S1_geneB\t+\t1\t36\t\n")
            script = os.path.join(os.path.dirname(__file__), '..', 'bin', 'check_cds.py')
            subprocess.run([sys.executable, script, genome, tsv],
                           cwd=d, check=True, capture_output=True)

            table = pd.read_csv(os.path.join(d, 'sample.AMAP.tsv'), sep='\t')
            with open(os.path.join(d, 'sample.AMAP.faa')) as f:
                faa = {l[1:].split(';')[0] for l in f if l.startswith('>')}
            with open(os.path.join(d, 'sample.AMAP.fna')) as f:
                fna = {l[1:].split(';')[0] for l in f if l.startswith('>')}

            assert set(table['id']) == faa
            assert faa == fna


# =============================================================================
# DEDUPLICATE RVS TESTS
# =============================================================================

class TestDeduplicateRvs:
    """Tests for deduplicate_rvs.py."""

    @staticmethod
    def _gene(id_, start, end, rvs, strand='+', pseudo='no', contig='c1'):
        return {'id': id_, 'strand': strand, 'start': str(start), 'end': str(end),
                'rvs': rvs, 'pseudogene': pseudo, 'contig': contig}

    def test_tie_goes_to_the_copy_bakta_also_called(self):
        """Between copies that tie, one Bakta called at the same coordinates
        wins: it used to be whichever came first in the genome."""
        from deduplicate_rvs import deduplicate
        genes = [self._gene('S1_g1', 266553, 266789, 'g1', '-', pseudo='yes'),
                 self._gene('S1_02737', 2770445, 2770681, 'g1', pseudo='yes')]
        kept, _ = deduplicate(genes, {(2770445, 2770681, '+', 'c1')})
        assert [g['id'] for g in kept] == ['S1_02737']

    def test_bakta_support_only_breaks_ties(self):
        """A longer copy still wins over a Bakta-confirmed shorter one."""
        from deduplicate_rvs import deduplicate
        genes = [self._gene('S1_a', 100, 399, 'g1'), self._gene('S1_b', 1000, 1201, 'g1')]
        kept, _ = deduplicate(genes, {(1000, 1201, '+', 'c1')})
        assert [g['id'] for g in kept] == ['S1_a']

    def test_further_copies_of_a_gene_are_kept(self):
        """E. coli Sakai carries K-12 prophage genes several times over, and
        every copy but one was deleted: genes missing from the annotation and
        from the duplication analysis, which is about exactly these."""
        from deduplicate_rvs import deduplicate
        genes = [self._gene('Sakai_01184', 1184948, 1185481, 'b1554', '-'),
                 self._gene('Sakai_01273', 1273009, 1273542, 'b1554'),
                 self._gene('Sakai_b1554', 1558513, 1559046, 'b1554')]
        kept, removed = deduplicate(genes)
        assert len(kept) == 3 and not removed

    def test_pseudogenes_and_pieces_are_not_copies(self):
        """A further copy must be intact, nearly as long as the best one, and
        clear of the copies kept; anything else is a remnant of the gene."""
        from deduplicate_rvs import deduplicate
        genes = [self._gene('S1_best', 1000, 1899, 'g1'),
                 self._gene('S1_pseudo', 5000, 5899, 'g1', pseudo='yes'),
                 self._gene('S1_piece', 9000, 9599, 'g1'),
                 self._gene('S1_overlapping', 1300, 2199, 'g1'),
                 self._gene('S1_copy', 20000, 20869, 'g1')]
        kept, _ = deduplicate(genes)
        assert [g['id'] for g in kept] == ['S1_best', 'S1_copy']

    def test_a_copy_on_another_reference_gene_is_not_kept(self):
        """In M. bovis an intact copy of Rv2177c lies inside Rv2424c, a
        pseudogene at its own locus."""
        from deduplicate_rvs import deduplicate
        genes = [self._gene('S1_02275', 2422727, 2423392, 'Rv2177c', '-'),
                 self._gene('S1_Rv2177c', 2693372, 2694037, 'Rv2177c', '-'),
                 self._gene('S1_Rv2424c', 2693372, 2694371, 'Rv2424c', '-', pseudo='yes')]
        kept, _ = deduplicate(genes)
        assert [g['id'] for g in kept] == ['S1_02275', 'S1_Rv2424c']

        # A gene on the other strand is another locus.
        genes[2]['strand'] = '+'
        kept, _ = deduplicate(genes)
        assert len(kept) == 3

    def test_copies_on_other_contigs_do_not_overlap(self):
        from deduplicate_rvs import deduplicate
        genes = [self._gene('S1_g1', 1000, 1899, 'g1', contig='c1'),
                 self._gene('S1_g1', 1000, 1899, 'g1', contig='c2')]
        kept, _ = deduplicate(genes)
        assert len(kept) == 2

    def test_synteny_breaks_the_remaining_ties(self):
        """In CDC1551 an 81-residue copy of Rv3022c also opens another PPE
        gene, 2.4 kb before its own locus, and that one came first."""
        from deduplicate_rvs import deduplicate
        genes = [self._gene('S1_Rv3022c', 3373737, 3373979, 'Rv3022c', '-', pseudo='yes'),
                 self._gene('S1_Rv3022c', 3376177, 3376419, 'Rv3022c', '-', pseudo='yes')]
        kept, _ = deduplicate(genes, distance={0: 2440.0, 1: 0.0})
        assert [g['start'] for g in kept] == ['3376177']

    def test_synteny_distance_follows_the_neighbours(self):
        """The sample is 5 kb ahead of the reference here; the copy that sits
        5 kb after the protein's reference position is the syntenic one."""
        from deduplicate_rvs import synteny_distance
        ref = {f'g{i}': (i * 1000, i * 1000 + 899) for i in range(1, 11)}
        genes = [self._gene(f'S1_g{i}', i * 1000 + 5000, i * 1000 + 5899, f'g{i}')
                 for i in range(1, 11) if i != 5]
        genes += [self._gene('S1_g5', 5000 + 5000, 5899 + 5000, 'g5'),
                  self._gene('S1_g5', 5000 + 40000, 5899 + 40000, 'g5')]
        dist = synteny_distance(genes, ref)
        assert dist[len(genes) - 2] < 100 < 30000 < dist[len(genes) - 1]

    def test_copies_get_ids_of_their_own(self):
        """Copies only miniprot placed are all named after their protein, and
        the master table keeps one row per id."""
        from deduplicate_rvs import unique_ids
        genes = [self._gene('S1_b4570', 1, 393, 'b4570'), self._gene('S1_b4570_2', 900, 1292, 'x'),
                 self._gene('S1_b4570', 2000, 2392, 'b4570'), self._gene('S1_b4570', 3000, 3392, 'b4570')]
        assert unique_ids(genes) == {2: 'S1_b4570_3', 3: 'S1_b4570_4'}

    def test_table_and_fasta_keep_the_same_ids(self):
        import subprocess
        import tempfile as tf
        with tf.TemporaryDirectory() as d:
            tsv, faa = os.path.join(d, 's.AMAP.tsv'), os.path.join(d, 's.AMAP.faa')
            with open(tsv, 'w') as f:
                f.write("id\tstrand\tstart\tend\trvs\tpseudogene\tcontig\n")
                f.write("S1_b4570\t+\t1\t393\tb4570\tno\tc1\n")
                f.write("S1_b4570\t-\t2000\t2392\tb4570\tno\tc1\n")
            with open(faa, 'w') as f:
                f.write(">S1_b4570;+;0-393\nMKT\n>S1_b4570;-;1999-2392\nMKT\n")
            script = os.path.join(os.path.dirname(__file__), '..', 'bin', 'deduplicate_rvs.py')
            subprocess.run([sys.executable, script, tsv, '--faa', faa], check=True, capture_output=True)
            table = pd.read_csv(tsv, sep='\t')
            with open(faa) as f:
                headers = [l[1:].split(';')[0] for l in f if l.startswith('>')]
        assert list(table['id']) == ['S1_b4570', 'S1_b4570_2'] == headers

    def test_removing_one_copy_keeps_the_others(self):
        """A multi-copy gene has one record per placement, all sharing an id.

        filter_fasta cut the header at the first ';', so removing one copy took
        every copy of that gene with it.
        """
        from deduplicate_rvs import filter_fasta, record_key
        gene = {'id': 'S1_b4570', 'strand': '+', 'start': '1207083', 'end': '1207478'}
        key = record_key(gene)
        assert key == 'S1_b4570;+;1207082-1207478'

        with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as f:
            f.write(">S1_b4570;+;1207082-1207478\nMKT\n")
            f.write(">S1_b4570;+;1649191-1649588\nMKT\n")
            f.write(">S1_b4570;+;1968105-1968502\nMKT\n")
            f.flush()
        filter_fasta(f.name, {key})
        with open(f.name) as fh:
            kept = [l.strip() for l in fh if l.startswith('>')]
        os.unlink(f.name)
        assert kept == ['>S1_b4570;+;1649191-1649588', '>S1_b4570;+;1968105-1968502']

    def test_missing_file_is_ignored(self):
        from deduplicate_rvs import filter_fasta
        filter_fasta('/nonexistent/path.faa', {'anything'})


# =============================================================================
# MERGE RTABS TESTS
# =============================================================================

class TestMergeRtabs:
    """Tests for merge_rtabs.py."""
    
    def test_merge_different_genes(self):
        from merge_rtabs import merge_rtab_files as merge_rtabs
        with tempfile.NamedTemporaryFile(mode='w', suffix='.Rtab', delete=False) as f1:
            f1.write("Gene\tsample1\ngeneA\t1\ngeneB\t1\n")
            f1.flush()
        with tempfile.NamedTemporaryFile(mode='w', suffix='.Rtab', delete=False) as f2:
            f2.write("Gene\tsample2\ngeneB\t1\ngeneC\t1\n")
            f2.flush()
        outf = tempfile.mktemp(suffix='.tsv')
        merge_rtabs([f1.name, f2.name], outf)
        result = pd.read_csv(outf, sep='\t', index_col=0)
        os.unlink(f1.name)
        os.unlink(f2.name)
        os.unlink(outf)
        assert result.shape == (3, 2)  # 3 genes, 2 samples
        assert result.loc['geneA', 'sample2'] == 0  # filled with 0
        assert result.loc['geneB', 'sample1'] == 1


# =============================================================================
# PAF PARSING TESTS
# =============================================================================

class TestReadPaf:
    """Tests for _read_paf: miniprot PAF has a variable number of tag columns."""

    def _write(self, content):
        f = tempfile.NamedTemporaryFile(mode='w', suffix='.paf', delete=False)
        f.write(content)
        f.close()
        return f.name

    def test_variable_tag_columns(self):
        """Rows with different tag counts must not abort the parse.

        miniprot emits fs:i only for frameshifted alignments, so the widest
        row is not necessarily the first one.
        """
        from classify.io import _read_paf
        path = self._write(
            "Rv0001_1-1524\t507\t0\t507\t+\tchr\t4411532\t0\t1524\t1524\t1524\t60\tAS:i:100\tms:i:100\n"
            "Rv0002_2052-3260\t402\t0\t402\t+\tchr\t4411532\t2051\t3260\t1209\t1209\t60"
            "\tAS:i:90\tms:i:90\tfs:i:1\tcg:Z:402M\n"
        )
        df = _read_paf(path)
        os.unlink(path)
        assert len(df) == 2
        assert list(df['Prot Id']) == ['Rv0001_1-1524', 'Rv0002_2052-3260']
        assert list(df['Start']) == [0, 2051]

    def test_complete_means_the_whole_protein_and_its_stop(self):
        """plcD aligns end to end, its final '*' on a stop codon (st:i:1).
        Rv3022c, which the reference annotates without a stop, aligns end to end
        with none (st:i:0), and Rv0002 misses its first 20 residues. Rows
        without fs:i carry st:i one column earlier."""
        from classify.io import read_miniprot
        path = self._write(
            "Rv1755c_1986854-1987696\t281\t0\t281\t-\tchr\t9000\t100\t943\t843\t843\t0"
            "\tAS:i:1491\tms:i:1491\tnp:i:281\tfs:i:0\tst:i:1\tcg:Z:281M\n"
            "Rv3022c_3380440-3380682\t81\t0\t81\t-\tchr\t9000\t2000\t2243\t243\t243\t0"
            "\tAS:i:411\tms:i:411\tnp:i:81\tfs:i:0\tst:i:0\tcg:Z:81M\n"
            "Rv0001_1-903\t301\t0\t301\t+\tchr\t9000\t3000\t3903\t903\t903\t0"
            "\tAS:i:1500\tms:i:1500\tnp:i:301\tst:i:1\tcg:Z:301M\n"
            "Rv0002_2052-3260\t403\t20\t403\t+\tchr\t9000\t5060\t6209\t1149\t1149\t0"
            "\tAS:i:1000\tms:i:1000\tnp:i:383\tst:i:1\tcg:Z:383M\n"
        )
        df = read_miniprot(path)
        os.unlink(path)
        assert dict(zip(df['Prot Id'], df['Complete'])) == {
            'Rv1755c': True, 'Rv3022c': False, 'Rv0001': True, 'Rv0002': False}

    def test_empty_paf(self):
        from classify.io import _read_paf
        path = self._write("")
        df = _read_paf(path)
        os.unlink(path)
        assert len(df) == 0
        assert 'Prot Id' in df.columns

    def test_comment_only_paf(self):
        from classify.io import _read_paf
        path = self._write("# no alignments\n")
        df = _read_paf(path)
        os.unlink(path)
        assert len(df) == 0
        assert 'Prot Id' in df.columns


# =============================================================================
# MASTER TABLE / MYCOBROWSER TESTS
# =============================================================================

class TestGenerateTsvGffGbk:
    """Tests for generate_tsv_gff_gbk.py."""

    def test_myco_columns_are_unique(self):
        """The empty Mycobrowser frame must have a unique 'Locus' merge key.

        name_myco already contains 'Locus'; prepending it again made
        pd.merge(on='Locus') raise and broke every --skip_mycobrowser run.
        """
        from generate_tsv_gff_gbk import columns
        _, _, _, name_myco = columns()
        myco_pd = pd.DataFrame(columns=name_myco)
        assert list(myco_pd.columns).count('Locus') == 1

        merged = pd.DataFrame({'simple_id': ['g1'], 'Locus': ['Rv0001']})
        out = pd.merge(merged, myco_pd, on='Locus', how='left')
        assert len(out) == 1

    def test_read_pandas_skips_header_row(self):
        from generate_tsv_gff_gbk import read_pandas, columns
        _, name_correctedtsv, _, _ = columns()
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("id\tstrand\tstart\tend\trvs\tpseudogene\n")
            f.write("L4_G1_Rv0001\t+\t1\t1524\tRv0001_1-1524\tno\n")
            f.flush()
            df = read_pandas(f.name, name_correctedtsv)
        os.unlink(f.name)
        assert len(df) == 1
        assert df.iloc[0]['id'] == 'L4_G1_Rv0001'
        assert pd.api.types.is_numeric_dtype(df['start'])

    def test_locus_keeps_multi_part_gene_names(self):
        """Splitting rvs on '_' truncated any gene name containing one.

        PE_PGRS1_703830-704057 became 'PE', so every PE/PPE gene of a sample
        collapsed onto one Mycobrowser join key.
        """
        from generate_tsv_gff_gbk import locus_from_rvs
        assert locus_from_rvs('PE_PGRS1_703830-704057') == 'PE_PGRS1'
        assert locus_from_rvs('PPE_MPTR1_100-200') == 'PPE_MPTR1'
        assert locus_from_rvs('ECOLI_thrA_337-2799') == 'ECOLI_thrA'

    def test_locus_unchanged_for_rv_numbers(self):
        """The common MTBC case must be untouched by the fix."""
        from generate_tsv_gff_gbk import locus_from_rvs
        assert locus_from_rvs('Rv0001_1-1524') == 'Rv0001'
        assert locus_from_rvs('Rv0609A_703830-704057') == 'Rv0609A'

    def test_locus_does_not_collide_across_pe_genes(self):
        from generate_tsv_gff_gbk import locus_from_rvs
        genes = ['PE_PGRS1_100-200', 'PE_PGRS2_300-400', 'PE_PGRS3_500-600']
        assert len({locus_from_rvs(g) for g in genes}) == 3

    def test_locus_takes_first_of_several(self):
        """A gene mapping to more than one reference keeps the first, as before."""
        from generate_tsv_gff_gbk import locus_from_rvs
        assert locus_from_rvs('Rv0001_1-1524, Rv0002_2052-3260') == 'Rv0001'

    def test_locus_passes_missing_values_through(self):
        from generate_tsv_gff_gbk import locus_from_rvs
        assert pd.isna(locus_from_rvs(np.nan))
        assert locus_from_rvs('') == ''

    def test_read_pandas_keeps_headerless_rows(self):
        """GET_OUTPUTS strips the header before calling; that row must survive."""
        from generate_tsv_gff_gbk import read_pandas, columns
        _, name_correctedtsv, _, _ = columns()
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("L4_G1_Rv0001\t+\t1\t1524\tRv0001_1-1524\tno\n")
            f.flush()
            df = read_pandas(f.name, name_correctedtsv)
        os.unlink(f.name)
        assert len(df) == 1
        assert df.iloc[0]['id'] == 'L4_G1_Rv0001'


# =============================================================================
# DIAMOND OUTPUT PARSING TESTS
# =============================================================================

class TestGetDuplications:
    """Tests for get_duplications.py."""

    def test_first_hit_is_not_eaten_as_header(self):
        """Diamond runs with --header, so all headers start with '#'.

        Without header=None pandas promotes the first alignment to column
        names and that duplication pair disappears.
        """
        pytest.importorskip('pyvis')
        pytest.importorskip('community')
        from get_duplications import read_and_create
        with tempfile.NamedTemporaryFile(mode='w', suffix='.diamond', delete=False) as f:
            f.write("# Fields: qseqid, sseqid, pident\n")
            for q, s, bits in [('gA', 'gB', 200), ('gC', 'gD', 190), ('gE', 'gF', 180)]:
                f.write(f"{q}\t{s}\t99.0\t100\t0\t0\t1\t100\t1\t100\t1e-50\t{bits}\t100\t100\t99\n")
            f.flush()
            df = read_and_create(f.name)
        os.unlink(f.name)
        assert len(df) == 3
        assert set(df['qseqid']) == {'gA', 'gC', 'gE'}


# =============================================================================
# BIOPYTHON COMPATIBILITY TESTS
# =============================================================================

class TestBiopythonCompat:
    """SeqFeature.strand was removed in Biopython 1.85; use location.strand."""

    def test_to_panaroo_reads_strand_from_location(self):
        pytest.importorskip('Bio')
        from pathlib import Path
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, SimpleLocation
        from to_panaroo import convert_gbk_to_panaroo_gff3

        rec = SeqRecord(Seq("ATG" * 20 + "TAA"), id="chr1", name="chr1", description="")
        rec.annotations['molecule_type'] = 'DNA'
        rec.features = [SeqFeature(SimpleLocation(0, 63, strand=-1), type="CDS",
                                   qualifiers={"locus_tag": ["g1"], "product": ["thing"]})]
        gbk = tempfile.mktemp(suffix='.gbk')
        gff = tempfile.mktemp(suffix='.gff')
        SeqIO.write([rec], gbk, "genbank")
        convert_gbk_to_panaroo_gff3(Path(gbk), Path(gff))
        with open(gff) as f:
            content = f.read()
        os.unlink(gbk)
        os.unlink(gff)
        assert "\tCDS\t1\t63\t.\t-\t" in content
        assert "ID=g1;Parent=g1_gene" in content


if __name__ == '__main__':
    pytest.main([__file__, '-v'])


# =============================================================================
# PARALOG PLACEMENTS (synteny.py)
# =============================================================================

class TestParalogPlacements:
    """drop_paralog_placements keeps each protein at its own locus.

    The frames mimic read_miniprot(..., ref_coords=True): anchors are proteins
    placed once, every 10 kb, and a 2 kb insertion halfway through the sample
    shifts everything after it, so reference coordinates cannot be compared
    with sample positions directly.
    """

    COLS = ['Prot Id', 'Prot len', 'Strand', 'Start', 'End', 'Match Nt', 'Num nt',
            'Ref Start', 'Ref End', 'Contig']

    @classmethod
    def _frame(cls, extra=(), contig='chr', insertion_at=200_000, shift=2_000, n=40):
        rows = []
        for i in range(n):
            ref = 10_000 * i + 1
            pos = ref + (shift if ref > insertion_at else 0)
            rows.append([f'g{i}', 300, '+', pos, pos + 899, 900, 900, ref, ref + 899, contig])
        rows += [list(r) + [contig] if len(r) == 9 else list(r) for r in extra]
        return pd.DataFrame(rows, columns=cls.COLS)

    @staticmethod
    def _placed(df, prot):
        return sorted(df.loc[df['Prot Id'] == prot, 'Start'])

    def test_returns_the_usual_columns(self):
        """-rvmissing.tsv and the report expect exactly these seven."""
        from classify.synteny import drop_paralog_placements
        out = drop_paralog_placements(self._frame())
        assert list(out.columns) == ['Prot Id', 'Prot len', 'Strand', 'Start', 'End',
                                     'Match Nt', 'Num nt', 'Contig']

    def test_paralog_on_another_genes_locus_is_dropped(self):
        """Rv1945 also lands inside Rv1148c; that copy is Rv1148c's locus."""
        from classify.synteny import drop_paralog_placements
        # P belongs at reference 305_001, i.e. sample 307_001 after the insertion,
        # and also aligns inside g5 at 50_001.
        df = self._frame(extra=[
            ['P', 300, '+', 307_101, 307_900, 800, 800, 305_101, 305_900],
            ['P', 300, '+', 50_100, 50_899, 780, 800, 305_101, 305_900],
        ])
        out = drop_paralog_placements(df)
        assert self._placed(out, 'P') == [307_101]
        assert self._placed(out, 'g5') == [50_001]

    def test_offsets_come_from_the_sample_not_the_reference(self):
        """After an insertion, a protein's reference coordinates point at
        another locus of the sample. Compared in absolute terms, a copy sitting
        exactly there would pass for the protein's own; the local offset says
        the protein is 5 kb further on.
        """
        from classify.synteny import drop_paralog_placements
        df = self._frame(shift=5_000)
        df = df[df['Prot Id'] != 'g30']
        extra = pd.DataFrame([
            # Q belongs at reference 295_001, sample 300_001 after the shift.
            ['Q', 300, '+', 300_001, 300_900, 900, 900, 295_001, 295_900, 'chr'],
            # P belongs at reference 300_101, sample 305_101; it also aligns
            # at 300_101, inside Q, where its reference coordinates point.
            ['P', 300, '+', 305_101, 305_800, 700, 700, 300_101, 300_800, 'chr'],
            ['P', 300, '+', 300_100, 300_799, 700, 700, 300_101, 300_800, 'chr'],
        ], columns=self.COLS)
        out = drop_paralog_placements(pd.concat([df, extra], ignore_index=True))
        assert self._placed(out, 'P') == [305_101]
        assert self._placed(out, 'Q') == [300_001]

    def test_identical_copies_annotated_differently_follow_position(self):
        """Rv1150 and Rv1041c: two identical copies of one element, annotated
        as the short protein at one and the long one at the other. Only
        position can tell which is which."""
        from classify.synteny import drop_paralog_placements
        df = self._frame(extra=[
            # short protein S belongs at 120_001, long protein L at 340_001
            ['S', 184, '+', 120_001, 120_552, 552, 552, 120_001, 120_552],
            ['S', 184, '+', 342_001, 342_552, 552, 552, 120_001, 120_552],
            ['L', 288, '+', 120_001, 120_864, 864, 864, 340_001, 340_864],
            ['L', 288, '+', 342_001, 342_864, 864, 864, 340_001, 340_864],
        ])
        out = drop_paralog_placements(df)
        assert self._placed(out, 'S') == [120_001]
        assert self._placed(out, 'L') == [342_001]

    def test_better_match_is_not_dropped_for_its_position(self):
        """E. coli Sakai has two of K-12's three ldr toxin copies. The second
        one sits where b4423 belongs, but it is a b4421, which differs from
        b4423 at two residues; position alone labelled it b4423."""
        from classify.synteny import drop_paralog_placements
        df = self._frame(extra=[
            # b4421-like belongs at 120_001; b4423-like at 125_001. The sample's
            # copy at 125_001 is an exact b4421 and a 94% b4423.
            ['b4421', 36, '+', 120_001, 120_108, 108, 108, 120_001, 120_108],
            ['b4421', 36, '+', 125_001, 125_108, 108, 108, 120_001, 120_108],
            ['b4423', 36, '+', 125_001, 125_108, 102, 108, 125_001, 125_108],
        ])
        out = drop_paralog_placements(df)
        assert self._placed(out, 'b4421') == [120_001, 125_001]
        assert self._placed(out, 'b4423') == [125_001]

    def test_overlapping_genes_in_different_frames_are_kept(self):
        """Rv1041c and Rv1042c really overlap, in different frames. A SNP in
        one of them must not cost it its place to the other."""
        from classify.synteny import drop_paralog_placements
        df = self._frame(extra=[
            ['A', 288, '-', 100_001, 100_864, 864, 864, 100_001, 100_864],
            ['B', 136, '-', 100_521, 100_928, 396, 408, 100_521, 100_928],
        ])
        out = drop_paralog_placements(df)
        assert self._placed(out, 'A') == [100_001]
        assert self._placed(out, 'B') == [100_521]

    def test_neighbours_brought_into_one_frame_are_both_kept(self):
        """Rv2879c and Rv2880c overlap in H37Rv in different frames; an indel
        in CDC1551 puts them in one. They are neighbours in the reference, not
        a paralog squatting on someone else's locus, so neither is dropped."""
        from classify.synteny import drop_paralog_placements
        df = self._frame(extra=[
            ['Rv2879c', 190, '-', 100_001, 100_570, 552, 570, 100_001, 100_570],
            ['Rv2880c', 276, '-', 100_268, 100_828, 558, 561, 100_269, 100_828],
        ])
        out = drop_paralog_placements(df)
        assert self._placed(out, 'Rv2879c') == [100_001]
        assert self._placed(out, 'Rv2880c') == [100_268]

    def test_paralog_next_door_after_a_deletion(self):
        """CDC1551 lacks Rv3427c-Rv3428c; Rv3425 then lands on Rv3429, close
        enough to home to look syntenic. The better match keeps the locus."""
        from classify.synteny import drop_paralog_placements
        df = self._frame(extra=[
            ['Rv3429', 179, '+', 155_001, 155_537, 537, 537, 155_001, 155_537],
            ['Rv3425', 176, '+', 155_001, 155_528, 465, 528, 150_001, 150_528],
        ])
        out = drop_paralog_placements(df)
        assert self._placed(out, 'Rv3429') == [155_001]
        assert self._placed(out, 'Rv3425') == []

    def test_copy_at_a_locus_nobody_claims_is_kept(self):
        """An IS6110 insertion at a site the reference lacks, or a duplication
        in the sample: no protein belongs there, so nothing is dropped."""
        from classify.synteny import drop_paralog_placements
        df = self._frame(extra=[
            ['IS', 300, '+', 150_001, 150_900, 900, 900, 150_001, 150_900],
            ['IS', 300, '+', 255_001, 255_900, 900, 900, 150_001, 150_900],
        ])
        out = drop_paralog_placements(df)
        assert self._placed(out, 'IS') == [150_001, 255_001]

    def test_reverse_oriented_sequence(self):
        """A reverse-complemented assembly maps positions backwards."""
        from classify.synteny import drop_paralog_placements
        L = 500_000
        rows = []
        for i in range(40):
            ref = 10_000 * i + 1
            rows.append([f'g{i}', 300, '-', L - ref - 899, L - ref, 900, 900, ref, ref + 899, 'chr'])
        # P belongs at reference 205_101; its paralog copy sits inside g5.
        own = L - 205_900
        rows.append(['P', 300, '-', own, own + 799, 800, 800, 205_101, 205_900, 'chr'])
        g5 = L - 50_001 - 899
        rows.append(['P', 300, '-', g5 + 100, g5 + 899, 780, 800, 205_101, 205_900, 'chr'])
        out = drop_paralog_placements(pd.DataFrame(rows, columns=self.COLS))
        assert self._placed(out, 'P') == [own]

    def test_contigs_are_kept_apart(self):
        """Offsets from one contig say nothing about another."""
        from classify.synteny import drop_paralog_placements
        a = self._frame(contig='c1')
        b = self._frame(contig='c2', shift=0)
        b['Prot Id'] = 'h' + b['Prot Id']
        extra = pd.DataFrame([
            ['P', 300, '+', 120_101, 120_900, 800, 800, 120_101, 120_900, 'c2'],
            ['P', 300, '+', 50_100, 50_899, 780, 800, 120_101, 120_900, 'c2'],
        ], columns=self.COLS)
        out = drop_paralog_placements(pd.concat([a, b, extra], ignore_index=True))
        assert self._placed(out, 'P') == [120_101]
        assert len(out) == len(a) + len(b) + 1

    def test_unplaced_proteins_pass_through(self):
        """miniprot -u writes a row for every protein it could not place, with
        '*' for target and strand and zeros for coordinates. They used to crash
        the interval tree on CDC1551, where Rv3427c is deleted."""
        from classify.synteny import drop_paralog_placements
        # The three CDC1551 lacks: three of them is what it takes for the '*'
        # rows to look like a sequence with enough anchors to predict from.
        df = self._frame(extra=[
            ['Rv0793', 102, '*', 1, 0, 0, 0, 886_719, 887_024, '*'],
            ['Rv3427c', 252, '*', 1, 0, 0, 0, 3_843_885, 3_844_640, '*'],
            ['Rv3428c', 411, '*', 1, 0, 0, 0, 3_844_738, 3_845_970, '*'],
        ])
        out = drop_paralog_placements(df)
        assert len(out) == len(df)
        assert set(out.loc[out['Strand'] == '*', 'Prot Id']) == {'Rv0793', 'Rv3427c', 'Rv3428c'}

    def test_without_reference_coordinates_nothing_changes(self):
        from classify.synteny import drop_paralog_placements
        df = self._frame(extra=[
            ['P', 300, '+', 307_101, 307_900, 800, 800, None, None],
            ['P', 300, '+', 50_101, 50_900, 780, 800, None, None],
        ])
        df['Ref Start'] = np.nan
        df['Ref End'] = np.nan
        assert len(drop_paralog_placements(df)) == len(df)


class TestStrandOfResolvedHit:
    """Resolving a conflict takes the strand of the hit being resolved."""

    def test_second_placement_keeps_its_own_strand(self):
        """esxM (Rv1792) is placed twice, its paralog first and on the other
        strand; it came out at its own coordinates on the paralog's strand."""
        from classify.overlaps import (build_overlap_tree, create_overlap_groups,
                                       resolve_single_overlaps)
        mp = pd.DataFrame([
            ['Rv1792', 99, '-', 1_160_855, 1_161_151, 294, 297],
            ['Rv1792', 99, '+', 2_030_347, 2_030_643, 297, 297],
        ], columns=['Prot Id', 'Prot len', 'Strand', 'Start', 'End', 'Match Nt', 'Num nt'])
        case4 = pd.DataFrame([['S1_01880', 2_030_401, 2_030_643, '+', 243]],
                             columns=['Locus Tag', 'Start', 'End', 'Strand', 'Length_nt'])
        groups = create_overlap_groups(build_overlap_tree(mp, case4, 'Prot Id', 'Locus Tag'),
                                       case4, mp)
        b_groups = create_overlap_groups(build_overlap_tree(case4, mp, 'Locus Tag', 'Prot Id'),
                                         case4, mp)
        combined = {}
        resolve_single_overlaps(groups, mp, combined, {v['Prot Id']: v for v in b_groups.values()})
        assert combined == {(2_030_347, 2_030_643, '+', ''): 'Rv1792'}

    def test_bakta_coordinates_keep_baktas_strand(self):
        """When the Bakta call wins, it keeps the strand Bakta gave it."""
        from classify.overlaps import (build_overlap_tree, create_overlap_groups,
                                       resolve_single_overlaps)
        mp = pd.DataFrame([['X', 100, '-', 10_000, 10_300, 300, 300]],
                          columns=['Prot Id', 'Prot len', 'Strand', 'Start', 'End', 'Match Nt', 'Num nt'])
        # Much longer than the hit, so the Bakta coordinates are chosen.
        case4 = pd.DataFrame([['S1_00001', 9_000, 10_300, '-', 1_301]],
                             columns=['Locus Tag', 'Start', 'End', 'Strand', 'Length_nt'])
        groups = create_overlap_groups(build_overlap_tree(mp, case4, 'Prot Id', 'Locus Tag'),
                                       case4, mp)
        combined = {}
        resolve_single_overlaps(groups, mp, combined, {})
        assert combined == {(9_000, 10_300, '-', ''): 'S1_00001'}

    def test_opposite_strands_do_not_compete(self):
        """Case 4 is defined strand by strand; the grouping that resolves it was
        not, and about 8% of its groups mixed an alignment on one strand with a
        call on the other."""
        from classify.overlaps import build_overlap_tree
        mp = pd.DataFrame([['X', 100, '-', 10_000, 10_300, 300, 300]],
                          columns=['Prot Id', 'Prot len', 'Strand', 'Start', 'End', 'Match Nt', 'Num nt'])
        case4 = pd.DataFrame([['S1_00001', 9_000, 10_300, '+', 1_301]],
                             columns=['Locus Tag', 'Start', 'End', 'Strand', 'Length_nt'])
        assert build_overlap_tree(mp, case4, 'Prot Id', 'Locus Tag') == {}


class TestAlignmentQuality:
    """Case 4 weighs how much of its protein an alignment covers, not only how
    long it is next to Bakta's call."""

    @staticmethod
    def _resolve(mp, case4):
        from classify.overlaps import (build_overlap_tree, create_overlap_groups,
                                       resolve_single_overlaps, resolve_multiple_overlaps)
        overlaps = build_overlap_tree(mp, case4, 'Prot Id', 'Locus Tag')
        single = create_overlap_groups({k: v for k, v in overlaps.items() if v['count'] == 1},
                                       case4, mp)
        multiple = create_overlap_groups({k: v for k, v in overlaps.items() if v['count'] > 1},
                                         case4, mp)
        b_groups = create_overlap_groups(build_overlap_tree(case4, mp, 'Locus Tag', 'Prot Id'),
                                         case4, mp)
        new_groups = {v['Prot Id']: v for v in b_groups.values()}
        combined = {}
        resolve_single_overlaps(single, mp, combined, new_groups)
        resolve_multiple_overlaps(multiple, mp, combined, new_groups)
        return combined

    @staticmethod
    def _mp(rows):
        return pd.DataFrame(rows, columns=['Prot Id', 'Prot len', 'Strand', 'Start', 'End',
                                           'Match Nt', 'Num nt', 'Complete'])

    @staticmethod
    def _case4(rows):
        return pd.DataFrame(rows, columns=['Locus Tag', 'Start', 'End', 'Strand', 'Length_nt'])

    def test_complete_alignment_beats_a_much_longer_call(self):
        """plcD on CDC1551: Bakta's call ends at the same stop and starts
        702 bp further upstream, and it replaced the gene."""
        mp = self._mp([['Rv1755c', 281, '-', 1_977_868, 1_978_710, 843, 843, True]])
        case4 = self._case4([['S1_01840', 1_977_868, 1_979_412, '-', 1_545]])
        assert self._resolve(mp, case4) == {(1_977_868, 1_978_710, '-', ''): 'Rv1755c'}

    def test_partial_alignment_still_gives_way(self):
        mp = self._mp([['Rv1755c', 281, '-', 1_977_868, 1_978_710, 843, 843, False]])
        case4 = self._case4([['S1_01840', 1_977_868, 1_979_412, '-', 1_545]])
        assert self._resolve(mp, case4) == {(1_977_868, 1_979_412, '-', ''): 'S1_01840'}

    def test_complete_alignment_also_wins_among_several_calls(self):
        """The same alignment touching a second call by its last four bases
        is resolved with the multiple overlaps, under the same rule."""
        mp = self._mp([['Rv1755c', 281, '-', 1_977_868, 1_978_710, 843, 843, True]])
        case4 = self._case4([['S1_01840', 1_977_868, 1_979_412, '-', 1_545],
                             ['S1_01839', 1_976_870, 1_977_871, '-', 1_002]])
        assert self._resolve(mp, case4) == {(1_977_868, 1_978_710, '-', ''): 'Rv1755c'}

    def test_call_spanning_several_reference_genes_does_not_win(self):
        """Rv3022c and Rv3021c are one reading frame in CDC1551, called as one
        1308 bp gene. Rv3022c shares four bases with PE29, so it is resolved
        with the multiple overlaps, where the call won on length and then lost
        to Rv3021c, and Rv3022c was gone."""
        mp = self._mp([['Rv3021c', 359, '-', 3_375_112, 3_376_176, 1_050, 1_065, False],
                       ['Rv3022c', 81, '-', 3_376_177, 3_376_419, 243, 243, False],
                       ['Rv3022A', 105, '-', 3_376_416, 3_376_730, 315, 315, True]])
        case4 = self._case4([['S1_03152', 3_375_112, 3_376_419, '-', 1_308],
                             ['S1_03153', 3_376_416, 3_376_733, '-', 318]])
        assert self._resolve(mp, case4) == {
            (3_375_112, 3_376_176, '-', ''): 'Rv3021c',
            (3_376_177, 3_376_419, '-', ''): 'Rv3022c',
            (3_376_416, 3_376_730, '-', ''): 'Rv3022A'}


# =============================================================================
# CONTIGS
# =============================================================================

class TestContigs:
    """Every contig of a draft assembly starts at 1, so coordinates alone do not
    say which contig a gene is on. Cut into 40 contigs, CDC1551 had 87% of its
    proteins translated from the wrong contig."""

    @staticmethod
    def _bakta(rows):
        return pd.DataFrame(rows, columns=['Locus Tag', 'Start', 'End', 'Strand', 'Length_nt', 'Contig'])

    @staticmethod
    def _mp(rows):
        return pd.DataFrame(rows, columns=['Prot Id', 'Prot len', 'Strand', 'Start', 'End',
                                           'Match Nt', 'Num nt', 'Contig'])

    def test_same_coordinates_on_another_contig_are_not_a_match(self):
        from classify.cases import classify_case2_case3, classify_case1_case4, find_rv_missing
        bk = self._bakta([['S1_00001', 100, 399, '+', 300, 'c2']])
        mp = self._mp([['Rv0001', 100, '+', 100, 399, 300, 300, 'c1']])
        case2, _ = classify_case2_case3(bk, mp)
        case4, case1 = classify_case1_case4(bk, mp)
        assert case2.empty and case4.empty
        assert list(case1['Locus Tag']) == ['S1_00001']
        assert list(find_rv_missing(mp, bk)['Prot Id']) == ['Rv0001']

    def test_genes_on_different_contigs_never_displace_each_other(self):
        from classify.overlaps import _insert_with_dedup
        combined = {(100, 900, '+', 'c1'): 'Rv0001'}
        _insert_with_dedup(combined, (150, 950, '+', 'c2'), 'S1_00002')
        assert combined == {(100, 900, '+', 'c1'): 'Rv0001', (150, 950, '+', 'c2'): 'S1_00002'}

    def test_bakta_names_must_match_the_genome(self):
        """Without --keep-contig-headers Bakta calls every sequence contig_N."""
        from classify_cases import match_bakta_contigs
        one = match_bakta_contigs(self._bakta([['S1_00001', 1, 300, '+', 300, 'contig_1']]),
                                  self._mp([['Rv0001', 100, '+', 1, 300, 300, 300, 'S1']]))
        assert set(one['Contig']) == {'S1'}
        with pytest.raises(SystemExit):
            match_bakta_contigs(
                self._bakta([['S1_00001', 1, 300, '+', 300, 'contig_1'],
                             ['S1_00002', 1, 300, '+', 300, 'contig_2']]),
                self._mp([['Rv0001', 100, '+', 1, 300, 300, 300, 'S1'],
                          ['Rv0002', 100, '+', 1, 300, 300, 300, 'node_2']]))

    def test_check_cds_translates_from_the_genes_own_contig(self):
        """The gene is on c2, but its coordinates also fit the shorter c1, which
        is what the old guess picked."""
        import subprocess
        import tempfile as tf
        with tf.TemporaryDirectory() as d:
            with open(os.path.join(d, 'g.fasta'), 'w') as f:
                f.write(">c1\n" + "GCT" * 20 + "\n")
                f.write(">c2\n" + "ATG" + "GCT" * 11 + "TAA" + "GCT" * 40 + "\n")
            with open(os.path.join(d, 'sample_corrected.tsv'), 'w') as f:
                f.write("S1_00001\t+\t1\t39\tRv0001\tc2\n")
            script = os.path.join(os.path.dirname(__file__), '..', 'bin', 'check_cds.py')
            subprocess.run([sys.executable, script, 'g.fasta', 'sample_corrected.tsv'],
                           cwd=d, check=True, capture_output=True)
            out = pd.read_csv(os.path.join(d, 'sample.AMAP.tsv'), sep='\t')
            with open(os.path.join(d, 'sample.AMAP.faa')) as fh:
                protein = fh.read().split('\n')[1]
        assert out.iloc[0]['pseudogene'] == 'no'
        assert out.iloc[0]['contig'] == 'c2'
        assert list(out.columns)[-1] == 'contig'
        assert protein == 'M' + 'A' * 11

    def test_dedup_overlap_only_counts_on_the_same_contig(self):
        from deduplicate_rvs import deduplicate
        genes = [
            {'id': 'S1_Rv0001', 'strand': '+', 'start': '100', 'end': '900',
             'rvs': 'Rv0001', 'pseudogene': 'no', 'contig': 'c1'},
            {'id': 'S1_00002', 'strand': '+', 'start': '100', 'end': '900',
             'rvs': '', 'pseudogene': 'no', 'contig': 'c2'},
        ]
        kept, _ = deduplicate(genes)
        assert [g['id'] for g in kept] == ['S1_Rv0001', 'S1_00002']

    def test_final_tables_put_each_gene_on_its_contig(self):
        """generate_tsv placed each gene on the smallest contig its coordinates
        fit in: 90% of a 40-contig draft ended up on the wrong one."""
        from generate_tsv_gff_gbk import generate_tsv
        import tempfile as tf
        merged = pd.DataFrame({
            'id': ['S1_00001', 'S1_00002'], 'strand': ['+', '+'], 'start': [1, 1], 'end': [30, 30],
            'rvs': ['Rv0001', 'Rv0002'], 'Preferred_name': ['-', '-'], 'Gene': ['a', 'b'],
            'Product_x': ['p', 'q'], 'Description': ['-', '-'], 'contig': ['big', 'small']})
        rest = pd.DataFrame(columns=["Sequence Id", "Type", "Start", "Stop", "Strand",
                                     "Locus Tag", "Gene", "Product", "DbXrefs"])
        with tf.TemporaryDirectory() as d:
            fasta = os.path.join(d, 'g.fasta')
            with open(fasta, 'w') as f:
                f.write(">big\n" + "A" * 900 + "\n>small\n" + "A" * 60 + "\n")
            out = generate_tsv(merged, rest, fasta, os.path.join(d, 'out.tsv'))
        placed = dict(zip(out['Locus Tag'], out['Sequence Id']))
        assert placed == {'S1_00001': 'big', 'S1_00002': 'small'}
