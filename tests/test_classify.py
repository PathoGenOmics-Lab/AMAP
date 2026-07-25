#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the AMAP classification module.

Tests cover:
- ID cleaning and formatting (ids.py)
- Overlap calculations (overlaps.py)
- Coordinate reconciliation (reconciliation.py)
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
from classify.reconciliation import (
    parse_ref_coords_from_id, _extract_gene_name, _gene_names_match
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
    """Tests for _is_reference_gene — detects embedded coordinates."""
    
    def test_rv_with_coords(self):
        assert _is_reference_gene('Rv0001_1-1524') is True
    
    def test_locus_tag_no_coords(self):
        assert _is_reference_gene('ABCDE_00042') is False
    
    def test_pe_pgrs_with_coords(self):
        assert _is_reference_gene('PE_PGRS1_703830-704057') is True
    
    def test_empty_string(self):
        assert _is_reference_gene('') is False


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
            (150, 950, '+'): 'Rv0001_150-950',
        }
        _insert_with_dedup(d, (120, 920, '+'), 'newGene')
        assert (120, 920, '+') not in d          # candidate rejected
        assert set(d.values()) == {'baktaGene_short', 'Rv0001_150-950'}


# =============================================================================
# RECONCILIATION MODULE
# =============================================================================

class TestParseRefCoords:
    """Tests for parse_ref_coords_from_id."""
    
    def test_simple_coords(self):
        start, end = parse_ref_coords_from_id('Rv0001_1-1524')
        assert start == 1
        assert end == 1524
    
    def test_large_coords(self):
        start, end = parse_ref_coords_from_id('Rv0609A_703830-704057')
        assert start == 703830
        assert end == 704057
    
    def test_no_coords(self):
        result = parse_ref_coords_from_id('Rv0001')
        assert result == (None, None)
    
    def test_prefixed_coords(self):
        start, end = parse_ref_coords_from_id('L8_G26689_Rv0001_1-1524')
        assert start == 1
        assert end == 1524


class TestExtractGeneName:
    """Tests for _extract_gene_name."""
    
    def test_with_coords(self):
        assert _extract_gene_name('Rv0001_1-1524') == 'Rv0001'
    
    def test_without_coords(self):
        assert _extract_gene_name('Rv0001') == 'Rv0001'
    
    def test_prefixed(self):
        # _extract_gene_name uses clean_id which strips coords, then takes last part
        result = _extract_gene_name('L8_G26689_Rv0001')
        assert 'Rv0001' in result


class TestGeneNamesMatch:
    """Tests for _gene_names_match — prevents substring false matches."""
    
    def test_exact_match(self):
        assert _gene_names_match('Rv3508', 'Rv3508') is True
    
    def test_substring_no_match(self):
        """Rv3508 should NOT match Rv3508c."""
        assert _gene_names_match('Rv3508', 'Rv3508c') is False
    
    def test_rv1_rv10_no_match(self):
        """Rv1 should NOT match Rv10."""
        assert _gene_names_match('Rv1', 'Rv10') is False
    
    def test_with_coords(self):
        assert _gene_names_match('Rv0001_1-1524', 'Rv0001') is True
    
    def test_completely_different(self):
        assert _gene_names_match('Rv0001', 'Rv9999') is False


class TestReconcileCoordinates:
    """Behaviour locked in around the gene-name indexed lookups.

    reconcile_coordinates used to scan every protein and every reference for
    each gene; it now indexes both by gene name up front. These cases pin the
    semantics that indexing must preserve.
    """

    @staticmethod
    def _mp(rows):
        return pd.DataFrame(rows, columns=['Prot Id', 'Start', 'End', 'Strand',
                                           'Match Nt', 'Num nt'])

    def test_replaces_with_reference_coordinates(self):
        from classify.reconciliation import reconcile_coordinates
        combined = {(1000, 1900, '+'): 'Rv0001_1010-1910'}
        mp = self._mp([['Rv0001_1010-1910', 1010, 1910, '+', 900, 900]])
        reconcile_coordinates(combined, mp)
        assert combined == {(1010, 1910, '+'): 'Rv0001_1010-1910'}

    def test_does_not_reconcile_similar_names(self):
        """Rv3508 must not be pulled onto Rv3508c's coordinates."""
        from classify.reconciliation import reconcile_coordinates
        combined = {(1000, 1900, '+'): 'Rv3508_1000-1900'}
        mp = self._mp([['Rv3508c_1010-1910', 1010, 1910, '+', 900, 900],
                       ['Rv3508_1000-1900', 1000, 1900, '+', 900, 900]])
        reconcile_coordinates(combined, mp)
        assert combined == {(1000, 1900, '+'): 'Rv3508_1000-1900'}

    def test_multi_copy_element_survives_paralog_check(self):
        """>2 miniprot hits means placement away from the embedded coords is expected."""
        from classify.reconciliation import reconcile_coordinates
        combined = {(500000, 500600, '+'): 'PE_PGRS1_1000-1600'}
        mp = self._mp([['PE_PGRS1_1000-1600', 1000, 1600, '+', 600, 600],
                       ['PE_PGRS1_1000-1600', 500000, 500600, '+', 600, 600],
                       ['PE_PGRS1_1000-1600', 900000, 900600, '+', 600, 600]])
        reconcile_coordinates(combined, mp)
        assert (500000, 500600, '+') in combined

    def test_bakta_validates_position(self):
        """A Bakta call at the exact coordinates rescues a mismatched placement."""
        from classify.reconciliation import reconcile_coordinates
        combined = {(500000, 500600, '+'): 'Rv0001_1000-1600'}
        mp = self._mp([['Rv0001_1000-1600', 500000, 500600, '+', 600, 600]])
        bk = pd.DataFrame([['L4_G1_00001', 500000, 500600, '+', 601]],
                          columns=['Locus Tag', 'Start', 'End', 'Strand', 'Length_nt'])
        without_bakta = dict(combined)
        reconcile_coordinates(without_bakta, mp, None)
        with_bakta = dict(combined)
        reconcile_coordinates(with_bakta, mp, bk)
        assert (500000, 500600, '+') not in without_bakta
        assert (500000, 500600, '+') in with_bakta


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
        assert result[(100, 200, '+')] == 'gene_fwd'
        assert result[(100, 200, '-')] == 'gene_rev'
    
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
        assert result[(100, 200, '+')] == 'first'


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
# CHECK CDS TESTS
# =============================================================================

class TestCheckCds:
    """Tests for check_cds.py."""

    @staticmethod
    def _run(tmpdir, gene_id, rvs):
        """Annotate one 36 bp gene with no stop codon, with a TAA just past it.

        Whether the gene is treated as a known reference pseudogene decides
        between keeping 1-36 and extending to 1-39.
        """
        import subprocess
        genome = os.path.join(tmpdir, 'g.fasta')
        with open(genome, 'w') as f:
            f.write(">chr1\n" + "ATG" + "GCT" * 11 + "TAA" + "GCT" * 20 + "\n")
        tsv = os.path.join(tmpdir, 'sample_corrected.tsv')
        with open(tsv, 'w') as f:
            f.write(f"{gene_id}\t+\t1\t36\t{rvs}\n")
        script = os.path.join(os.path.dirname(__file__), '..', 'bin', 'check_cds.py')
        subprocess.run([sys.executable, script, genome, tsv],
                       cwd=tmpdir, check=True, capture_output=True)
        out = pd.read_csv(os.path.join(tmpdir, 'sample.AMAP.tsv'), sep='\t')
        return int(out.iloc[0]['start']), int(out.iloc[0]['end']), out.iloc[0]['pseudogene']

    def test_reference_pseudogene_kept_for_mtbc(self):
        import tempfile as tf
        with tf.TemporaryDirectory() as d:
            assert self._run(d, 'S1_Rv0001', 'Rv0001_1-36') == (1, 36, 'yes')

    def test_reference_pseudogene_kept_for_other_organisms(self):
        """The check was gated on the ID starting with 'Rv'.

        A gene called thrA behaves identically to one called Rv0001, so it must
        not be extended just because of its name.
        """
        import tempfile as tf
        with tf.TemporaryDirectory() as d:
            assert self._run(d, 'S1_thrA', 'thrA_1-36') == (1, 36, 'yes')

    def test_reference_pseudogene_kept_with_several_references(self):
        """rvs is comma separated; any of them matching the position counts."""
        import tempfile as tf
        with tf.TemporaryDirectory() as d:
            assert self._run(d, 'S1_multi', 'Rv0001_1-36, Rv0002_500-600') == (1, 36, 'yes')

    def test_gene_without_reference_coordinates_is_extended(self):
        """No embedded coordinates means nothing vouches for the boundary."""
        import tempfile as tf
        with tf.TemporaryDirectory() as d:
            start, end, _ = self._run(d, 'S1_00042', '')
            assert (start, end) == (1, 39)


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
