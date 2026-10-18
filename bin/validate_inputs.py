#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AMAP Input Validation — pre-flight checks for FASTA files and pipeline inputs.

Validates:
1. FASTA format (parseable, non-empty sequences)
2. Sequence quality (N-content, minimum length, GC%)
3. Assembly metrics (total size, number of contigs, N50)
4. Potential contamination flags (abnormal size, extreme GC)
5. Database paths (Bakta DB, eggNOG DB)

Usage:
    validate_inputs.py --tsv samples.tsv [--genomes staged.tsv] [--bakta-db /path] [--eggnog-db /path] [--strict]

The sample sheet says where each genome is on the host, absolute or relative
to where Nextflow was launched. A container sees only the folders Nextflow
mounts for it, and a task runs in a directory of its own, so in the pipeline
VALIDATE_INPUTS stages every genome as an input and passes --genomes, a table
of where each sample's genome was staged, and each one is read from there.
Messages still name the sheet's path.

Exit codes:
    0 = all checks passed
    1 = warnings only (non-strict mode)
    2 = critical errors found
"""

import argparse
import csv
import sys
import os
from pathlib import Path
from collections import Counter

# Avoid requiring Biopython for validation — parse FASTA manually
def parse_fasta(path):
    """Minimal FASTA parser. Yields (header, sequence) tuples."""
    header = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_parts)
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
    if header is not None:
        yield header, ''.join(seq_parts)


def calc_n50(lengths):
    """Calculate N50 from a list of contig lengths."""
    if not lengths:
        return 0
    sorted_lens = sorted(lengths, reverse=True)
    total = sum(sorted_lens)
    running = 0
    for l in sorted_lens:
        running += l
        if running >= total / 2:
            return l
    return 0


def validate_fasta(path, sample_id, strict=False, shown_as=None):
    """Validate a single FASTA file. Returns (errors, warnings).

    `shown_as` is the name messages give the file, when it is read from
    somewhere other than the path the user knows it by."""
    errors = []
    warnings = []
    name = shown_as or path

    # Check file exists and is readable
    if not os.path.isfile(path):
        return [f"FASTA file not found: {name}"], []
    
    file_size = os.path.getsize(path)
    if file_size == 0:
        return [f"FASTA file is empty: {name}"], []

    # Parse sequences
    try:
        contigs = list(parse_fasta(path))
    except Exception as e:
        return [f"Cannot parse FASTA {name}: {e}"], []

    if len(contigs) == 0:
        return [f"No sequences found in {name}"], []

    # Basic metrics
    n_contigs = len(contigs)
    lengths = [len(seq) for _, seq in contigs]
    total_size = sum(lengths)
    n50 = calc_n50(lengths)
    max_len = max(lengths)
    min_len = min(lengths)

    # Sequence composition
    all_bases = ''.join(seq for _, seq in contigs)
    base_counts = Counter(all_bases)
    n_count = base_counts.get('N', 0)
    n_pct = (n_count / total_size * 100) if total_size > 0 else 0
    gc_count = base_counts.get('G', 0) + base_counts.get('C', 0)
    atgc_count = sum(base_counts.get(b, 0) for b in 'ATGC')
    gc_pct = (gc_count / atgc_count * 100) if atgc_count > 0 else 0

    # Non-ATGCN characters
    valid_bases = set('ATGCNRYSWKMBDHV')
    invalid = set(all_bases) - valid_bases
    if invalid:
        warnings.append(f"[{sample_id}] Non-standard characters in sequence: {invalid}")

    # --- CRITICAL CHECKS ---
    
    # Minimum genome size (prokaryote: ~500 kb minimum for Mycoplasma, ~4 Mb for MTB)
    if total_size < 100_000:
        errors.append(f"[{sample_id}] Genome too small: {total_size:,} bp (min 100 kb)")
    elif total_size < 500_000:
        warnings.append(f"[{sample_id}] Very small genome: {total_size:,} bp (expected >500 kb for most bacteria)")

    # Maximum genome size (prokaryote: ~15 Mb for Streptomyces, ~4.4 Mb for MTB)
    if total_size > 20_000_000:
        errors.append(f"[{sample_id}] Genome too large: {total_size:,} bp (max 20 Mb for prokaryote)")
    elif total_size > 10_000_000:
        warnings.append(f"[{sample_id}] Large genome: {total_size:,} bp (unusual for MTBC)")

    # N-content
    if n_pct > 10:
        errors.append(f"[{sample_id}] Excessive N-content: {n_pct:.1f}% ({n_count:,} Ns)")
    elif n_pct > 1:
        warnings.append(f"[{sample_id}] Elevated N-content: {n_pct:.1f}% ({n_count:,} Ns)")

    # GC content sanity (MTBC is ~65%, but pipeline supports any prokaryote)
    if gc_pct < 20 or gc_pct > 80:
        errors.append(f"[{sample_id}] Extreme GC%: {gc_pct:.1f}% (expected 20-80% for bacteria)")
    elif gc_pct < 25 or gc_pct > 75:
        warnings.append(f"[{sample_id}] Unusual GC%: {gc_pct:.1f}%")

    # Empty sequences
    empty_seqs = [h for h, s in contigs if len(s) == 0]
    if empty_seqs:
        errors.append(f"[{sample_id}] Empty sequences: {', '.join(empty_seqs[:5])}")

    # Short contigs
    short = sum(1 for l in lengths if l < 200)
    if short > 0 and n_contigs > 1:
        warnings.append(f"[{sample_id}] {short}/{n_contigs} contigs shorter than 200 bp")

    # --- ASSEMBLY QUALITY WARNINGS ---
    
    # Fragmented assembly
    if n_contigs > 500:
        warnings.append(f"[{sample_id}] Highly fragmented: {n_contigs} contigs (N50={n50:,} bp)")
    elif n_contigs > 100:
        warnings.append(f"[{sample_id}] Fragmented assembly: {n_contigs} contigs")

    # Duplicate headers
    headers = [h for h, _ in contigs]
    if len(headers) != len(set(headers)):
        dupes = [h for h in headers if headers.count(h) > 1]
        errors.append(f"[{sample_id}] Duplicate contig headers: {set(dupes)}")

    return errors, warnings


def validate_database(db_path, db_name):
    """Validate that a database directory exists and has expected contents."""
    errors = []
    warnings = []

    if not db_path:
        return errors, warnings  # Optional DB

    db_path = Path(db_path)
    
    if not db_path.exists():
        errors.append(f"{db_name} path does not exist: {db_path}")
        return errors, warnings

    if not db_path.is_dir():
        errors.append(f"{db_name} path is not a directory: {db_path}")
        return errors, warnings

    if db_name == 'Bakta DB':
        # Bakta DB should have specific files
        expected = ['db.sqlite3', 'version.json']  # Bakta light DB
        found = os.listdir(db_path)
        for exp in expected:
            # Check recursively (Bakta DB can have nested structure)
            if not any(exp in f for f in found):
                warnings.append(f"{db_name}: expected file '{exp}' not found in {db_path}")

    elif db_name == 'eggNOG DB':
        # eggNOG DB should have bacteria.dmnd at minimum
        expected = ['bacteria.dmnd']
        found = os.listdir(db_path)
        for exp in expected:
            if exp not in found:
                warnings.append(f"{db_name}: '{exp}' not found in {db_path} (will be downloaded)")

    return errors, warnings


def validate_tsv(tsv_path):
    """Validate the input TSV file format."""
    errors = []
    
    if not os.path.isfile(tsv_path):
        return [f"Input TSV not found: {tsv_path}"], []

    try:
        with open(tsv_path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            fields = reader.fieldnames

            if not fields:
                return [f"Empty TSV file: {tsv_path}"], []
            
            if 'sampleId' not in fields:
                errors.append(f"TSV missing 'sampleId' column. Found: {fields}")
            if 'fasta' not in fields:
                errors.append(f"TSV missing 'fasta' column. Found: {fields}")

            rows = list(reader)
            if len(rows) == 0:
                errors.append(f"TSV has headers but no data rows")
            
            # Check for duplicate sample IDs
            sample_ids = [r.get('sampleId', '') for r in rows]
            dupes = [s for s in sample_ids if sample_ids.count(s) > 1]
            if dupes:
                errors.append(f"Duplicate sampleIds in TSV: {set(dupes)}")

            # Check for empty sample IDs
            empty_ids = [i for i, s in enumerate(sample_ids, 2) if not s.strip()]
            if empty_ids:
                errors.append(f"Empty sampleId on line(s): {empty_ids}")

            return errors, rows

    except Exception as e:
        return [f"Cannot parse TSV {tsv_path}: {e}"], []


def read_genome_table(handle):
    """{sampleId: path} from a table with the columns sampleId and fasta.

    VALIDATE_INPUTS writes it with the path each genome was staged at. Nothing
    in it is quoted: a sample ID or a file name may hold a quote character."""
    reader = csv.DictReader(handle, delimiter='\t', quoting=csv.QUOTE_NONE)
    return {row['sampleId']: row['fasta'] for row in reader}


def main():
    parser = argparse.ArgumentParser(
        description="AMAP pre-flight input validation.")
    parser.add_argument("--tsv", required=True, help="Input samples TSV.")
    parser.add_argument("--genomes", type=argparse.FileType('r'), default=None,
                        help="Where to read each sample's genome instead of the "
                             "path in the TSV: a table with the columns sampleId "
                             "and fasta, or '-' to read it from standard input.")
    parser.add_argument("--bakta-db", default=None, help="Bakta database path.")
    parser.add_argument("--eggnog-db", default=None, help="eggNOG database path.")
    parser.add_argument("--strict", action="store_true",
                        help="Treat warnings as errors (fail on any issue).")
    parser.add_argument("--quiet", action="store_true",
                        help="Only print errors, not warnings.")
    args = parser.parse_args()

    all_errors = []
    all_warnings = []

    # 1. Validate TSV
    print("=" * 60)
    print("AMAP Input Validation")
    print("=" * 60)

    tsv_errors, rows = validate_tsv(args.tsv)
    all_errors.extend(tsv_errors)
    
    if tsv_errors:
        for e in tsv_errors:
            print(f"  ✗ {e}")
        print(f"\n{'=' * 60}")
        print(f"FAILED: {len(all_errors)} critical error(s)")
        sys.exit(2)

    n_samples = len(rows)
    print(f"\n  Samples: {n_samples}")
    print(f"  TSV: {args.tsv}")

    # 2. Validate databases
    print(f"\n--- Database Checks ---")
    for db_path, db_name in [(args.bakta_db, 'Bakta DB'), (args.eggnog_db, 'eggNOG DB')]:
        db_err, db_warn = validate_database(db_path, db_name)
        all_errors.extend(db_err)
        all_warnings.extend(db_warn)
        if db_path:
            status = "✓" if not db_err else "✗"
            print(f"  {status} {db_name}: {db_path}")
            for w in db_warn:
                print(f"    ⚠ {w}")
            for e in db_err:
                print(f"    ✗ {e}")
        else:
            print(f"  - {db_name}: not specified")

    # 3. Validate each FASTA
    print(f"\n--- FASTA Validation ({n_samples} samples) ---")
    
    n_ok = 0
    n_warn = 0
    n_fail = 0
    staged = read_genome_table(args.genomes) if args.genomes else {}

    for row in rows:
        sample_id = row.get('sampleId', 'unknown')
        fasta_path = row.get('fasta', '')

        errors, warnings = validate_fasta(staged.get(sample_id, fasta_path), sample_id,
                                          args.strict, shown_as=fasta_path)
        all_errors.extend(errors)
        all_warnings.extend(warnings)

        if errors:
            n_fail += 1
            for e in errors:
                print(f"  ✗ {e}")
        elif warnings:
            n_warn += 1
            if not args.quiet:
                for w in warnings:
                    print(f"  ⚠ {w}")
        else:
            n_ok += 1

    # 4. Summary
    print(f"\n{'=' * 60}")
    print(f"Results: {n_ok} OK, {n_warn} warnings, {n_fail} failed")
    print(f"Total: {len(all_errors)} error(s), {len(all_warnings)} warning(s)")
    print(f"{'=' * 60}")

    if all_errors:
        print(f"\n{'=' * 60}")
        print(f"  ✗ FAILED: {len(all_errors)} critical error(s)")
        print(f"")
        print(f"  Fix the issues above and re-run validation:")
        print(f"    nextflow run PathoGenOmics-Lab/AMAP --validate_only --tsv {args.tsv}")
        print(f"{'=' * 60}")
        sys.exit(2)
    elif all_warnings and args.strict:
        print(f"\n  ✗ FAILED (strict mode): {len(all_warnings)} warning(s)")
        print(f"  Re-run without --strict to proceed, or fix the warnings.")
        sys.exit(1)
    elif all_warnings:
        print(f"\n  ✓ PASSED with {len(all_warnings)} warning(s)")
        print(f"  Warnings are non-blocking but worth reviewing.")
        sys.exit(0)
    else:
        print(f"\n  ✓ PASSED — all {n_samples} samples look good")
        sys.exit(0)


if __name__ == "__main__":
    main()
