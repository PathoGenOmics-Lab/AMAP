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
import json
import sys
import os
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


# EGGNOG_DB (modules/local/eggnog.nf) takes a smaller bacteria.dmnd for a
# truncated download and builds it again. The two must agree.
EGGNOG_DMND_MIN_BYTES = 100_000_000


def validate_bakta_db(path, shown_as=None):
    """Check a Bakta database directory. Returns (errors, warnings, note),
    the note saying which database it is.

    `shown_as` is the name messages give it, when it is read from somewhere
    other than the path the user gave."""
    name = shown_as or path
    if not os.path.exists(path):
        return [f"Bakta DB not found: {name}"], [], None
    if not os.path.isdir(path):
        return [f"Bakta DB is not a directory: {name}"], [], None

    version_file = os.path.join(path, 'version.json')
    if not os.path.isfile(version_file):
        # `bakta_db download --output X` puts the database in X/db or X/db-light
        inner = [d for d in ('db', 'db-light')
                 if os.path.isfile(os.path.join(path, d, 'version.json'))]
        if inner:
            return [f"Bakta DB: {name} is the folder the database was downloaded to; "
                    f"point --bakta_db at {os.path.join(name, inner[0])}"], [], None
        return [f"Bakta DB: {name} has no version.json, so it is not a Bakta database"], [], None
    try:
        with open(version_file) as f:
            version = json.load(f)
        note = f"{version['major']}.{version['minor']}"
    except (OSError, ValueError, KeyError, TypeError) as e:
        return [f"Bakta DB: cannot read the version.json in {name}: {e}"], [], None
    if version.get('type'):
        note = f"{version['type']} {note}"

    if not os.path.isfile(os.path.join(path, 'bakta.db')):
        return [f"Bakta DB: {name} has no bakta.db"], [], note
    return [], [], note


def validate_eggnog_db(path, shown_as=None):
    """Check eggNOG's database as EGGNOG_DB and eggNOG-mapper will find it.
    Returns (errors, warnings, note).

    It need not exist yet: EGGNOG_DB downloads and builds it where there is no
    bacteria.dmnd. Once eggnog_db.ok is in it, though, Nextflow takes it for
    built and never runs EGGNOG_DB again (its storeDir), so whatever is
    missing then stays missing."""
    name = shown_as or path
    if os.path.exists(path) and not os.path.isdir(path):
        return [f"eggNOG DB is not a directory: {name}"], [], None

    dmnd = os.path.join(path, 'bacteria.dmnd')
    dmnd_bytes = os.path.getsize(dmnd) if os.path.isfile(dmnd) else None
    dmnd_ok = dmnd_bytes is not None and dmnd_bytes >= EGGNOG_DMND_MIN_BYTES
    # eggNOG-mapper's annotations, which it reads from --data_dir
    annotations_ok = os.path.isfile(os.path.join(path, 'eggnog.db'))
    if dmnd_ok and annotations_ok:
        return [], [], None

    if os.path.isfile(os.path.join(path, 'eggnog_db.ok')):
        if dmnd_bytes is None:
            missing = "there is no bacteria.dmnd"
        elif not dmnd_ok:
            missing = f"bacteria.dmnd is only {dmnd_bytes:,} bytes"
        else:
            missing = "there is no eggnog.db"
        # Built again from scratch only when bacteria.dmnd is missing or short
        delete = "eggnog_db.ok" if not dmnd_ok else "eggnog_db.ok and bacteria.dmnd"
        return [f"eggNOG DB: {name} has eggnog_db.ok, so EGGNOG_DB will not run again, "
                f"but {missing}; delete {delete} to have it built again"], [], None
    if dmnd_bytes is None:
        return [], [], "not built yet: EGGNOG_DB will download and build it, several gigabytes"
    if not dmnd_ok:
        return [], [f"eggNOG DB: bacteria.dmnd in {name} is only {dmnd_bytes:,} bytes; "
                    f"EGGNOG_DB will take it for a truncated download and build it again"], None
    return [f"eggNOG DB: {name} has no eggnog.db, which eggNOG-mapper needs, and EGGNOG_DB "
            f"will not download it while bacteria.dmnd is there; delete bacteria.dmnd to "
            f"have it built again"], [], None


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
    parser.add_argument("--bakta-db-staged", default=None,
                        help="Where to read the Bakta database instead, when the "
                             "pipeline staged it; messages keep naming --bakta-db.")
    parser.add_argument("--eggnog-db", default=None,
                        help="eggNOG database path. It need not exist yet: "
                             "EGGNOG_DB builds it.")
    parser.add_argument("--eggnog-db-staged", default=None,
                        help="Where to read the eggNOG database instead, when the "
                             "pipeline staged it; messages keep naming --eggnog-db.")
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
    databases = [(args.bakta_db, args.bakta_db_staged, 'Bakta DB', validate_bakta_db),
                 (args.eggnog_db, args.eggnog_db_staged, 'eggNOG DB', validate_eggnog_db)]
    for db_path, db_staged, db_name, check in databases:
        if not db_path:
            print(f"  - {db_name}: not specified")
            continue
        db_err, db_warn, note = check(db_staged or db_path, shown_as=db_path)
        all_errors.extend(db_err)
        all_warnings.extend(db_warn)
        status = "✓" if not db_err else "✗"
        print(f"  {status} {db_name}: {db_path}" + (f" ({note})" if note else ""))
        for w in db_warn:
            print(f"    ⚠ {w}")
        for e in db_err:
            print(f"    ✗ {e}")

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
        print(f"  Fix the warnings, or re-run without --strict_validation")
        print(f"  (--strict when running this script yourself).")
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
