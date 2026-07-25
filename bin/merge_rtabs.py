#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge multiple pangene Rtab files into a unified presence/absence matrix.

Each per-sample Rtab has potentially different gene families (rows) and
different reference columns. This script performs an outer join on the
Gene column, filling absent entries with 0.

Usage:
    merge_rtabs.py --output merged.Rtab file1.Rtab file2.Rtab ...
"""

import argparse
import sys
import pandas as pd
from pathlib import Path


def merge_rtab_files(rtab_paths, output_path):
    """Merge multiple Rtab files into a single presence/absence matrix.
    
    Each Rtab is a tab-separated file with:
    - First column: gene family name
    - Remaining columns: samples/references (0/1 values)
    
    The merge performs an outer join on gene family names, filling
    missing entries with 0 (absent).
    """
    if not rtab_paths:
        # Write empty Rtab with just a header placeholder
        with open(output_path, 'w') as f:
            f.write("Gene\n")
        return 0

    merged = None
    n_files = 0

    for rtab_path in rtab_paths:
        path = Path(rtab_path)
        if not path.exists() or path.stat().st_size == 0:
            continue

        try:
            df = pd.read_csv(path, sep='\t', index_col=0)
        except Exception as e:
            print(f"WARNING: Could not read {path}: {e}", file=sys.stderr)
            continue

        if df.empty:
            continue

        n_files += 1

        if merged is None:
            merged = df
        else:
            # Outer join on gene index, fill missing with 0
            merged = merged.join(df, how='outer', rsuffix=f'_dup{n_files}')

    if merged is None:
        with open(output_path, 'w') as f:
            f.write("Gene\n")
        return 0

    # Fill NaN (genes absent in some samples) with 0
    merged = merged.fillna(0).astype(int)

    # Write output — use index_label='Gene' so R reads it correctly
    # (without a named index column, R treats the first data column as row.names,
    #  breaking downstream pa_matrix = pangene_raw[, -1] logic)
    merged.to_csv(output_path, sep='\t', index_label='Gene')
    print(f"Merged {n_files} Rtab files: {merged.shape[0]} gene families × {merged.shape[1]} samples")
    return merged.shape[0]


def main():
    parser = argparse.ArgumentParser(
        description="Merge per-sample pangene Rtab files into a unified matrix.")
    parser.add_argument("rtab_files", nargs='*', help="Input Rtab files to merge.")
    parser.add_argument("--output", "-o", required=True, help="Output merged Rtab file.")
    args = parser.parse_args()

    merge_rtab_files(args.rtab_files, args.output)


if __name__ == "__main__":
    main()
