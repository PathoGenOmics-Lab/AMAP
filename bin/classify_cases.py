#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AMAP CDS Classification — CLI entry point.

Classifies Bakta and Miniprot records into correction cases,
resolves overlaps, and produces a corrected gene coordinate set.
"""

import argparse
import os
import sys

# Allow importing from bin/classify/ when run as a script
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from classify.io import read_bakta, read_miniprot, read_miniprot_unfiltered, write_summary
from classify.cases import (
    classify_case2_case3, classify_case1_case4,
    find_rv_missing, filter_case1
)
from classify.overlaps import (
    build_overlap_tree, create_overlap_groups,
    resolve_single_overlaps, resolve_multiple_overlaps
)
from classify.reconciliation import (
    reconcile_coordinates, recover_missing_genes,
    deduplicate, write_corrected_tsv
)
from classify.plotting import generate_overlap_plot
from classify.ids import format_cluster_id

__author__ = 'Paula Ruiz-Rodriguez'
__credits__ = ['Paula Ruiz-Rodriguez', 'PathoGenOmics']
__license__ = ['GPL 3']
__version__ = '1.0.0'
__maintainer__ = '@paururo'
__email__ = 'paula.ruiz.rodriguez@csic.es'
__status__ = 'finished'
__lab__ = 'PathoGenOmics (I2SysBio)'


def create_coords_dict(df, id_col):
    """Create a coordinate→ID dict. First occurrence wins for duplicate coordinates."""
    result = {}
    for _, row in df.iterrows():
        key = (row['Start'], row['End'])
        if key not in result:
            result[key] = row[id_col]
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Classify Bakta and Miniprot records and resolve overlaps.")
    parser.add_argument('-b', '--bakta', required=True, help="Bakta TSV file.")
    parser.add_argument('-p', '--miniprot', required=True, help="Miniprot PAF file.")
    parser.add_argument('-bb', '--badbakta', required=True, help="Miniprot 'bad bakta' PAF.")
    parser.add_argument('-c', '--cds', required=True, help="Reference CDS PAF for Miniprot.")
    args = parser.parse_args()

    base_name = os.path.basename(args.bakta).split('.')[0]

    # --- 1. Read inputs ---
    bakta_df = read_bakta(args.bakta)
    miniprot_df = read_miniprot(args.miniprot)
    bakta_df_bad = read_miniprot_unfiltered(args.badbakta)

    # --- 2. Classify cases ---
    case2, case3 = classify_case2_case3(bakta_df, miniprot_df, base_name)
    case4, case1 = classify_case1_case4(bakta_df, miniprot_df, base_name)
    rv_missing = find_rv_missing(miniprot_df, bakta_df, base_name)
    case1_keep, case1_remove = filter_case1(case1, bakta_df_bad, args.cds)

    # --- 3. Write summary ---
    write_summary(base_name,
                  [case2, case3, case4, case1, rv_missing, bakta_df, case1_remove, case1_keep],
                  ["case2", "case3", "case4", "case1_initial", "rv_missing",
                   "total_cds", "case1_removed", "case1_kept"])

    # --- 4. Build overlap groups ---
    mp_case4_overlaps = build_overlap_tree(miniprot_df, case4, 'Prot Id', 'Locus Tag')
    case4_mp_overlaps = build_overlap_tree(case4, miniprot_df, 'Locus Tag', 'Prot Id')

    single_overlaps = {k: v for k, v in mp_case4_overlaps.items() if v['count'] == 1}
    multiple_overlaps = {k: v for k, v in mp_case4_overlaps.items() if v['count'] > 1}

    single_groups = create_overlap_groups(single_overlaps, case4, miniprot_df)
    multiple_groups = create_overlap_groups(multiple_overlaps, case4, miniprot_df)

    # --- 5. Generate overlap plots ---
    generate_overlap_plot(single_groups, f"{base_name}_single_overlaps.html")
    generate_overlap_plot(multiple_groups, f"{base_name}_multiple_overlaps.html")

    # --- 6. Build combined coordinate dict ---
    combined = {}
    combined.update(create_coords_dict(case2, 'Locus Tag'))
    combined.update(create_coords_dict(case3, 'Locus Tag'))
    combined.update(create_coords_dict(rv_missing, 'Prot Id'))
    combined.update(create_coords_dict(case1_keep, 'Locus Tag'))

    # --- 7. Resolve overlaps ---
    b_overlap_groups = create_overlap_groups(case4_mp_overlaps, case4, miniprot_df)
    new_groups = {v['Prot Id']: v for _, v in b_overlap_groups.items()}

    resolve_single_overlaps(single_groups, miniprot_df, combined, new_groups)
    resolve_multiple_overlaps(multiple_groups, miniprot_df, combined)

    # --- 8. Reconciliation + recovery + dedup ---
    reconcile_coordinates(combined, miniprot_df, bakta_df)
    recover_missing_genes(combined, miniprot_df)
    deduplicate(combined)

    # --- 9. Write final output ---
    write_corrected_tsv(combined, miniprot_df, bakta_df, base_name)


if __name__ == "__main__":
    main()
