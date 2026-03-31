"""
Case classification: split genes into case1–4 + rv_missing.

- Case 1: Bakta CDS with no overlap to Miniprot (potentially novel)
- Case 2: Perfect coordinate match between Bakta and Miniprot
- Case 3: Same coordinates but different alignment length
- Case 4: Overlapping but mismatched coordinates
- rv_missing: Miniprot hits with no overlap to Bakta
"""

import pandas as pd
from intervaltree import IntervalTree


def classify_case2_case3(bakta_df, miniprot_df, base_name=None):
    """Split into case 2 (perfect match) and case 3 (coordinate match, length mismatch)."""
    merged = pd.merge(bakta_df, miniprot_df, on=['Start', 'End'])
    merged['Difference'] = merged['Length_nt'] - merged['Num nt']

    case2 = merged[merged['Difference'] == 0]
    case3 = merged[merged['Difference'] != 0]

    case2 = case2[["Locus Tag", "Start", "End", "Strand_x", 'Length_nt']].drop_duplicates()

    if base_name:
        case3.to_csv(f"{base_name}-case3.tsv", sep='\t', index=False)
        case2.to_csv(f"{base_name}-case2.tsv", sep='\t', index=False)

    return case2, case3


def classify_case1_case4(bakta_df, miniprot_df, base_name=None):
    """Split non-matching Bakta CDS into case 1 (no overlap) and case 4 (overlap)."""
    bakta_df = bakta_df.copy()
    bakta_df['Length_nt'] = bakta_df['End'] - bakta_df['Start']
    merged = pd.merge(bakta_df, miniprot_df, on=['Start', 'End'], how='inner')

    common_cols = [c for c in bakta_df.columns if c in merged.columns]
    indicator = bakta_df.merge(merged, on=common_cols, how='left', indicator=True)
    not_in_merged = indicator[indicator['_merge'] == 'left_only'].drop('_merge', axis=1)
    not_in_merged = not_in_merged[["Locus Tag", "Start", "End", "Strand", 'Length_nt']]

    # Build IntervalTree for fast overlap queries
    tree = IntervalTree()
    for _, row in miniprot_df.iterrows():
        if row['Start'] < row['End']:
            tree[row['Start']:row['End']] = True

    not_in_merged = not_in_merged.copy()
    not_in_merged['Overlaps'] = not_in_merged.apply(
        lambda r: len(tree[r['Start']:r['End']]) > 0, axis=1)

    case4 = not_in_merged[not_in_merged['Overlaps']][["Locus Tag", "Start", "End", "Strand", 'Length_nt']]
    case1 = not_in_merged[~not_in_merged['Overlaps']][["Locus Tag", "Start", "End", "Strand", 'Length_nt']]

    if base_name:
        case4.to_csv(f"{base_name}-case4.tsv", sep='\t', index=False)
        case1.to_csv(f"{base_name}-case1.tsv", sep='\t', index=False)

    return case4, case1


def find_rv_missing(miniprot_df, bakta_df, base_name=None):
    """Find Miniprot hits with no overlap to any Bakta CDS (missing reference genes)."""
    tree = IntervalTree()
    for _, row in bakta_df.iterrows():
        if row['Start'] < row['End']:
            tree[row['Start']:row['End']] = True

    mp = miniprot_df.copy()
    mp['Overlaps'] = mp.apply(lambda r: len(tree[r['Start']:r['End']]) > 0, axis=1)

    no_overlap = mp[~mp['Overlaps']]
    no_overlap = no_overlap[~((no_overlap['Match Nt'] == 0) & (no_overlap['Num nt'] == 0))]

    if base_name:
        no_overlap.to_csv(f"{base_name}-rvmissing.tsv", sep='\t', index=False)

    return no_overlap


def filter_case1(case1, bakta_df_bad, cds_path):
    """
    Refine case 1: separate into keep (truly novel) and remove (likely false positives).
    
    A gene is removed if:
    - ≥80% overlap with a known problematic annotation, OR
    - ≥80% match identity against the reference
    """
    from .io import read_miniprot_unfiltered

    # Build overlap tree from bad-bakta alignments
    bad_tree = IntervalTree()
    for _, row in bakta_df_bad.iterrows():
        if row['Start'] < row['End']:
            bad_tree[row['Start']:row['End']] = True

    def calc_overlap(start, end):
        total = 0
        for iv in bad_tree[start:end]:
            total += max(0, min(end, iv.end) - max(start, iv.begin))
        return total

    case1 = case1.copy()
    case1['Overlap_Length'] = case1.apply(lambda r: calc_overlap(r['Start'], r['End']), axis=1)
    case1['Overlap_Percentage'] = (case1['Overlap_Length'] / case1['Length_nt']) * 100

    high_overlap = case1[case1['Overlap_Percentage'] >= 80]
    low_overlap = case1[case1['Overlap_Percentage'] < 80]

    # Check low-overlap genes against reference
    locus_ref = read_miniprot_unfiltered(cds_path)
    tags_to_check = low_overlap['Locus Tag'].tolist()
    matches = locus_ref[locus_ref['Prot Id'].isin(tags_to_check)].copy()
    matches['match_per'] = 100 - ((matches['Num nt'] - matches['Match Nt']) / matches['Num nt']) * 100
    matches['match_per'] = matches['match_per'].fillna(0)

    remove_ids = set(matches[matches['match_per'] >= 80]['Prot Id'])
    keep_ids = set(matches[matches['match_per'] < 80]['Prot Id'])

    case1_remove2 = low_overlap[low_overlap['Locus Tag'].isin(remove_ids)]
    case1_keep = low_overlap[low_overlap['Locus Tag'].isin(keep_ids)]

    # Genes not found in reference at all → keep (potentially novel)
    unchecked = low_overlap[~low_overlap['Locus Tag'].isin(remove_ids | keep_ids)]
    case1_keep = pd.concat([case1_keep, unchecked])

    case1_remove = pd.concat([high_overlap, case1_remove2])
    return case1_keep, case1_remove
