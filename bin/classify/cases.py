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

from .io import PAF_COLS


def _with_contig(df):
    """A frame built without contigs (tests, single-sequence callers) is one contig."""
    return df if 'Contig' in df.columns else df.assign(Contig='')


def _trees(df):
    """IntervalTrees of df's features, one per (contig, strand)."""
    trees = {}
    for row in df.itertuples(index=False):
        if row.Start < row.End:
            tree = trees.setdefault((row.Contig, '-' if row.Strand == '-' else '+'), IntervalTree())
            tree[row.Start:row.End + 1] = True
    return trees


def _overlaps_any(trees, contig, strand, start, end):
    tree = trees.get((contig, '-' if strand == '-' else '+'))
    return tree is not None and len(tree[start:end + 1]) > 0


def _overlapping(df, trees):
    """Whether each feature of df overlaps one in trees, as a boolean mask.

    Typed bool so that it stays a mask when df is empty: a column built from
    an empty list is float64, and indexing a frame with it selects columns by
    name. With every Bakta call in case 2 or 3, none is left for case 1 or
    case 4, and the split stopped with a KeyError.
    """
    return pd.Series([_overlaps_any(trees, r.Contig, r.Strand, r.Start, r.End)
                      for r in df.itertuples(index=False)], index=df.index, dtype=bool)


def classify_case2_case3(bakta_df, miniprot_df, base_name=None):
    """Split into case 2 (perfect match) and case 3 (coordinate match, length mismatch)."""
    # Merge on Contig, Start, End and Strand: same coordinates on another strand
    # or another contig are another gene.
    bakta_for_merge = _with_contig(bakta_df)
    mp_for_merge = _with_contig(miniprot_df)
    merged = pd.merge(bakta_for_merge, mp_for_merge, on=['Contig', 'Start', 'End', 'Strand'])
    merged['Difference'] = merged['Length_nt'] - merged['Num nt']

    case2 = merged[merged['Difference'] == 0]
    case3 = merged[merged['Difference'] != 0]

    case2 = case2[["Locus Tag", "Start", "End", "Strand", 'Length_nt', 'Contig']].drop_duplicates()

    if base_name:
        case3.to_csv(f"{base_name}-case3.tsv", sep='\t', index=False)
        case2.to_csv(f"{base_name}-case2.tsv", sep='\t', index=False)

    return case2, case3


def classify_case1_case4(bakta_df, miniprot_df, base_name=None):
    """Split non-matching Bakta CDS into case 1 (no overlap) and case 4 (overlap)."""
    bakta_df = _with_contig(bakta_df).copy()
    miniprot_df = _with_contig(miniprot_df)
    bakta_df['Length_nt'] = bakta_df['End'] - bakta_df['Start'] + 1
    # Merge on Contig, Start, End AND Strand to avoid cross-strand false exclusions.
    # Without Strand: a gene at (100,200,+) matching miniprot at (100,200,-)
    # would be excluded here AND from case2/3 → silently lost from annotation.
    merged = pd.merge(bakta_df, miniprot_df, on=['Contig', 'Start', 'End', 'Strand'], how='inner')

    common_cols = [c for c in bakta_df.columns if c in merged.columns]
    indicator = bakta_df.merge(merged, on=common_cols, how='left', indicator=True)
    not_in_merged = indicator[indicator['_merge'] == 'left_only'].drop('_merge', axis=1)
    cols = ["Locus Tag", "Start", "End", "Strand", 'Length_nt', 'Contig']
    not_in_merged = not_in_merged[cols]

    overlaps = _overlapping(not_in_merged, _trees(miniprot_df))
    case4 = not_in_merged[overlaps]
    case1 = not_in_merged[~overlaps]

    if base_name:
        case4.to_csv(f"{base_name}-case4.tsv", sep='\t', index=False)
        case1.to_csv(f"{base_name}-case1.tsv", sep='\t', index=False)

    return case4, case1


def find_rv_missing(miniprot_df, bakta_df, base_name=None):
    """Find Miniprot hits with no overlap to any Bakta CDS on the same contig and strand."""
    mp = _with_contig(miniprot_df)
    no_overlap = mp[~_overlapping(mp, _trees(_with_contig(bakta_df)))]
    no_overlap = no_overlap[~((no_overlap['Match Nt'] == 0) & (no_overlap['Num nt'] == 0))]

    if base_name:
        # The report reads this file by position: keep its seven columns.
        no_overlap[[c for c in PAF_COLS if c in no_overlap.columns]].to_csv(
            f"{base_name}-rvmissing.tsv", sep='\t', index=False)

    return no_overlap


def filter_case1(case1, bakta_df_bad, cds_path):
    """
    Refine case 1: separate into keep (truly novel) and remove (likely false positives).
    
    A gene is removed if:
    - ≥80% overlap with a known problematic annotation, OR
    - ≥80% match identity against the reference
    """
    from .io import read_miniprot_unfiltered

    # Build overlap trees from bad-bakta alignments, one per contig
    bad_trees = {}
    for row in _with_contig(bakta_df_bad).itertuples(index=False):
        if row.Start < row.End:
            bad_trees.setdefault(row.Contig, IntervalTree())[row.Start:row.End + 1] = True

    def calc_overlap(contig, start, end):
        total = 0
        for iv in bad_trees.get(contig, IntervalTree())[start:end+1]:
            # iv.begin = Start (1-based), iv.end = End+1 (half-open)
            # Use end+1 and iv.end directly for half-open interval math
            total += max(0, min(end + 1, iv.end) - max(start, iv.begin))
        return total

    case1 = _with_contig(case1).copy()
    case1['Overlap_Length'] = [calc_overlap(r.Contig, r.Start, r.End)
                               for r in case1.itertuples(index=False)]
    case1['Overlap_Percentage'] = (case1['Overlap_Length'] / case1['Length_nt']) * 100

    high_overlap = case1[case1['Overlap_Percentage'] >= 80]
    low_overlap = case1[case1['Overlap_Percentage'] < 80]

    # Check low-overlap genes against reference
    locus_ref = read_miniprot_unfiltered(cds_path)
    tags_to_check = low_overlap['Locus Tag'].tolist()
    matches = locus_ref[locus_ref['Prot Id'].isin(tags_to_check)].copy()
    matches['match_per'] = matches.apply(
        lambda r: (r['Match Nt'] / r['Num nt']) * 100 if r['Num nt'] > 0 else 0, axis=1
    )

    remove_ids = set(matches[matches['match_per'] >= 80]['Prot Id'])
    keep_ids = set(matches[matches['match_per'] < 80]['Prot Id'])

    case1_remove2 = low_overlap[low_overlap['Locus Tag'].isin(remove_ids)]
    case1_keep = low_overlap[low_overlap['Locus Tag'].isin(keep_ids)]

    # Genes not found in reference at all → keep (potentially novel)
    unchecked = low_overlap[~low_overlap['Locus Tag'].isin(remove_ids | keep_ids)]
    case1_keep = pd.concat([case1_keep, unchecked])

    case1_remove = pd.concat([high_overlap, case1_remove2])
    return case1_keep, case1_remove
