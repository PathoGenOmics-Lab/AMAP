"""
Overlap detection and resolution for case 4 (conflicting coordinates).

Builds IntervalTrees, finds single/multiple overlap groups, and resolves
conflicts by comparing Bakta vs Miniprot alignment lengths.
"""

import pandas as pd
from intervaltree import IntervalTree


def get_overlap_percentage(range1, range2):
    """Calculate overlap percentage relative to the larger range."""
    intersection = max(0, min(range1[1], range2[1]) - max(range1[0], range2[0]))
    max_length = max(range1[1] - range1[0], range2[1] - range2[0])
    return (intersection / max_length) * 100 if max_length > 0 else 0


def get_subset_overlap_percentage(range1, range2):
    """Calculate overlap percentage relative to the smaller range."""
    intersection = max(0, min(range1[1], range2[1]) - max(range1[0], range2[0]))
    min_length = min(range1[1] - range1[0], range2[1] - range2[0])
    return (intersection / min_length) * 100 if min_length > 0 else 0


def _range_length(r):
    return r[1] - r[0]


def build_overlap_tree(tree_df, overlap_df, tree_id_col, overlap_id_col):
    """Build an IntervalTree and find overlaps between two DataFrames."""
    tree = IntervalTree()
    for idx, row in tree_df.iterrows():
        if row['Start'] < row['End']:
            tree[row['Start']:row['End']] = (idx, row[tree_id_col])

    overlap_counts = {}
    for _, row in overlap_df.iterrows():
        for iv in tree[row['Start']:row['End']]:
            tree_idx, tree_id = iv.data
            if tree_idx not in overlap_counts:
                overlap_counts[tree_idx] = {'count': 1, 'Prot Id': tree_id,
                                            'Locus Tags': [row[overlap_id_col]]}
            else:
                overlap_counts[tree_idx]['count'] += 1
                overlap_counts[tree_idx]['Locus Tags'].append(row[overlap_id_col])
    return overlap_counts


def create_overlap_groups(filtered_overlaps, case4_df, miniprot_df):
    """Create a dict of overlapping gene groups for resolution."""
    overlaps_df = pd.DataFrame.from_dict(filtered_overlaps, orient='index')
    groups = {}
    for i, row in overlaps_df.iterrows():
        case4_rows = case4_df[case4_df['Locus Tag'].isin(row['Locus Tags'])]
        groups[i] = {
            'case4': [(r['Start'], r['End']) for _, r in case4_rows.iterrows()],
            'miniprot_df': [(miniprot_df.at[i, 'Start'], miniprot_df.at[i, 'End'])],
            'Prot Id': row['Prot Id'],
            'Locus Tags': row['Locus Tags']
        }
    return groups


def resolve_single_overlaps(single_groups, miniprot_df, combined_dict, new_groups):
    """
    Resolve 1-to-1 overlaps between Bakta and Miniprot.
    
    Rules:
    - Prefer Miniprot (reference) by default
    - Use Bakta only if substantially longer (>550 bp difference)
    - Remove existing entries with ≥70% overlap if the new one is longer
    """
    for key, value in single_groups.items():
        bakta_len = value['case4'][0][1] - value['case4'][0][0]
        # Use genomic span (End - Start) for fair comparison, not Match Nt
        # which only counts aligned bases and is always <= span
        mp_row = miniprot_df.loc[miniprot_df['Prot Id'] == value['Prot Id']].iloc[0]
        rv_len = mp_row['End'] - mp_row['Start']

        locus_tag = value['Locus Tags'][0]
        n_overlaps = len(new_groups[locus_tag]['Locus Tags']) if locus_tag in new_groups else 0

        if n_overlaps > 1:
            new_coord, new_id = value['miniprot_df'][0], value['Prot Id']
        elif rv_len > bakta_len:
            new_coord, new_id = value['miniprot_df'][0], value['Prot Id']
        elif rv_len < bakta_len and abs(bakta_len - rv_len) > 550:
            new_coord, new_id = value['case4'][0], value['Locus Tags'][0]
        else:
            new_coord, new_id = value['miniprot_df'][0], value['Prot Id']

        _insert_with_dedup(combined_dict, new_coord, new_id)


def resolve_multiple_overlaps(multiple_groups, miniprot_df, combined_dict):
    """
    Resolve N-to-1 overlaps (multiple Bakta CDS vs one Miniprot hit).
    
    Pick the best overlapping Bakta CDS, then compare to Miniprot by length.
    """
    for key, value in multiple_groups.items():
        miniprot_range = value['miniprot_df'][0]
        prot_id = value['Prot Id']

        best = (None, None, -1)
        for i, r in enumerate(value['case4']):
            pct = get_overlap_percentage(r, miniprot_range)
            if pct > best[2]:
                tag = value['Locus Tags'][i] if i < len(value['Locus Tags']) else "Unknown"
                best = (r, tag, pct)

        if not best[0]:
            continue

        if _range_length(miniprot_range) >= _range_length(best[0]) or \
           (_range_length(best[0]) - _range_length(miniprot_range)) < 500:
            chosen_coord, chosen_id = miniprot_range, prot_id
        else:
            chosen_coord, chosen_id = best[0], best[1]

        _insert_with_dedup(combined_dict, chosen_coord, chosen_id, threshold=550)


def _insert_with_dedup(combined_dict, new_coord, new_id, threshold=550):
    """Insert a coordinate into the dict, resolving overlaps with existing entries."""
    coords_to_add = {new_coord: new_id}
    coords_to_remove = set()

    for exist_coord, exist_id in list(combined_dict.items()):
        if new_coord == exist_coord:
            continue
        if new_coord[1] <= exist_coord[0] or new_coord[0] >= exist_coord[1]:
            continue  # no overlap

        pct = get_overlap_percentage(new_coord, exist_coord)
        if pct >= 70:
            len_new = _range_length(new_coord)
            len_exist = _range_length(exist_coord)
            if len_new > len_exist and (len_new - len_exist) > threshold:
                coords_to_remove.add(exist_coord)
            elif (len_exist - len_new) < threshold and "Rv" not in str(exist_id):
                coords_to_remove.add(exist_coord)
            else:
                coords_to_add.clear()
                break

    for c in coords_to_remove:
        combined_dict.pop(c, None)
    combined_dict.update(coords_to_add)
