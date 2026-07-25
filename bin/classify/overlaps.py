"""
Overlap detection and resolution for case 4 (conflicting coordinates).

Builds IntervalTrees, finds single/multiple overlap groups, and resolves
conflicts by comparing Bakta vs Miniprot alignment lengths.

Coordinate convention: all coordinates are 1-based inclusive (GFF/GenBank style).
Length = End - Start + 1. IntervalTree uses half-open [Start, End+1).
"""

from typing import Dict, List, Tuple, Optional

import pandas as pd
from intervaltree import IntervalTree

# Type aliases
Coord = Tuple[int, int]           # (Start, End) pair
Coord3 = Tuple[int, int, str]     # (Start, End, Strand) triple
CoordDict = Dict[Coord3, str]     # Gene coordinate dictionary


def get_overlap_percentage(range1: Coord, range2: Coord) -> float:
    """Calculate overlap percentage relative to the larger range.
    
    Coordinates are 1-based inclusive: length = End - Start + 1.
    Overlap of [100,200] and [150,300] = min(200,300) - max(100,150) + 1 = 51 bp.
    """
    intersection = max(0, min(range1[1], range2[1]) - max(range1[0], range2[0]) + 1)
    max_length = max(range1[1] - range1[0] + 1, range2[1] - range2[0] + 1)
    return (intersection / max_length) * 100 if max_length > 0 else 0


def get_subset_overlap_percentage(range1, range2):
    """Calculate overlap percentage relative to the smaller range.
    
    Coordinates are 1-based inclusive: length = End - Start + 1.
    """
    intersection = max(0, min(range1[1], range2[1]) - max(range1[0], range2[0]) + 1)
    min_length = min(range1[1] - range1[0] + 1, range2[1] - range2[0] + 1)
    return (intersection / min_length) * 100 if min_length > 0 else 0


def _range_length(r):
    """Length of a 1-based inclusive coordinate range."""
    return r[1] - r[0] + 1


def build_overlap_tree(tree_df, overlap_df, tree_id_col, overlap_id_col):
    """Build an IntervalTree and find overlaps between two DataFrames."""
    tree = IntervalTree()
    for idx, row in tree_df.iterrows():
        if row['Start'] < row['End']:
            tree[row['Start']:row['End']+1] = (idx, row[tree_id_col])

    overlap_counts = {}
    for _, row in overlap_df.iterrows():
        for iv in tree[row['Start']:row['End']+1]:
            tree_idx, tree_id = iv.data
            if tree_idx not in overlap_counts:
                overlap_counts[tree_idx] = {'count': 1, 'Prot Id': tree_id,
                                            'Locus Tags': [row[overlap_id_col]]}
            else:
                overlap_counts[tree_idx]['count'] += 1
                overlap_counts[tree_idx]['Locus Tags'].append(row[overlap_id_col])
    return overlap_counts


def create_overlap_groups(filtered_overlaps, case4_df, miniprot_df):
    """Create a dict of overlapping gene groups for resolution.
    
    Note: `filtered_overlaps` keys may come from either miniprot_df or case4_df
    indices (depending on which was the tree vs query in build_overlap_tree).
    When the key is a case4 index, it may not exist in miniprot_df — in that
    case we use (0, 0) as placeholder coords (callers that use case4 indices
    only need 'Prot Id' and 'Locus Tags', never the coordinates).
    """
    overlaps_df = pd.DataFrame.from_dict(filtered_overlaps, orient='index')
    groups = {}
    for i, row in overlaps_df.iterrows():
        case4_rows = case4_df[case4_df['Locus Tag'].isin(row['Locus Tags'])]
        # Safely access miniprot_df coords — index may not exist when
        # filtered_overlaps was built with case4 as the tree source
        if i in miniprot_df.index:
            mp_coords = [(miniprot_df.at[i, 'Start'], miniprot_df.at[i, 'End'])]
        else:
            mp_coords = [(0, 0)]
        groups[i] = {
            'case4': [(r['Start'], r['End']) for _, r in case4_rows.iterrows()],
            'miniprot_df': mp_coords,
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
        bakta_len = value['case4'][0][1] - value['case4'][0][0] + 1
        mp_row = miniprot_df.loc[miniprot_df['Prot Id'] == value['Prot Id']].iloc[0]
        rv_len = mp_row['End'] - mp_row['Start'] + 1
        strand = mp_row.get('Strand', '+')

        locus_tag = value['Locus Tags'][0]
        n_overlaps = len(new_groups[locus_tag]['Locus Tags']) if locus_tag in new_groups else 0

        if n_overlaps > 1:
            new_coord, new_id = (*value['miniprot_df'][0], strand), value['Prot Id']
        elif rv_len > bakta_len:
            new_coord, new_id = (*value['miniprot_df'][0], strand), value['Prot Id']
        elif rv_len < bakta_len and abs(bakta_len - rv_len) > 550:
            new_coord, new_id = (*value['case4'][0], strand), value['Locus Tags'][0]
        else:
            new_coord, new_id = (*value['miniprot_df'][0], strand), value['Prot Id']

        _insert_with_dedup(combined_dict, new_coord, new_id)


def resolve_multiple_overlaps(multiple_groups, miniprot_df, combined_dict):
    """
    Resolve N-to-1 overlaps (multiple Bakta CDS vs one Miniprot hit).
    
    Pick the best overlapping Bakta CDS, then compare to Miniprot by length.
    """
    for key, value in multiple_groups.items():
        miniprot_range = value['miniprot_df'][0]
        prot_id = value['Prot Id']

        # Get strand from miniprot
        mp_row = miniprot_df.loc[miniprot_df['Prot Id'] == prot_id]
        strand = mp_row.iloc[0]['Strand'] if not mp_row.empty else '+'

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
            chosen_coord, chosen_id = (*miniprot_range, strand), prot_id
        else:
            chosen_coord, chosen_id = (*best[0], strand), best[1]

        _insert_with_dedup(combined_dict, chosen_coord, chosen_id, threshold=550)


def _is_reference_gene(gene_id):
    """Check if a gene ID looks like a reference gene (has embedded coordinates).
    
    Reference proteins from extract_ref_proteins.py always have the format
    GENE_ID_START-END (e.g., Rv0001_1-1524, GENE_001_1234-5678).
    This replaces the hardcoded "Rv" check to support organism-agnostic mode.
    """
    import re
    return bool(re.search(r'_\d+-\d+$', str(gene_id)))


def _insert_with_dedup(combined_dict, new_coord, new_id, threshold=550):
    """Insert a coordinate into the dict, resolving overlaps with existing entries.
    
    Coordinates are (Start, End, Strand) 3-tuples. Overlap checks only
    compare genes on the SAME strand — opposite-strand genes never conflict.
    """
    new_start, new_end = new_coord[0], new_coord[1]
    new_strand = new_coord[2] if len(new_coord) > 2 else '+'
    coords_to_add = {new_coord: new_id}
    coords_to_remove = set()

    for exist_coord, exist_id in list(combined_dict.items()):
        if new_coord == exist_coord:
            continue
        exist_start, exist_end = exist_coord[0], exist_coord[1]
        exist_strand = exist_coord[2] if len(exist_coord) > 2 else '+'
        # Only compare genes on the same strand
        if new_strand != exist_strand:
            continue
        if new_end <= exist_start or new_start >= exist_end:
            continue  # no overlap

        pct = get_overlap_percentage((new_start, new_end), (exist_start, exist_end))
        if pct >= 70:
            len_new = _range_length((new_start, new_end))
            len_exist = _range_length((exist_start, exist_end))
            if len_new > len_exist and (len_new - len_exist) > threshold:
                coords_to_remove.add(exist_coord)
            elif (len_exist - len_new) < threshold and not _is_reference_gene(exist_id):
                coords_to_remove.add(exist_coord)
            else:
                # The existing entry wins, so the new coordinate is discarded.
                # Any removals queued so far were only justified as replacements
                # for new_coord; applying them now would delete those genes
                # without inserting anything in their place.
                coords_to_add.clear()
                coords_to_remove.clear()
                break

    for c in coords_to_remove:
        combined_dict.pop(c, None)
    combined_dict.update(coords_to_add)
