"""
Overlap detection and resolution for case 4 (conflicting coordinates).

Builds IntervalTrees, finds single/multiple overlap groups, and resolves
conflicts by comparing Bakta vs Miniprot alignment lengths.

Coordinate convention: all coordinates are 1-based inclusive (GFF/GenBank style).
Length = End - Start + 1. IntervalTree uses half-open [Start, End+1).
A gene's key is (Start, End, Strand, Contig): the contig goes last so that code
reading the first three keeps working, and features only ever overlap when they
share a contig.
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


def _contig(coord):
    return coord[3] if len(coord) > 3 else ''


def build_overlap_tree(tree_df, overlap_df, tree_id_col, overlap_id_col):
    """For each feature of tree_df, the features of overlap_df overlapping it on
    the same strand of the same contig.

    Case 4 is defined strand by strand, but this grouping was not: about 8% of
    the groups (101 of 1255 on H37Rv) put an alignment on one strand together
    with a Bakta call on the other, and a gene on the other strand counted as a
    competing call.
    """
    def where(row):
        return row.get('Contig', ''), '-' if row.get('Strand', '+') == '-' else '+'

    trees = {}
    for idx, row in tree_df.iterrows():
        if row['Start'] < row['End']:
            tree = trees.setdefault(where(row), IntervalTree())
            tree[row['Start']:row['End']+1] = (idx, row[tree_id_col])

    overlap_counts = {}
    for _, row in overlap_df.iterrows():
        tree = trees.get(where(row))
        if tree is None:
            continue
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
            mp_contig = miniprot_df.at[i, 'Contig'] if 'Contig' in miniprot_df.columns else ''
        else:
            mp_coords = [(0, 0)]
            mp_contig = ''
        groups[i] = {
            'case4': [(r['Start'], r['End']) for _, r in case4_rows.iterrows()],
            # Kept apart from 'case4' because callers unpack those as pairs.
            'case4_strands': [r.get('Strand', '+') for _, r in case4_rows.iterrows()],
            'case4_contigs': [r.get('Contig', '') for _, r in case4_rows.iterrows()],
            'miniprot_df': mp_coords,
            'contig': mp_contig,
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
    reference_ids = set(miniprot_df['Prot Id'])
    for key, value in single_groups.items():
        bakta_len = value['case4'][0][1] - value['case4'][0][0] + 1
        # The group key is the index of the Miniprot hit being resolved. Looking
        # the protein up by name instead returns its first placement, which for
        # a protein with several copies can sit elsewhere on the other strand:
        # esxM (Rv1792) came out at its own coordinates on the wrong strand.
        mp_row = miniprot_df.loc[key]
        rv_len = mp_row['End'] - mp_row['Start'] + 1
        strand = mp_row.get('Strand', '+')
        contig = mp_row.get('Contig', '')
        bakta_strand = value['case4_strands'][0] if value.get('case4_strands') else strand
        bakta_contig = value['case4_contigs'][0] if value.get('case4_contigs') else contig

        locus_tag = value['Locus Tags'][0]
        n_overlaps = len(new_groups[locus_tag]['Locus Tags']) if locus_tag in new_groups else 0

        # A much longer Bakta call is taken as a sign the alignment fragmented,
        # but only an alignment that stops short of the end of its protein can
        # have. One that reaches it, stop codon included, is the whole gene,
        # and a longer call on it ends at the same stop and merely starts
        # further upstream: plcD (Rv1755c) on CDC1551, 843 bp in H37Rv, 1545 bp
        # in Bakta's call.
        complete = bool(mp_row.get('Complete', False))

        if n_overlaps > 1:
            new_coord, new_id = (*value['miniprot_df'][0], strand, contig), value['Prot Id']
        elif rv_len > bakta_len:
            new_coord, new_id = (*value['miniprot_df'][0], strand, contig), value['Prot Id']
        elif rv_len < bakta_len and abs(bakta_len - rv_len) > 550 and not complete:
            new_coord, new_id = (*value['case4'][0], bakta_strand, bakta_contig), value['Locus Tags'][0]
        else:
            new_coord, new_id = (*value['miniprot_df'][0], strand, contig), value['Prot Id']

        _insert_with_dedup(combined_dict, new_coord, new_id, reference_ids=reference_ids)


def resolve_multiple_overlaps(multiple_groups, miniprot_df, combined_dict, new_groups=None):
    """
    Resolve N-to-1 overlaps (multiple Bakta CDS vs one Miniprot hit).
    
    Pick the best overlapping Bakta CDS, then compare to Miniprot by length.
    """
    reference_ids = set(miniprot_df['Prot Id'])
    new_groups = new_groups or {}
    for key, value in multiple_groups.items():
        miniprot_range = value['miniprot_df'][0]
        prot_id = value['Prot Id']

        # Strand of the hit being resolved, not of the protein's first placement
        # (see resolve_single_overlaps).
        strand = miniprot_df.at[key, 'Strand'] if key in miniprot_df.index else '+'
        contig = value.get('contig', '')
        case4_strands = value.get('case4_strands') or [strand] * len(value['case4'])
        case4_contigs = value.get('case4_contigs') or [contig] * len(value['case4'])

        best = (None, None, -1, strand, contig)
        for i, r in enumerate(value['case4']):
            pct = get_overlap_percentage(r, miniprot_range)
            if pct > best[2]:
                tag = value['Locus Tags'][i] if i < len(value['Locus Tags']) else "Unknown"
                best = (r, tag, pct, case4_strands[i], case4_contigs[i])

        if not best[0]:
            continue

        # As in resolve_single_overlaps: only a fragmentary alignment loses to a
        # much longer Bakta call, and never to one that spans other reference
        # genes as well, which is several genes read as one. Rv3022c and
        # Rv3021c are one reading frame in CDC1551, and Rv3022c ends up here
        # rather than with the single overlaps because it shares four bases
        # with PE29.
        complete = (bool(miniprot_df.at[key, 'Complete'])
                    if key in miniprot_df.index and 'Complete' in miniprot_df.columns else False)
        spans_others = best[1] in new_groups and len(new_groups[best[1]]['Locus Tags']) > 1
        if complete or spans_others or _range_length(miniprot_range) >= _range_length(best[0]) or \
           (_range_length(best[0]) - _range_length(miniprot_range)) < 500:
            chosen_coord, chosen_id = (*miniprot_range, strand, contig), prot_id
        else:
            chosen_coord, chosen_id = (*best[0], best[3], best[4]), best[1]

        _insert_with_dedup(combined_dict, chosen_coord, chosen_id, threshold=550,
                           reference_ids=reference_ids)


def _is_reference_gene(gene_id, reference_ids=frozenset()):
    """Whether a gene in the dict is a reference protein's placement.

    This used to look for coordinates embedded in the ID, but clean_id strips
    them before any ID reaches the dict (since the March 2026 refactor), so no
    gene ever counted as a reference gene and any call of similar length could
    displace one: Bakta's version of Rv1614, one codon longer, replaced it.
    The callers pass the reference protein IDs instead.
    """
    return gene_id in reference_ids


def _insert_with_dedup(combined_dict, new_coord, new_id, threshold=550,
                       reference_ids=frozenset()):
    """Insert a coordinate into the dict, resolving overlaps with existing entries.
    
    Coordinates are (Start, End, Strand, Contig). Overlap checks only compare
    genes on the SAME strand of the same contig; anything else never conflicts.
    A reference gene already in the dict only gives way to a call that is
    longer by more than the threshold.
    """
    new_start, new_end = new_coord[0], new_coord[1]
    new_strand = new_coord[2] if len(new_coord) > 2 else '+'
    new_contig = _contig(new_coord)
    coords_to_add = {new_coord: new_id}
    coords_to_remove = set()

    for exist_coord, exist_id in list(combined_dict.items()):
        if new_coord == exist_coord:
            continue
        exist_start, exist_end = exist_coord[0], exist_coord[1]
        exist_strand = exist_coord[2] if len(exist_coord) > 2 else '+'
        # Only compare genes on the same strand of the same contig
        if new_strand != exist_strand or new_contig != _contig(exist_coord):
            continue
        if new_end <= exist_start or new_start >= exist_end:
            continue  # no overlap

        pct = get_overlap_percentage((new_start, new_end), (exist_start, exist_end))
        if pct >= 70:
            len_new = _range_length((new_start, new_end))
            len_exist = _range_length((exist_start, exist_end))
            if len_new > len_exist and (len_new - len_exist) > threshold:
                coords_to_remove.add(exist_coord)
            elif (len_exist - len_new) < threshold and not _is_reference_gene(exist_id, reference_ids):
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
