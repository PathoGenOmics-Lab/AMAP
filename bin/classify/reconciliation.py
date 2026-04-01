"""
Post-resolution steps: reconciliation pass, missing gene recovery,
deduplication, and final corrected TSV output.
"""

import re
import pandas as pd

from .ids import format_cluster_id
from .overlaps import get_subset_overlap_percentage


def parse_ref_coords_from_id(prot_id):
    """Extract 1-based start and end from IDs like 'Rv2427A_2725571-2726087'."""
    m = re.search(r'_(\d+)-(\d+)$', str(prot_id))
    if m:
        return int(m.group(1)), int(m.group(2))
    return None, None


def reconcile_coordinates(combined_dict, miniprot_df, bakta_df=None):
    """
    Replace coordinates with true reference coordinates when available.
    
    For each gene in the dict, check if the protein ID contains embedded
    coordinates (e.g., Rv1234_12345-67890). If the embedded coordinates
    differ from the current ones but overlap significantly, replace them.
    
    Safety checks:
    1. Only reconcile when lengths are comparable (within 2×) to prevent
       paralog-driven coordinate replacement.
    2. Remove genes placed at wrong positions (embedded coords don't match
       actual position), UNLESS they are multi-copy elements (>2 Miniprot
       hits), where placement at different positions is expected.
    """
    # Build ref coord map
    ref_coords = {}
    for prot_id in miniprot_df['Prot Id'].unique():
        s, e = parse_ref_coords_from_id(prot_id)
        if s and e:
            ref_coords[prot_id] = (s, e)

    # Count hits per protein to identify multi-copy elements
    hit_counts = miniprot_df.groupby('Prot Id').size()

    to_remove = set()
    to_add = {}
    for coord, gene_id in list(combined_dict.items()):
        # Check if THIS gene is a paralog placement: if its own ID has embedded
        # coordinates that don't match its current position, it was placed at
        # a paralogous region by Miniprot and should be removed.
        # Exception: skip check if ANY protein hitting this position in Miniprot
        # has embedded coords that match (proving this IS a real gene locus).
        own_s, own_e = parse_ref_coords_from_id(str(gene_id))
        if own_s is not None and own_e is not None:
            own_ref_len = own_e - own_s
            own_overlap = max(0, min(coord[1], own_e) - max(coord[0], own_s))
            if own_ref_len > 0 and own_overlap / own_ref_len < 0.5:
                # gene_id's embedded coords don't match → check if any OTHER
                # Miniprot hit at this position validates it as a real locus
                mp_at_pos = miniprot_df[
                    (miniprot_df['Start'] == coord[0]) & (miniprot_df['End'] == coord[1])
                ]
                position_validated = False
                for _, mp_row in mp_at_pos.iterrows():
                    mp_s, mp_e = parse_ref_coords_from_id(mp_row['Prot Id'])
                    if mp_s is not None and mp_e is not None:
                        mp_len = mp_e - mp_s
                        mp_ov = max(0, min(coord[1], mp_e) - max(coord[0], mp_s))
                        if mp_len > 0 and mp_ov / mp_len >= 0.5:
                            position_validated = True
                            break
                # Also validate if Bakta has a gene at these exact coordinates
                # (Bakta gene caller is independent of Miniprot, so a gene
                # with matching coords confirms a real coding region)
                if not position_validated and bakta_df is not None:
                    bk_at_pos = bakta_df[
                        (bakta_df['Start'] == coord[0]) &
                        (bakta_df['End'] == coord[1])
                    ]
                    if not bk_at_pos.empty:
                        position_validated = True
                if not position_validated:
                    to_remove.add(coord)
                    continue

        for ref_id, ref_coord in ref_coords.items():
            if ref_id in str(gene_id) or str(gene_id) in ref_id:
                if ref_coord != coord:
                    overlap = max(0, min(coord[1], ref_coord[1]) - max(coord[0], ref_coord[0]))
                    min_len = min(coord[1] - coord[0], ref_coord[1] - ref_coord[0])
                    max_len = max(coord[1] - coord[0], ref_coord[1] - ref_coord[0])
                    # Only reconcile if: significant overlap AND lengths are
                    # comparable (ref not >2× longer/shorter than current).
                    # This prevents paralog-driven coordinate replacement
                    # (e.g., Rv1003720new replacing Rv3508 coords).
                    if min_len > 0 and overlap / min_len > 0.5 and max_len < 2 * min_len:
                        to_remove.add(coord)
                        to_add[ref_coord] = ref_id
                        break

    # When removing paralog placements, check if a Bakta gene at that position
    # should be restored (e.g., PE_PGRS genes where Miniprot fails to align
    # the correct protein but Bakta calls the gene correctly).
    bakta_restore = {}
    if bakta_df is not None:
        for c in to_remove:
            if c not in to_add:  # Only for pure removals, not replacements
                bk_overlap = bakta_df[
                    (bakta_df['Start'] <= c[0] + 50) &
                    (bakta_df['End'] >= c[1] - 50)
                ]
                if not bk_overlap.empty:
                    bk = bk_overlap.iloc[0]
                    bk_coord = (bk['Start'], bk['End'])
                    if bk_coord not in combined_dict and bk_coord not in to_add:
                        bakta_restore[bk_coord] = bk['Locus Tag']

    for c in to_remove:
        combined_dict.pop(c, None)
    for c, v in to_add.items():
        if c not in combined_dict:
            combined_dict[c] = v
    for c, v in bakta_restore.items():
        if c not in combined_dict:
            combined_dict[c] = v


def recover_missing_genes(combined_dict, miniprot_df):
    """
    Add reference genes that are completely absent from the coordinate dict.
    
    Safety: for each candidate, verify that the Miniprot alignment position
    matches the coordinates embedded in the protein ID. If a protein named
    'Rv1150_1278269-1278820' is placed at position 1164572-1165123, it
    means Miniprot mapped it to a paralogous region — skip it.
    """
    all_ref_ids = set(miniprot_df['Prot Id'].unique())
    included_ids = set(combined_dict.values())
    missing = all_ref_ids - included_ids

    if not missing:
        return

    df_missing = miniprot_df[miniprot_df['Prot Id'].isin(missing)]
    if df_missing.empty:
        return

    try:
        best = df_missing.groupby('Prot Id').apply(
            lambda g: g.loc[g['Match Nt'].idxmax()], include_groups=False).reset_index()
    except TypeError:
        best = df_missing.groupby('Prot Id').apply(
            lambda g: g.loc[g['Match Nt'].idxmax()]).reset_index(drop=True)

    for _, row in best.iterrows():
        coord = (row['Start'], row['End'])
        if coord not in combined_dict:
            # Verify placement matches embedded coordinates in Prot Id.
            # If the protein has embedded coords (e.g., Rv1150_1278269-1278820),
            # the alignment should be near those coords, not at a paralog.
            ref_s, ref_e = parse_ref_coords_from_id(row['Prot Id'])
            if ref_s is not None and ref_e is not None:
                ref_len = ref_e - ref_s
                placement_overlap = max(0, min(coord[1], ref_e) - max(coord[0], ref_s))
                if ref_len > 0 and placement_overlap / ref_len < 0.5:
                    # Placement doesn't match embedded coords → paralog mapping
                    continue
            combined_dict[coord] = row['Prot Id']


def deduplicate(combined_dict):
    """Remove duplicate gene entries that are subsets of other entries."""
    all_vals = list(combined_dict.values())
    duplicates = {x for x in all_vals if all_vals.count(x) > 1}
    if not duplicates:
        return

    coords_by_val = {
        val: [k for k, v in combined_dict.items() if v == val]
        for val in duplicates
    }

    to_delete = set()
    for dup_val, dup_coords in coords_by_val.items():
        for dc in dup_coords:
            for coord, value in combined_dict.items():
                if value not in duplicates:
                    if get_subset_overlap_percentage(dc, coord) == 100:
                        to_delete.add(dc)

    for c in to_delete:
        combined_dict.pop(c, None)


def write_corrected_tsv(combined_dict, miniprot_df, bakta_df, base_name):
    """
    Write the final corrected TSV with gene IDs, strands, coordinates, and Rv mappings.
    
    Filters out entries with start=0 or end=0, assigns strands from
    miniprot/bakta matches, and formats IDs.
    """
    # Sort and filter
    sorted_dict = dict(sorted(combined_dict.items()))
    filtered = {
        (s, e): v for (s, e), v in sorted_dict.items()
        if s != 0 and e != 0
    }

    singles = []
    multis = []
    dup_counter = 1

    for (start, end), value in filtered.items():
        strand = _find_strand(start, end, miniprot_df, bakta_df)

        mp_matches = miniprot_df[(miniprot_df['Start'] == start) & (miniprot_df['End'] == end)]
        bk_matches = bakta_df[(bakta_df['Start'] == start) & (bakta_df['End'] == end)]

        if not bk_matches.empty:
            cluster_id = format_cluster_id(bk_matches.iloc[0]['Locus Tag'], base_name)
            rv_ids = ", ".join(mp_matches['Prot Id'].unique()) if not mp_matches.empty else ""
            singles.append([cluster_id, strand, start, end, rv_ids])
        else:
            rv_ids = mp_matches['Prot Id'].unique() if not mp_matches.empty else []
            if len(rv_ids) == 0:
                cluster_id = format_cluster_id(value, base_name)
                singles.append([cluster_id, strand, start, end, ""])
            elif len(rv_ids) < 2:
                cluster_id = format_cluster_id(rv_ids[0], base_name)
                singles.append([cluster_id, strand, start, end, ", ".join(rv_ids)])
            else:
                cluster_id = f"{base_name}_dup_{dup_counter}"
                dup_counter += 1
                multis.append([cluster_id, strand, start, end, ", ".join(rv_ids)])

    final = sorted(singles + multis, key=lambda x: x[2])

    output_path = f"{base_name}_corrected.tsv"
    with open(output_path, 'w') as f:
        for row in final:
            f.write('\t'.join(map(str, row)) + '\n')

    print(f"Processing complete. Final gene list saved to {output_path}")
    print(f"Total genes in final set: {len(final)}")


def _find_strand(start, end, miniprot_df, bakta_df):
    """Determine strand for a coordinate pair from miniprot or bakta."""
    mp = miniprot_df[(miniprot_df['Start'] == start) & (miniprot_df['End'] == end)]
    if not mp.empty:
        return mp.iloc[0]['Strand']

    bk = bakta_df[(bakta_df['Start'] == start) & (bakta_df['End'] == end)]
    if not bk.empty:
        return bk.iloc[0]['Strand']

    # Try overlapping
    mp_ov = miniprot_df[(miniprot_df['Start'] <= end) & (miniprot_df['End'] >= start)]
    if not mp_ov.empty:
        return mp_ov.iloc[0]['Strand']

    bk_ov = bakta_df[(bakta_df['Start'] <= end) & (bakta_df['End'] >= start)]
    if not bk_ov.empty:
        return bk_ov.iloc[0]['Strand']

    print(f"WARNING: No strand found for gene at {start}-{end}, defaulting to '+'")
    return "+"
