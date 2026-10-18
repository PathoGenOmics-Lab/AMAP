"""
Post-resolution steps: missing gene recovery, deduplication, and the final
corrected TSV.
"""

from collections import Counter

import pandas as pd
from intervaltree import IntervalTree

from .ids import format_cluster_id
from .overlaps import get_subset_overlap_percentage, _contig


def _frame(coord):
    """Reading frame of a (Start, End, Strand, ...) coordinate, from its 5' end."""
    return (coord[1] if len(coord) > 2 and coord[2] == '-' else coord[0]) % 3


def recover_missing_genes(combined_dict, miniprot_df, min_cover=0.5):
    """
    Add back reference proteins that no gene in the dict accounts for.

    A placement is accounted for when a gene already sits on it: same contig
    and strand, the same reading frame, covering at least min_cover of it, and
    not known to be another reference gene. A protein is missing when none of
    its placements is, and then its best one goes in.

    'Missing' used to mean that the protein's name was not among the dict's
    values. Genes matched through case 2 and 3 are stored under Bakta locus
    tags, so every protein looked missing (4034 of 4034 on H37Rv), and what
    kept that from duplicating every gene was an exact-coordinate check that a
    placement one base off went straight past. The frame matters because genes
    that genuinely overlap sit in different frames (Rv0724A and Rv0725c), and
    one does not account for the other. Neither does a different reference
    gene in the same frame: E. coli nests hokC inside mokC and b4706 inside
    b2537, and both genes of each pair are real.
    """
    reference_ids = set(miniprot_df['Prot Id'])
    placed = miniprot_df[(miniprot_df['Start'] < miniprot_df['End'])
                         & ~((miniprot_df['Match Nt'] == 0) & (miniprot_df['Num nt'] == 0))]
    contig = placed['Contig'] if 'Contig' in placed.columns else pd.Series('', index=placed.index)
    coords = {idx: (r['Start'], r['End'], r.get('Strand', '+'), contig[idx])
              for idx, r in placed.iterrows()}

    # The reference genes each gene of the dict stands for: itself, or for a
    # Bakta call, the proteins placed at exactly its coordinates. None known
    # means it can be any of them.
    at_coord = {}
    for idx, c in coords.items():
        at_coord.setdefault(c, set()).add(placed.at[idx, 'Prot Id'])

    def labels(coord, gene_id):
        return {gene_id} if gene_id in reference_ids else at_coord.get(coord, set())

    trees = {}

    def where(coord):
        return _contig(coord), coord[2] if len(coord) > 2 else '+'

    def add_to_tree(coord, gene_id):
        tree = trees.setdefault(where(coord), IntervalTree())
        tree[coord[0]:coord[1] + 1] = (_frame(coord), labels(coord, gene_id))

    def accounted(coord, prot_id):
        length = coord[1] - coord[0] + 1
        for iv in trees.get(where(coord), IntervalTree())[coord[0]:coord[1] + 1]:
            frame, genes = iv.data
            if frame != _frame(coord) or (genes and prot_id not in genes):
                continue
            if min(coord[1] + 1, iv.end) - max(coord[0], iv.begin) >= min_cover * length:
                return True
        return False

    for coord, gene_id in combined_dict.items():
        add_to_tree(coord, gene_id)

    present = set()
    for idx, c in coords.items():
        prot_id = placed.at[idx, 'Prot Id']
        if c in combined_dict or accounted(c, prot_id):
            present.add(prot_id)

    # Best placements first, so a missing protein goes in where it fits best.
    for idx in placed.sort_values('Match Nt', ascending=False).index:
        prot_id = placed.at[idx, 'Prot Id']
        coord = coords[idx]
        if prot_id in present or coord in combined_dict or accounted(coord, prot_id):
            continue
        combined_dict[coord] = prot_id
        add_to_tree(coord, prot_id)
        present.add(prot_id)


def deduplicate(combined_dict):
    """Remove duplicate gene entries that are subsets of other entries.
    
    When a gene appears at multiple coordinates, remove entries that are
    100% subsets of a DIFFERENT gene's coordinates. But always keep at
    least ONE entry per gene to prevent total gene loss.
    """
    # Counter instead of list.count() per element, which was O(genes^2).
    val_counts = Counter(combined_dict.values())
    duplicates = {v for v, n in val_counts.items() if n > 1}
    if not duplicates:
        return

    # Single pass instead of one full scan per duplicated gene.
    coords_by_val = {}
    for k, v in combined_dict.items():
        if v in duplicates:
            coords_by_val.setdefault(v, []).append(k)

    # The set of non-duplicated coordinates is the same for every comparison,
    # so build it once rather than re-filtering combined_dict in the inner loop.
    non_dup_coords = [k for k, v in combined_dict.items() if v not in duplicates]

    to_delete = set()
    for dup_val, dup_coords in coords_by_val.items():
        # Collect which coords of this gene are subsets of other genes
        subset_coords = set()
        for dc in dup_coords:
            for coord in non_dup_coords:
                if _contig(dc) == _contig(coord) and get_subset_overlap_percentage(dc, coord) == 100:
                    subset_coords.add(dc)
                    break

        # Only delete subset coords if at least one coord survives
        # This prevents total gene loss when ALL coords are subsets
        surviving = [c for c in dup_coords if c not in subset_coords]
        if surviving:
            to_delete.update(subset_coords)
        elif subset_coords:
            # All are subsets — keep the largest one
            best = max(dup_coords, key=lambda c: c[1] - c[0])
            to_delete.update(subset_coords - {best})

    for c in to_delete:
        combined_dict.pop(c, None)


def write_corrected_tsv(combined_dict, miniprot_df, bakta_df, base_name):
    """
    Write the final corrected TSV with gene IDs, strands, coordinates, and Rv mappings.
    
    Filters out entries with start=0 or end=0, assigns strands from
    miniprot/bakta matches, and formats IDs.
    """
    # Contigs in the order the assembly has them (Bakta lists them that way),
    # then by start coordinate; filter invalid entries.
    rank = {}
    if 'Contig' in bakta_df.columns:
        rank = {c: i for i, c in enumerate(dict.fromkeys(bakta_df['Contig']))}

    def order(contig, start, end=0):
        return (rank.get(contig, len(rank)), contig, start, end)

    sorted_dict = dict(sorted(combined_dict.items(),
                              key=lambda x: order(_contig(x[0]), x[0][0], x[0][1])))
    filtered = {
        k: v for k, v in sorted_dict.items()
        if k[0] != 0 and k[1] != 0
    }

    mp_contig = miniprot_df['Contig'] if 'Contig' in miniprot_df.columns else ''
    bk_contig = bakta_df['Contig'] if 'Contig' in bakta_df.columns else ''

    singles = []
    multis = []
    dup_counter = 1

    for coord_key, value in filtered.items():
        start, end = coord_key[0], coord_key[1]
        contig = _contig(coord_key)
        # Use strand from key if available, otherwise derive from data
        strand = coord_key[2] if len(coord_key) > 2 else _find_strand(start, end, miniprot_df, bakta_df)

        mp_matches = miniprot_df[(miniprot_df['Start'] == start) & (miniprot_df['End'] == end)
                                 & (mp_contig == contig)]
        bk_matches = bakta_df[(bakta_df['Start'] == start) & (bakta_df['End'] == end)
                              & (bk_contig == contig)]

        if not bk_matches.empty:
            cluster_id = format_cluster_id(bk_matches.iloc[0]['Locus Tag'], base_name)
            rv_ids = ", ".join(mp_matches['Prot Id'].unique()) if not mp_matches.empty else ""
            singles.append([cluster_id, strand, start, end, rv_ids, contig])
        else:
            rv_ids = mp_matches['Prot Id'].unique() if not mp_matches.empty else []
            if len(rv_ids) == 0:
                cluster_id = format_cluster_id(value, base_name)
                singles.append([cluster_id, strand, start, end, "", contig])
            elif len(rv_ids) < 2:
                cluster_id = format_cluster_id(rv_ids[0], base_name)
                singles.append([cluster_id, strand, start, end, ", ".join(rv_ids), contig])
            else:
                cluster_id = f"{base_name}_dup_{dup_counter}"
                dup_counter += 1
                multis.append([cluster_id, strand, start, end, ", ".join(rv_ids), contig])

    # The contig is the sixth column: check_cds.py translates each gene from
    # it, where it used to guess the contig from the coordinates alone.
    final = sorted(singles + multis, key=lambda x: order(x[5], x[2]))

    output_path = f"{base_name}_corrected.tsv"
    with open(output_path, 'w') as f:
        for row in final:
            f.write('\t'.join(map(str, row)) + '\n')

    print(f"Processing complete. Final gene list saved to {output_path}")
    print(f"Total genes in final set: {len(final)}")


def _find_strand(start, end, miniprot_df, bakta_df):
    """Determine strand for a coordinate pair from miniprot or bakta.
    
    Priority: exact match > best overlapping (by overlap length).
    """
    mp = miniprot_df[(miniprot_df['Start'] == start) & (miniprot_df['End'] == end)]
    if not mp.empty:
        return mp.iloc[0]['Strand']

    bk = bakta_df[(bakta_df['Start'] == start) & (bakta_df['End'] == end)]
    if not bk.empty:
        return bk.iloc[0]['Strand']

    # Try overlapping — pick the one with the largest overlap
    def _best_overlap_strand(df):
        ov = df[(df['Start'] <= end) & (df['End'] >= start)].copy()
        if ov.empty:
            return None
        ov['_ov_len'] = ov.apply(
            lambda r: min(end, r['End']) - max(start, r['Start']) + 1, axis=1)
        return ov.loc[ov['_ov_len'].idxmax(), 'Strand']

    mp_strand = _best_overlap_strand(miniprot_df)
    if mp_strand is not None:
        return mp_strand

    bk_strand = _best_overlap_strand(bakta_df)
    if bk_strand is not None:
        return bk_strand

    print(f"WARNING: No strand found for gene at {start}-{end}, defaulting to '+'")
    return "+"
