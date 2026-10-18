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
    recover_missing_genes, deduplicate, write_corrected_tsv
)
from classify.plotting import generate_overlap_plot
from classify.ids import format_cluster_id
from classify.synteny import drop_paralog_placements
from classify.exact import place_exact_matches
from classify.pieces import add_contig_end_pieces

__author__ = 'Paula Ruiz-Rodriguez'
__credits__ = ['Paula Ruiz-Rodriguez', 'PathoGenOmics']
__license__ = ['GPL 3']
__version__ = '1.0.1'
__maintainer__ = '@paururo'
__email__ = 'paula.ruiz.rodriguez@csic.es'
__status__ = 'finished'
__lab__ = 'PathoGenOmics (I2SysBio)'


def create_coords_dict(df, id_col):
    """Create a (Start, End, Strand, Contig)→ID dict.
    
    Keys include Strand to prevent collisions between genes at the same
    coordinates on opposite strands (e.g., antisense overlapping genes), and
    Contig because every contig of a draft assembly starts at 1.
    First occurrence wins for duplicate coordinates.
    """
    result = {}
    for _, row in df.iterrows():
        strand = row.get('Strand', '+')
        key = (row['Start'], row['End'], strand, row.get('Contig', ''))
        if key not in result:
            result[key] = row[id_col]
    return result


def match_bakta_contigs(bakta_df, miniprot_df):
    """Make sure Bakta and miniprot name the sequences alike.

    Bakta renames every sequence to contig_N unless it runs with
    --keep-contig-headers, and miniprot keeps the FASTA headers. With nothing
    in common every comparison would be between different contigs and the
    whole classification would come out wrong without a word. One sequence
    each can only be the same sequence; anything else is a hard error. With
    no CDS from Bakta there is no name to compare.
    """
    bakta_contigs = set(bakta_df['Contig'])
    mp_contigs = set(miniprot_df.loc[miniprot_df['Contig'] != '*', 'Contig'])
    if not bakta_contigs or not mp_contigs or bakta_contigs & mp_contigs:
        return bakta_df
    if len(bakta_contigs) == 1 and len(mp_contigs) == 1:
        return bakta_df.assign(Contig=next(iter(mp_contigs)))
    sys.exit(f"ERROR: Bakta names the sequences {sorted(bakta_contigs)[:3]}... but the "
             f"genome has {sorted(mp_contigs)[:3]}...; run Bakta with --keep-contig-headers.")


def main():
    parser = argparse.ArgumentParser(
        description="Classify Bakta and Miniprot records and resolve overlaps.")
    parser.add_argument('-b', '--bakta', required=True, help="Bakta TSV file.")
    parser.add_argument('-p', '--miniprot', required=True, help="Miniprot PAF file.")
    parser.add_argument('-bb', '--badbakta', required=True, help="Miniprot 'bad bakta' PAF.")
    parser.add_argument('-c', '--cds', required=True, help="Reference CDS PAF for Miniprot.")
    parser.add_argument('-g', '--genome', help="Genome FASTA. With --ref-proteome, the proteins "
                                               "miniprot did not place exactly are looked for in it.")
    parser.add_argument('-r', '--ref-proteome', help="Reference proteome the miniprot PAF was made from.")
    parser.add_argument('--pieces', help="Miniprot PAF of what contig ends cut off the placements "
                                         "in the first one (contig_end_pieces.py).")
    parser.add_argument('--sample', help="The sample's ID in the sample sheet, for the summary. "
                                         "Default: the name of the Bakta table.")
    parser.add_argument('--lineage', help="The lineage typing called the sample, for the summary. "
                                          "Default: Untyped.")
    args = parser.parse_args()

    base_name = os.path.basename(args.bakta).split('.')[0]

    # --- 1. Read inputs ---
    bakta_df = read_bakta(args.bakta)
    miniprot_df = read_miniprot(args.miniprot, ref_coords=True)
    # miniprot cannot seed an alignment for the smallest proteins, or puts them
    # on a near-identical paralog; place the ones the genome encodes exactly,
    # or for the shortest, with a few substitutions.
    if args.genome and args.ref_proteome:
        miniprot_df = place_exact_matches(miniprot_df, args.genome, args.ref_proteome)
    # A gene a contig end cut is placed by its larger piece alone; the others
    # come from a second pass over what that placement lacks.
    if args.pieces:
        miniprot_df = add_contig_end_pieces(miniprot_df, args.pieces)
    # A protein with a paralog is placed at both loci; keep only the placements
    # that do not land on a locus another protein syntenically owns.
    miniprot_df = drop_paralog_placements(miniprot_df)
    bakta_df_bad = read_miniprot_unfiltered(args.badbakta)
    bakta_df = match_bakta_contigs(bakta_df, miniprot_df)

    # --- 2. Classify cases ---
    case2, case3 = classify_case2_case3(bakta_df, miniprot_df, base_name)
    case4, case1 = classify_case1_case4(bakta_df, miniprot_df, base_name)
    rv_missing = find_rv_missing(miniprot_df, bakta_df, base_name)
    case1_keep, case1_remove = filter_case1(case1, bakta_df_bad, args.cds)

    # --- 3. Write summary ---
    write_summary(base_name,
                  [case2, case3, case4, case1, rv_missing, bakta_df, case1_remove, case1_keep],
                  ["case2", "case3", "case4", "case1_initial", "rv_missing",
                   "total_cds", "case1_removed", "case1_kept"],
                  sample=args.sample, lineage=args.lineage)

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
    n_before = 0
    for df, id_col in [(case2, 'Locus Tag'), (case3, 'Locus Tag'),
                        (rv_missing, 'Prot Id'), (case1_keep, 'Locus Tag')]:
        new_entries = create_coords_dict(df, id_col)
        n_before += len(new_entries)
        combined.update(new_entries)
    if len(combined) < n_before:
        n_lost = n_before - len(combined)
        print(f"WARNING: {n_lost} genes lost to (Start, End, Strand, Contig) collisions, "
              f"i.e. listed twice at the same place.")

    # --- 7. Resolve overlaps ---
    # NOTE: case4_mp_overlaps stores case4 indices, but create_overlap_groups
    # uses them to index miniprot_df for coordinates. The coordinates in
    # b_overlap_groups['miniprot_df'] are therefore wrong. However, they are
    # never accessed — only 'Locus Tags' (count of overlapping Miniprot hits
    # per case4 gene) is used via new_groups. Safe but fragile.
    b_overlap_groups = create_overlap_groups(case4_mp_overlaps, case4, miniprot_df)
    new_groups = {v['Prot Id']: v for _, v in b_overlap_groups.items()}

    resolve_single_overlaps(single_groups, miniprot_df, combined, new_groups)
    resolve_multiple_overlaps(multiple_groups, miniprot_df, combined, new_groups)

    # --- 8. Recovery + dedup ---
    recover_missing_genes(combined, miniprot_df)
    deduplicate(combined)

    # --- 9. Write final output ---
    write_corrected_tsv(combined, miniprot_df, bakta_df, base_name)


if __name__ == "__main__":
    main()
