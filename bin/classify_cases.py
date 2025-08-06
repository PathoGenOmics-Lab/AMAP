#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import os
from intervaltree import Interval, IntervalTree
import plotly.graph_objects as go

__author__ = 'Paula Ruiz-Rodriguez'
__credits__ = ['Paula Ruiz-Rodriguez', 'PathoGenOmics']
__license__ = ['GPL 3']
__version__ = '0.1.0'
__maintainer__ = '@paururo'
__email__ = 'paula.ruiz-rodriguez@uv.es'
__status__ = 'finished'
__lab__ = 'PathoGenOmics (I2SysBio)'

def clean_id(gene_id):
    """
    Cleans a gene ID based on the 4-part rule.
    If the ID split by '_' has 4 elements, it removes the last one.
    Example: 'L8_G26689_Rv0609A_703830-704057' becomes 'L8_G26689_Rv0609A'.
    """
    gene_id = str(gene_id)
    parts = gene_id.split('_')
    if len(parts) == 4:
        # If there are exactly 4 parts, join the first three.
        return '_'.join(parts[:-1])
    else:
        # Otherwise, return the original ID.
        return gene_id

def format_final_id(gene_id):
    """
    Formats a gene ID to keep only the first 3 parts.
    Example: 'L8_G26689_Rv0609A_703830-704057' becomes 'L8_G26689_Rv0609A'.
    If it has fewer than 3 parts, it returns the original ID.
    """
    parts = str(gene_id).split('_')
    if len(parts) >= 3:
        # Join the first three parts.
        return '_'.join(parts[:3])
    # Otherwise, return the ID as is.
    return gene_id

def read_to_pandaBakta(bakta_path):
    """
    Reads and processes a genomic data file in TSV format from Bakta.
    """
    bakta_df = pd.read_csv(bakta_path, comment='#', sep='\t', header=None)
    bakta_df.columns = ["Sequence Id", "Type", "Start", "End", "Strand", "Locus Tag", "Gene", "Product", "DbXrefs"]
    bakta_df['Locus Tag'] = bakta_df['Locus Tag'].apply(clean_id)
    bakta_cds = bakta_df[bakta_df['Type'] == 'cds'].copy()
    bakta_cds['Length_nt'] = bakta_cds['End'] - bakta_cds['Start']
    last_bakta = bakta_cds[["Locus Tag", "Start", "End", "Strand", 'Length_nt']]
    return last_bakta

def custom_filter(group):
    """
    Applies a custom filter to a pandas group to find perfect matches.
    """
    if any(group['Num nt'] - group['Match Nt'] == 0):
        return group[group['Num nt'] - group['Match Nt'] == 0]
    else:
        return group

def read_to_pandaMiniprot(miniprot_path):
    """
    Reads and processes a miniprot data file in TSV format.
    """
    miniprot_df = pd.read_csv(miniprot_path, sep='\t', header=None)
    miniprot_df.columns = [
        "Prot Id", "Prot len", "Prot Start", "ProtStop", "Strand", "fasta Tag", "fasta len",
        "Start", "End", "Match Nt", "Num nt", "mapping qual", "AS", "MS", "NP", "FS", "ST",
        "DA", "DO", "CG", "CS"
    ]
    miniprot_df['Prot Id'] = miniprot_df['Prot Id'].apply(clean_id)
    miniprot_df['Start'] = miniprot_df['Start'] + 1
    last_miniprot = miniprot_df[["Prot Id", "Prot len", "Strand", "Start", "End", "Match Nt", "Num nt"]]
    # Apply the filter to keep only the best alignments per protein ID
    filtered_miniprot = last_miniprot.groupby('Prot Id').apply(custom_filter, include_groups=False).reset_index()
    return filtered_miniprot

def read_to_pandaBadbakta(miniprot_path):
    """
    Reads a miniprot file, intended for 'bad' or alternative alignments, into a DataFrame.
    """
    miniprot_df = pd.read_csv(miniprot_path, sep='\t', header=None)
    miniprot_df.columns = [
        "Prot Id", "Prot len", "Prot Start", "ProtStop", "Strand", "fasta Tag", "fasta len",
        "Start", "End", "Match Nt", "Num nt", "mapping qual", "AS", "MS", "NP", "FS", "ST",
        "DA", "DO", "CG", "CS"
    ]
    miniprot_df['Prot Id'] = miniprot_df['Prot Id'].apply(clean_id)
    miniprot_df['Start'] = miniprot_df['Start'] + 1
    last_miniprot = miniprot_df[["Prot Id", "Prot len", "Strand", "Start", "End", "Match Nt", "Num nt"]]
    return last_miniprot

def process_case2_and_case3(bakta_df, miniprot_df, base_name=None):
    """
    Merges two dataframes and separates them into 'case 2' (perfect coordinate match)
    and 'case 3' (imperfect coordinate match).
    """
    case2_merged = pd.merge(bakta_df, miniprot_df, on=['Start', 'End'])
    case2_merged['Difference'] = (case2_merged['Length_nt'] - case2_merged['Num nt'])

    case_2 = case2_merged[case2_merged['Difference'] == 0]
    case3 = case2_merged[case2_merged['Difference'] != 0]

    case_2 = case_2[["Locus Tag", "Start", "End", "Strand_x", 'Length_nt']].drop_duplicates()

    if base_name:
        case3.to_csv(f"{base_name}-case3.tsv", sep='\t', index=False)
        case_2.to_csv(f"{base_name}-case2.tsv", sep='\t', index=False)

    return case_2, case3

def process_cases(bakta_df, miniprot_df, base_name=None):
    """
    Processes and separates data into 'case 1' (no overlap) and 'case 4' (coordinate overlap).
    """
    bakta_df['Length_nt'] = bakta_df['End'] - bakta_df['Start']
    merged_df = pd.merge(bakta_df, miniprot_df, on=['Start', 'End'], how='inner')

    common_columns = [col for col in bakta_df.columns if col in merged_df.columns]
    merged_indicator = bakta_df.merge(merged_df, on=common_columns, how='left', indicator=True)
    rows_not_in_merged = merged_indicator[merged_indicator['_merge'] == 'left_only'].drop('_merge', axis=1)
    rows_not_in_merged = rows_not_in_merged[["Locus Tag", "Start", "End", "Strand", 'Length_nt']]

    # Check for any positional overlap
    rows_not_in_merged['Overlaps'] = False
    for index, row in rows_not_in_merged.iterrows():
        has_overlap = ((miniprot_df['Start'] <= row['End']) & (miniprot_df['End'] >= row['Start'])).any()
        rows_not_in_merged.at[index, 'Overlaps'] = has_overlap

    case4 = rows_not_in_merged[rows_not_in_merged['Overlaps'] == True]
    case1 = rows_not_in_merged[rows_not_in_merged['Overlaps'] == False]

    case1 = case1[["Locus Tag", "Start", "End", "Strand", 'Length_nt']]
    case4 = case4[["Locus Tag", "Start", "End", "Strand", 'Length_nt']]

    if base_name:
        case4.to_csv(f"{base_name}-case4.tsv", sep='\t', index=False)
        case1.to_csv(f"{base_name}-case1.tsv", sep='\t', index=False)

    return case4, case1

def get_non_overlapping_rows(miniprot_df, bakta_df, base_name=None):
    """
    Identifies rows in miniprot_df that do not overlap with any row in bakta_df.
    These are considered 'rv_missing' from the Bakta annotation.
    """
    miniprot_df['Overlaps'] = False
    for index, row in miniprot_df.iterrows():
        has_overlap = ((bakta_df['Start'] <= row['End']) & (bakta_df['End'] >= row['Start'])).any()
        miniprot_df.at[index, 'Overlaps'] = has_overlap

    no_overlap_with_bakta = miniprot_df[miniprot_df['Overlaps'] == False]
    # Filter out zero-length alignments
    no_overlap_with_bakta = no_overlap_with_bakta[~((no_overlap_with_bakta['Match Nt'] == 0) & (no_overlap_with_bakta['Num nt'] == 0))]

    if base_name:
        no_overlap_with_bakta.to_csv(f"{base_name}-rvmissing.tsv", sep='\t', index=False)

    return no_overlap_with_bakta

def create_and_save_summary(base_name, data_frames, case_types):
    """
    Creates and saves a summary DataFrame to a TSV file.
    """
    summary_data = {
        "Type": case_types,
        "Number": [len(df) for df in data_frames],
        "Sequence": base_name,
        "Lineage": base_name.split("_")[0]
    }
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(f"{base_name}-summary.tsv", sep='\t', index=False)

def calculate_overlap_length(row1, df2):
    """
    Calculates the total length of overlaps between a given row and all rows in another DataFrame.
    """
    def overlap_len(start1, end1, start2, end2):
        return max(0, min(end1, end2) - max(start1, start2))

    total_overlap = df2.apply(lambda row2: overlap_len(row1['Start'], row1['End'], row2['Start'], row2['End']), axis=1).sum()
    return total_overlap

def get_case1(case1, bakta_df_bad, cds_path):
    """
    Further processes 'case 1' data to resolve potential misclassifications based on overlap percentages.
    """
    case1['Overlap_Length'] = case1.apply(lambda row: calculate_overlap_length(row, bakta_df_bad), axis=1)
    case1['Overlap_Percentage'] = (case1['Overlap_Length'] / case1['Length_nt']) * 100

    case1_remove = case1[case1['Overlap_Percentage'] >= 80]
    case1_check = case1[case1['Overlap_Percentage'] < 80]

    locus_ref_df = read_to_pandaBadbakta(cds_path)
    locus_tags_to_check = case1_check['Locus Tag'].tolist()

    matching_rows = locus_ref_df[locus_ref_df['Prot Id'].isin(locus_tags_to_check)].copy()
    matching_rows['match_per'] = 100 - ((matching_rows['Num nt'] - matching_rows['Match Nt']) / matching_rows['Num nt']) * 100
    matching_rows['match_per'] = matching_rows['match_per'].fillna(0)

    case1_remove_too = matching_rows[matching_rows['match_per'] >= 80]
    case1_keep = matching_rows[matching_rows['match_per'] < 80]

    locus_case1_remove_too = case1_remove_too['Prot Id'].tolist()
    locus_case1_keep = case1_keep['Prot Id'].tolist()

    case1_remove2 = case1_check[case1_check['Locus Tag'].isin(locus_case1_remove_too)]
    case1_keep2 = case1_check[case1_check['Locus Tag'].isin(locus_case1_keep)]

    final_case1_remove = pd.concat([case1_remove, case1_remove2])
    return case1_keep2, final_case1_remove

def get_overlap_tree(tree_df, overlap_df, tree_id_col, overlap_id_col):
    """
    Builds an Interval Tree and finds overlaps between two DataFrames.
    """
    tree = IntervalTree()
    for idx, row in tree_df.iterrows():
        start, end = row['Start'], row['End']
        if start <= end:
            tree[start:end] = (idx, row[tree_id_col])

    overlap_counts = {}
    for idx, row in overlap_df.iterrows():
        overlaps = tree[row['Start']:row['End']]
        for overlap in overlaps:
            tree_idx, tree_id = overlap.data
            if tree_idx not in overlap_counts:
                overlap_counts[tree_idx] = {'count': 1, 'Prot Id': tree_id, 'Locus Tags': [row[overlap_id_col]]}
            else:
                overlap_counts[tree_idx]['count'] += 1
                overlap_counts[tree_idx]['Locus Tags'].append(row[overlap_id_col])
    return overlap_counts

def create_overlap_groups(filtered_overlaps, case4_df, miniprot_df):
    """Creates a dictionary of overlapping gene groups for visualization."""
    overlaps_df = pd.DataFrame.from_dict(filtered_overlaps, orient='index')
    overlap_groups = {}
    for i, row in overlaps_df.iterrows():
        prot_id = row['Prot Id']
        locus_tags = row['Locus Tags']

        case4_rows = case4_df[case4_df['Locus Tag'].isin(locus_tags)]
        case4_coords = [(r['Start'], r['End']) for _, r in case4_rows.iterrows()]

        miniprot_coords = [(miniprot_df.at[i, 'Start'], miniprot_df.at[i, 'End'])]

        overlap_groups[i] = {
            'case4': case4_coords,
            'miniprot_df': miniprot_coords,
            'Prot Id': prot_id,
            'Locus Tags': locus_tags
        }
    return overlap_groups

def generate_plot(overlap_groups, output_path):
    """Generates an interactive HTML plot of gene overlaps."""
    fig = go.Figure()
    y_pos_counter = 0

    for i, (key, group) in enumerate(overlap_groups.items()):
        prot_id = group['Prot Id']
        locus_tags = group['Locus Tags']

        case4_x, case4_y, case4_text = [], [], []
        for idx, (start, end) in enumerate(group['case4']):
            case4_x.extend([start, end, None])
            case4_y.extend([y_pos_counter, y_pos_counter, None])
            tag = locus_tags[idx] if idx < len(locus_tags) else "Unknown"
            case4_text.extend([tag, tag, None])

        fig.add_trace(go.Scatter(
            x=case4_x, y=case4_y, mode='lines+markers',
            line=dict(color='blue', width=6), marker=dict(size=10),
            name=f"Bakta: {locus_tags[0] if locus_tags else ''}", hoverinfo='text', text=case4_text
        ))

        miniprot_x, miniprot_y, miniprot_text = [], [], []
        for start, end in group['miniprot_df']:
            miniprot_x.extend([start, end, None])
            miniprot_y.extend([y_pos_counter, y_pos_counter, None])
            miniprot_text.extend([prot_id, prot_id, None])

        fig.add_trace(go.Scatter(
            x=miniprot_x, y=miniprot_y, mode='lines+markers',
            line=dict(color='red', width=10), marker=dict(size=14),
            name=f"Miniprot: {prot_id}", hoverinfo='text', text=miniprot_text, opacity=0.7
        ))

        y_pos_counter += 1

    fig.update_layout(
        title='Visualization of Overlapping Gene Annotations',
        xaxis_title='Genomic Coordinate',
        yaxis_title='Overlap Group',
        showlegend=False
    )
    fig.write_html(output_path)


def main():
    parser = argparse.ArgumentParser(description="Classifies Bakta and Miniprot records and resolves overlaps.")
    parser.add_argument('-b', '--bakta', required=True, help="Path to the Bakta TSV file.")
    parser.add_argument('-p', '--miniprot', required=True, help="Path to the Miniprot output TSV file.")
    parser.add_argument('-bb', '--badbakta', required=True, help="Path to the Miniprot 'bad bakta' alignments.")
    parser.add_argument('-c', '--cds', required=True, help="Path to the reference CDS file for Miniprot.")
    args = parser.parse_args()

    base_name = os.path.basename(args.bakta).split('.')[0]

    bakta_df = read_to_pandaBakta(args.bakta)
    miniprot_df = read_to_pandaMiniprot(args.miniprot)
    bakta_df_bad = read_to_pandaBadbakta(args.badbakta)

    case_2, case3 = process_case2_and_case3(bakta_df, miniprot_df, base_name)
    case4, case1 = process_cases(bakta_df, miniprot_df, base_name)
    rv_missing = get_non_overlapping_rows(miniprot_df, bakta_df, base_name)
    case1_keep, case1_remove = get_case1(case1, bakta_df_bad, args.cds)

    create_and_save_summary(
        base_name,
        [case_2, case3, case4, case1, rv_missing, bakta_df, case1_remove, case1_keep],
        ["case2", "case3", "case4", "case1_initial", "rv_missing", "total_cds", "case1_removed", "case1_kept"]
    )

    miniprot_case4_overlaps = get_overlap_tree(miniprot_df, case4, 'Prot Id', 'Locus Tag')
    case4_miniprot_overlaps = get_overlap_tree(case4, miniprot_df, 'Locus Tag', 'Prot Id')

    single_overlaps = {k: v for k, v in miniprot_case4_overlaps.items() if v['count'] == 1}
    multiple_overlaps = {k: v for k, v in miniprot_case4_overlaps.items() if v['count'] > 1}

    single_overlap_groups = create_overlap_groups(single_overlaps, case4, miniprot_df)
    multiple_overlap_groups = create_overlap_groups(multiple_overlaps, case4, miniprot_df)

    generate_plot(single_overlap_groups, f"{base_name}_single_overlaps.html")
    generate_plot(multiple_overlap_groups, f"{base_name}_multiple_overlaps.html")

    def create_coords_dict(df, id_col):
        return {(row['Start'], row['End']): row[id_col] for _, row in df.iterrows()}

    coords_case2 = create_coords_dict(case_2, 'Locus Tag')
    coords_case3 = create_coords_dict(case3, 'Locus Tag')
    coords_rv_missing = create_coords_dict(rv_missing, 'Prot Id')
    coords_case1_kept = create_coords_dict(case1_keep, 'Locus Tag')

    combined_coordinates_dict = {
        **coords_case2, **coords_case3, **coords_rv_missing, **coords_case1_kept
    }

    def get_overlap_percentage(range1, range2):
        max_start = max(range1[0], range2[0])
        min_end = min(range1[1], range2[1])
        intersection = max(0, min_end - max_start)
        len_range1 = range1[1] - range1[0]
        len_range2 = range2[1] - range2[0]
        max_length = max(len_range1, len_range2)
        return (intersection / max_length) * 100 if max_length > 0 else 0

    b_overlap_groups = create_overlap_groups(case4_miniprot_overlaps, case4, miniprot_df)
    new_groups = {v['Prot Id']: v for k, v in b_overlap_groups.items()}

    for key, value in single_overlap_groups.items():
        bakta_match_nt = value['case4'][0][1] - value['case4'][0][0]
        rv_match_nt = miniprot_df.loc[miniprot_df['Prot Id'] == value['Prot Id']].iloc[0]['Match Nt']

        locus_tag = value['Locus Tags'][0]
        overlaps = 0
        if locus_tag in new_groups:
             overlaps = len(new_groups[locus_tag]['Locus Tags'])

        new_coord, new_id = None, None
        if overlaps > 1:
            new_coord = value['miniprot_df'][0]
            new_id = value['Prot Id']
        elif rv_match_nt > bakta_match_nt:
            new_coord = value['miniprot_df'][0]
            new_id = value['Prot Id']
        elif rv_match_nt < bakta_match_nt and abs(bakta_match_nt - rv_match_nt) > 550:
            new_coord = value['case4'][0]
            new_id = value['Locus Tags'][0]
        else:
            new_coord = value['miniprot_df'][0]
            new_id = value['Prot Id']

        coords_to_add = {new_coord: new_id}
        coords_to_remove = set()

        for exist_coord, exist_id in combined_coordinates_dict.items():
            if new_coord != exist_coord and new_coord[1] > exist_coord[0] and new_coord[0] < exist_coord[1]:
                percentage = get_overlap_percentage(new_coord, exist_coord)
                if percentage >= 70:
                    len_new_coord = new_coord[1] - new_coord[0]
                    len_exist_coord = exist_coord[1] - exist_coord[0]
                    if (len_new_coord > len_exist_coord) and (len_new_coord - len_exist_coord) > 550:
                        coords_to_remove.add(exist_coord)
                    elif (len_exist_coord - len_new_coord) < 550 and "Rv" not in str(exist_id):
                        coords_to_remove.add(exist_coord)
                    else:
                        coords_to_add.clear()
                        break

        for coord in coords_to_remove:
            if coord in combined_coordinates_dict:
                del combined_coordinates_dict[coord]

        combined_coordinates_dict.update(coords_to_add)

    def get_range_length(r):
        return r[1] - r[0]

    for key, value in multiple_overlap_groups.items():
        case4_ranges = value['case4']
        locus_tags = value['Locus Tags']
        miniprot_range = value['miniprot_df'][0]
        prot_id = value['Prot Id']

        best_overlap = (None, None, -1)
        for i, r in enumerate(case4_ranges):
            percentage = get_overlap_percentage(r, miniprot_range)
            if percentage > best_overlap[2]:
                locus_tag = locus_tags[i] if i < len(locus_tags) else "Unknown"
                best_overlap = (r, locus_tag, percentage)

        if best_overlap[0]:
            len_best_overlap = get_range_length(best_overlap[0])
            len_miniprot = get_range_length(miniprot_range)

            if len_miniprot >= len_best_overlap or (len_best_overlap - len_miniprot) < 500:
                largest_range = miniprot_range
                largest_range_id = prot_id
            else:
                largest_range = best_overlap[0]
                largest_range_id = best_overlap[1]

            coords_to_add = {largest_range: largest_range_id}
            coords_to_remove = set()
            new_coord = largest_range

            for exist_coord, exist_id in combined_coordinates_dict.items():
                if new_coord != exist_coord and new_coord[1] > exist_coord[0] and new_coord[0] < exist_coord[1]:
                    percentage = get_overlap_percentage(new_coord, exist_coord)
                    if percentage >= 70:
                        if (get_range_length(new_coord) > get_range_length(exist_coord)) and (get_range_length(new_coord) - get_range_length(exist_coord)) > 550:
                            coords_to_remove.add(exist_coord)
                        else:
                            coords_to_add.clear()
                            break

            for coord in coords_to_remove:
                if coord in combined_coordinates_dict:
                    del combined_coordinates_dict[coord]

            combined_coordinates_dict.update(coords_to_add)

    def get_best_match_filter(group):
        return group.loc[group['Match Nt'].idxmax()]

    all_reference_ids = set(miniprot_df['Prot Id'].unique())
    included_ids = set(combined_coordinates_dict.values())
    missing_ids = all_reference_ids - included_ids

    if missing_ids:
        df_missing_genes = miniprot_df[miniprot_df['Prot Id'].isin(missing_ids)]
        if not df_missing_genes.empty:
            df_best_match_for_missing = df_missing_genes.groupby('Prot Id').apply(get_best_match_filter, include_groups=False).reset_index()
            for _, row in df_best_match_for_missing.iterrows():
                coords_tuple = (row['Start'], row['End'])
                prot_id = row['Prot Id']
                if coords_tuple not in combined_coordinates_dict:
                    combined_coordinates_dict[coords_tuple] = prot_id

    all_values = list(combined_coordinates_dict.values())
    duplicates = {x for x in all_values if all_values.count(x) > 1}

    coords_by_value = {val: [k for k, v in combined_coordinates_dict.items() if v == val] for val in duplicates}

    def get_subset_overlap_percentage(range1, range2):
        max_start = max(range1[0], range2[0])
        min_end = min(range1[1], range2[1])
        intersection = max(0, min_end - max_start)
        min_length = min(get_range_length(range1), get_range_length(range2))
        return (intersection / min_length) * 100 if min_length > 0 else 0

    coords_to_delete_final = set()
    for dup_val, dup_coords in coords_by_value.items():
        for dup_coord in dup_coords:
            for coord, value in combined_coordinates_dict.items():
                if value not in duplicates:
                    if get_subset_overlap_percentage(dup_coord, coord) == 100:
                        coords_to_delete_final.add(dup_coord)

    for coord in coords_to_delete_final:
        if coord in combined_coordinates_dict:
            del combined_coordinates_dict[coord]

    sorted_dict = dict(sorted(combined_coordinates_dict.items(), key=lambda item: item[0]))

    # MODIFICATION: Filter out entries where start or end is 0
    filtered_sorted_dict = {
        (start, end): value
        for (start, end), value in sorted_dict.items()
        if start != 0 and end != 0
    }

    final_list = []
    dup_counter = 1

    temp_list_single = []
    temp_list_multi = []

    # MODIFICATION: Iterate over the newly filtered dictionary
    for (start, end), value in filtered_sorted_dict.items():
        miniprot_matches = miniprot_df[(miniprot_df['Start'] == start) & (miniprot_df['End'] == end)]
        bakta_matches = bakta_df[(bakta_df['Start'] == start) & (bakta_df['End'] == end)]

        strand = ""
        if not miniprot_matches.empty:
            strand = miniprot_matches.iloc[0]['Strand']
        elif not bakta_matches.empty:
            strand = bakta_matches.iloc[0]['Strand']

        if not bakta_matches.empty:
            original_id = bakta_matches.iloc[0]['Locus Tag']
            # Apply formatting to the final ID
            cluster_id = format_final_id(original_id)
            rv_ids = ", ".join(miniprot_matches['Prot Id'].unique())
            temp_list_single.append([cluster_id, strand, start, end, rv_ids])
        else:
            rv_ids = miniprot_matches['Prot Id'].unique()
            joined_rvs = ", ".join(rv_ids)
            if len(rv_ids) == 0:
                # Apply formatting to the final ID
                cluster_id = format_final_id(value)
                temp_list_single.append([cluster_id, strand, start, end, ""])
            elif len(rv_ids) < 2:
                # Apply formatting to the final ID
                cluster_id = format_final_id(f"{base_name}_{rv_ids[0]}")
                temp_list_single.append([cluster_id, strand, start, end, joined_rvs])
            else:
                cluster_id = f"{base_name}_dup_{dup_counter}"
                dup_counter += 1
                temp_list_multi.append([cluster_id, strand, start, end, joined_rvs])

    final_list = temp_list_single + temp_list_multi

    with open(f"{base_name}_corrected.tsv", 'w') as output_file:
        for line_data in sorted(final_list, key=lambda x: x[2]):
            output_file.write('\t'.join(map(str, line_data)) + '\n')

    print(f"Processing complete. Final gene list saved to {base_name}_corrected.tsv")
    print(f"Total genes in final set: {len(final_list)}")

if __name__ == "__main__":
    main()
