#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import os 
from intervaltree import Interval, IntervalTree
import plotly.graph_objects as go

def read_to_pandaBakta(bakta):
    """
    Read and process a genomic data file in TSV format.

    Args:
    bakta (str): File path to the TSV file.

    Process:
    - Reads the file, ignoring lines starting with '#'.
    - Assigns column names for easier data handling.
    - Filters for rows where 'Type' is 'cds'.
    - Calculates the length of coding sequences.

    Returns:
    pandas.DataFrame: Contains data for 'Locus Tag', 'Start', 'End', 'Strand', and coding sequence length.
    """
    bakta_df = pd.read_csv(bakta, comment='#', sep='\t', header=None)
    bakta_df.columns = ["Sequence Id", "Type", "Start", "End", "Strand", "Locus Tag", "Gene", "Product", "DbXrefs"]
    bakta_cds = bakta_df.loc[bakta_df['Type'] == 'cds']
    bakta_cds['Length_nt'] = bakta_cds['End'] - bakta_cds['Start']
    last_bakta = bakta_cds[["Locus Tag", "Start", "End", "Strand", 'Length_nt']]
    return last_bakta

def custom_filter(group):
    """
    Apply a custom filter to a pandas group.

    Process:
    - Checks if any row in the group has a 'Num nt' - 'Match Nt' difference equal to 0.
    - If such rows exist, filters and retains only those rows.
    - If no such rows exist, returns the entire group unaltered.

    Args:
    group (pandas.DataFrame): A DataFrame group from a GroupBy operation.

    Returns:
    pandas.DataFrame: The filtered or unaltered group based on the condition.
    """
    # Check if any row meets the condition ('Num nt' - 'Match Nt' == 0)
    if any(group['Num nt'] - group['Match Nt'] == 0):
        # Return only rows where the condition is met
        return group[group['Num nt'] - group['Match Nt'] == 0]
    else:
        # Return the entire group unaltered
        return group

def read_to_pandaMiniprot(miniprot):
    """
    Read and process a miniprot data file in TSV format.

    Args:
    miniprot (str): File path to the miniprot TSV file.

    Process:
    - Reads the file into a DataFrame, with tab-separated values and no header.
    - Assigns specific column names for ease of data handling.
    - Adjusts the 'Start' column by adding 1 to each value.
    - Selects a subset of columns for the final DataFrame.
    - Applies a custom filter function on groups of data, grouped by 'Prot Id'.

    Returns:
    pandas.DataFrame: Contains processed and filtered data, with columns such as 'Prot Id', 'Prot len', 'Strand', etc.
    """
    miniprot_df = pd.read_csv(miniprot, sep='\t', header=None)
    miniprot_df.columns = ["Prot Id", "Prot len", "Prot Start", "ProtStop", "Strand", "fasta Tag", "fasta len", "Start", "End", "Match Nt", "Num nt", "mapping qual", "AS","MS","NP","FS","ST","DA","DO","CG","CS"]
    miniprot_df['Start'] = miniprot_df['Start'] + 1
    last_miniprot = miniprot_df[["Prot Id", "Prot len", "Strand", "Start", "End", "Match Nt", "Num nt"]]
    #filtered_miniprot = last_miniprot
    filtered_miniprot = last_miniprot.groupby('Prot Id').apply(custom_filter).reset_index(drop=True)
    
    return filtered_miniprot

def read_to_pandaBadbakta(miniprot):
    """
    Read and process a miniprot data file in TSV format.

    Args:
    miniprot (str): File path to the miniprot TSV file.

    Process:
    - Reads the file into a DataFrame, with tab-separated values and no header.
    - Assigns specific column names for ease of data handling.
    - Adjusts the 'Start' column by adding 1 to each value.
    - Selects a subset of columns for the final DataFrame.
    - Applies a custom filter function on groups of data, grouped by 'Prot Id'.

    Returns:
    pandas.DataFrame: Contains processed and filtered data, with columns such as 'Prot Id', 'Prot len', 'Strand', etc.
    """
    miniprot_df = pd.read_csv(miniprot, sep='\t', header=None)
    miniprot_df.columns = ["Prot Id", "Prot len", "Prot Start", "ProtStop", "Strand", "fasta Tag", "fasta len", "Start", "End", "Match Nt", "Num nt", "mapping qual", "AS","MS","NP","FS","ST","DA","DO","CG","CS"]
    miniprot_df['Start'] = miniprot_df['Start'] + 1
    last_miniprot = miniprot_df[["Prot Id", "Prot len", "Strand", "Start", "End", "Match Nt", "Num nt"]]
    
    return last_miniprot

def process_case2_and_case3(bakta_df, miniprot_df, base_name=None):
    """
    Merge two dataframes and separate into case_2 and case3 based on a calculated difference.

    Args:
    bakta_df (pandas.DataFrame): DataFrame containing 'Start', 'End', and 'Length_nt' columns.
    miniprot_df (pandas.DataFrame): DataFrame containing 'Start', 'End', and 'Num nt' columns.
    base_name (str, optional): Base name for saving the case data to TSV files. Default is None.

    Returns:
    tuple: Two DataFrames corresponding to case_2 and case3.
    """
    # Merge DataFrames on 'Start' and 'End'
    case2 = pd.merge(bakta_df, miniprot_df, on=['Start', 'End'])

    case2['Diferencia'] = (case2['Length_nt'] - case2['Num nt']) + 1
    # Create case_2 where 'Diferencia' is 0
    case_2 = case2[case2['Diferencia'] == 0]

    # Create case3 where 'Diferencia' is not 0
    case3 = case2[case2['Diferencia'] != 0]
    
    case_2 = case_2[["Locus Tag", "Start", "End", "Strand_x", 'Length_nt']].drop_duplicates()
    # Optionally save the results
    if base_name:
        case3.to_csv(base_name + "-case3.tsv", sep='\t', index=False)
        case_2.to_csv(base_name + "-case2.tsv", sep='\t', index=False)

    return case_2, case3

def process_cases(bakta_df, miniprot_df, base_name=None):
    """
    Process and separate data into case1 and case4 based on overlap conditions.

    Args:
    bakta_df (pandas.DataFrame): DataFrame to be processed.
    miniprot_df (pandas.DataFrame): DataFrame to compare against for overlaps.
    base_name (str, optional): Base name for saving the case data to TSV files. Default is None.

    Returns:
    tuple: Two DataFrames corresponding to case4 and case1.
    """
    # Calculate length for bakta_df
    bakta_df['Length_nt'] = bakta_df['End'] - bakta_df['Start']

    # Merge DataFrames on common columns
    case2 = pd.merge(bakta_df, miniprot_df, on=['Start', 'End'], how='inner')

    # Identify rows not in case2
    common_columns = [col for col in bakta_df.columns if col in case2.columns]
    merged_indicator = bakta_df.merge(case2, on=common_columns, how='left', indicator=True)
    rows_not_in_case2 = merged_indicator[merged_indicator['_merge'] == 'left_only'].drop('_merge', axis=1)
    rows_not_in_case2 = rows_not_in_case2[["Locus Tag", "Start", "End", "Strand", 'Length_nt']]

    # Check for overlaps
    rows_not_in_case2['Solapa'] = False
    for index, row in rows_not_in_case2.iterrows():
        solapa = ((miniprot_df['Start'] <= row['End']) & (miniprot_df['End'] >= row['Start'])).any()
        rows_not_in_case2.at[index, 'Solapa'] = solapa

    # Separate into case4 and case1
    case4 = rows_not_in_case2[rows_not_in_case2['Solapa'] == True]
    case1 = rows_not_in_case2[rows_not_in_case2['Solapa'] == False]
    case1 = case1[["Locus Tag", "Start", "End", "Strand", 'Length_nt']]
    case4 = case4[["Locus Tag", "Start", "End", "Strand", 'Length_nt']]

    # Optionally save the results
    if base_name:
        case4.to_csv(base_name + "-case4.tsv", sep='\t', index=False)
        case1.to_csv(base_name + "-case1.tsv", sep='\t', index=False)

    return case4, case1


def get_non_overlapping_rows(miniprot_df, bakta_df, base_name=None):
    """
    Identifies rows in miniprot_df that do not overlap with any row in bakta_df.

    Args:
    miniprot_df (pandas.DataFrame): DataFrame containing 'Start' and 'End' columns to compare.
    bakta_df (pandas.DataFrame): DataFrame containing 'Start' and 'End' columns to compare against.
    base_name (str, optional): Base name for saving the non-overlapping rows to a TSV file. Default is None.

    Returns:
    no_solapan_con_bakta (pandas.DataFrame): DataFrame containing rows from miniprot_df that do not overlap with bakta_df.
    """
    # Initialize a column to flag overlap
    miniprot_df['Solapa'] = False

    # Check for overlap for each row
    for index, row in miniprot_df.iterrows():
        solapa = ((bakta_df['Start'] <= row['End']) & (bakta_df['End'] >= row['Start'])).any()
        miniprot_df.at[index, 'Solapa'] = solapa

    # Select rows that do not overlap
    no_solapan_con_bakta = miniprot_df[miniprot_df['Solapa'] == False]
    no_solapan_con_bakta = no_solapan_con_bakta[~((no_solapan_con_bakta['Match Nt'] == 0) & (no_solapan_con_bakta['Num nt'] == 0))]
    # Optionally save the result
    if base_name:
        no_solapan_con_bakta.to_csv(base_name + "-rvmissing.tsv", sep='\t', index=False)

    return no_solapan_con_bakta

def create_and_save_summary(base_name, data_frames, case_types):
    """
    Creates and saves a summary DataFrame as a TSV file based on provided data frames and case types.

    Args:
    base_name (str): Base name used for identifying the sequence and lineage, and as part of the output file name.
    data_frames (list of pandas.DataFrame): List of data frames representing different cases or datasets.
    case_types (list of str): List of case type names corresponding to each data frame.

    Process:
    - Constructs a summary dictionary with case types, their counts, the base sequence name, and the lineage.
    - Converts this dictionary into a pandas DataFrame.
    - Saves the DataFrame to a TSV file named based on the base_name.

    Returns:
    None: The function doesn't return anything. It saves the summary DataFrame to a file.
    """
    summary_data = {
        "Type": case_types,
        "Number": [len(df) for df in data_frames],
        "Sequence": base_name,
        "Lineage": base_name.split("_")[0]
    }
    summary_df = pd.DataFrame(summary_data)

    # Save the summary DataFrame to a TSV file
    summary_df.to_csv(base_name + "-summary.tsv", sep='\t', index=False)

def calculate_overlap_length(row1, df2):
    """
    Calculate the total length of overlaps between row1 in df1 and all rows in df2.
    """
    def overlap_length(start1, end1, start2, end2):
        return max(0, min(end1, end2) - max(start1, start2))
    total_overlap_length = df2.apply(lambda row2: overlap_length(row1['Start'], row1['End'], row2['Start'], row2['End']), axis=1)

    return total_overlap_length.sum()


def get_case1(case1, bakta_df_bad, cds):
    
    case1['Overlap_Length'] = case1.apply(lambda row: calculate_overlap_length(row, bakta_df_bad), axis=1)
    case1['Overlap_Percentage'] = (case1['Overlap_Length'] / case1['Length_nt']) * 100
    # Filter for non-exact overlaps
    # Filter for Overlap_Percentage >= 80
    case1_remove = case1[case1['Overlap_Percentage'] >= 80]
    # Filter for Overlap_Percentage < 80
    case1_check = case1[case1['Overlap_Percentage'] < 80]
    locus_toref = read_to_pandaBadbakta(cds)
    locus_tags = case1_check['Locus Tag'].tolist()
    # Filter df_other for rows where 'Prot Id' matches any value in locus_tags
    matching_rows = locus_toref[locus_toref['Prot Id'].isin(locus_tags)]
    matching_rows['match_per'] = 100 - ( (matching_rows['Num nt'] - matching_rows['Match Nt']) / matching_rows['Num nt']) * 100
    matching_rows['match_per'] = matching_rows['match_per'].fillna(0)
    case1_remove_too = matching_rows[matching_rows['match_per'] >= 80]
    case1_keep = matching_rows[matching_rows['match_per'] < 80]
    locus_case1_remove_too = case1_remove_too['Prot Id'].tolist()
    locus_case1_keep = case1_keep['Prot Id'].tolist()
    case1_remove2 = case1_check[case1_check['Locus Tag'].isin(locus_case1_remove_too)]
    case1_keep2 = case1_check[case1_check['Locus Tag'].isin(locus_case1_keep)]
    final_case1_remove = pd.concat([case1_remove, case1_remove2])
    return case1_keep2, final_case1_remove

def grupos_solapamientos(filtered_overlaps,case4, miniprot_df):
    overlaps_df = pd.DataFrame.from_dict(filtered_overlaps, orient='index')

    # Crear grupos de solapamientos
    grupos_solapamiento = {}
    for i, row in overlaps_df.iterrows():
    # Índices de case4 y miniprot_df que pertenecen a este grupo
        indices_case4_grupo = case4[case4['Locus Tag'].isin(row['Locus Tags'])].index
        indices_miniprot_df_grupo = [i]  # El índice de miniprot_df ya está dado

    # Rangos de case4 y miniprot_df para este grupo
        rangos_case4 = [(case4.at[idx, 'Start'], case4.at[idx, 'End']) for idx in indices_case4_grupo]
        rangos_miniprot_df = [(miniprot_df.at[idx, 'Start'], miniprot_df.at[idx, 'End']) for idx in indices_miniprot_df_grupo]

    # Asignar a la estructura de grupos
        #grupos_solapamiento[i] = {'case4': rangos_case4, 'miniprot_df': rangos_miniprot_df}
        prot_id = overlaps_df.at[i, 'Prot Id']
        locus_tags = row['Locus Tags']
        grupos_solapamiento[i] = {
        'case4': rangos_case4, 
        'miniprot_df': rangos_miniprot_df,
        'Prot Id': prot_id, 
        'Locus Tags': locus_tags
        }
    return grupos_solapamiento

def grupos_solapamientosMOD(filtered_overlaps, miniprot_df, case4):
    grupos_solapamiento = {}

    for idx, overlap_data in filtered_overlaps.items():
        prot_id = overlap_data['Prot Id']
        # Usar un conjunto (set) para evitar duplicados
        locus_tags_solapados = set(overlap_data['Locus Tags'])

        # Rangos para los Locus Tags en miniprot_df
        rangos_miniprot_df = []
        for tag in locus_tags_solapados:
            filtered_miniprot = miniprot_df[miniprot_df['Prot Id'] == tag]
            for _, row in filtered_miniprot.iterrows():
                rango = (row['Start'], row['End'])
                if rango not in rangos_miniprot_df:
                    rangos_miniprot_df.append(rango)

        # Rangos para el Prot Id en case4
        rangos_case4 = []
        filtered_case4 = case4[case4['Locus Tag'] == prot_id]
        for _, row in filtered_case4.iterrows():
            rango_case4 = (row['Start'], row['End'])
            if rango_case4 not in rangos_case4:
                rangos_case4.append(rango_case4)

        # Crear entrada en el diccionario de solapamientos
        if rangos_miniprot_df or rangos_case4:
            grupos_solapamiento[idx] = {
                'Prot Id': prot_id,
                'Locus Tags': list(locus_tags_solapados),
                'miniprot_df': rangos_miniprot_df,
                'case4': rangos_case4
            }

    return grupos_solapamiento

def generate_plot(sgrupos_solapamiento, salida):
    fig = go.Figure()

    grosor_rojo = 10  # Grosor para miniprot_df
    grosor_azul = 6  # Grosor para case4

    for i, (clave, grupo) in enumerate(sgrupos_solapamiento.items()):
        prot_id = grupo['Prot Id']
        locus_tags = grupo['Locus Tags']

        case4_x, case4_y, case4_text = [], [], []
        for idx, (start, end) in enumerate(grupo['case4']):
            case4_x.extend([start, end, None])
            case4_y.extend([i, i, None])
            case4_text.extend([locus_tags[idx], locus_tags[idx], None])
        
        fig.add_trace(go.Scatter(
        x=case4_x,
        y=case4_y,
        mode='lines+markers',
        line=dict(color='blue', width=grosor_azul),
        marker=dict(color='blue', size=14),  # Tamaño reducido de los marcadores
        text=case4_text,
        hoverinfo='text',
        opacity=0.5,  # Aumento de la transparencia
        name="bakta-" + str(case4_y[0])
        ))

        miniprot_x, miniprot_y, miniprot_text = [], [], []
        
        for start, end in grupo['miniprot_df']:
            miniprot_x.extend([start, end, None])
            miniprot_y.extend([i, i, None])
            miniprot_text.extend([prot_id, prot_id, None])
        fig.add_trace(go.Scatter(
            x=miniprot_x,
            y=miniprot_y,
            mode='lines+markers',
            line=dict(color='red', width=grosor_rojo),
            marker=dict(color='red', size=14),  # Tamaño reducido de los marcadores
            text=miniprot_text,
            hoverinfo='text',
            opacity=0.4,
            name=miniprot_text[0]
        ))

    fig.update_layout(
    title='Visualización de Grupos de Solapamientos con Marcadores de Extremos',
    xaxis_title='Rango',
    yaxis_title='Grupo de Solapamiento',
    hovermode='closest'
    )

    fig.write_html(salida)

def get_overlaptree(miniprot_df, case4):
    """
    Builds an Interval Tree from 'miniprot_df' and counts overlaps with intervals in 'case4'.

    Parameters:
    - miniprot_df (DataFrame): Contains 'Start', 'End', and 'Prot Id' columns for intervals.
    - case4 (DataFrame): Contains 'Start', 'End', and 'Locus Tag' columns to check for overlaps.

    Returns:
    - overlap_counts (dict): A dictionary mapping indices from 'miniprot_df' to a dict with:
        - 'count': Number of overlaps.
        - 'Prot Id': Associated 'Prot Id'.
        - 'Locus Tags': List of overlapping 'Locus Tags' from 'case4'.

    Notes:
    - Invalid intervals (start > end) in 'miniprot_df' are reported as errors.
    - Assumes numeric 'Start' and 'End' values.
    """
    tree = IntervalTree()
    for idx, row in miniprot_df.iterrows():
        start = row['Start']
        end = row['End']
        if start <= end:
            tree[start:end] = (idx, row['Prot Id'])
        else:
            print(f"Intervalo inválido: Start ({start}) es mayor que End ({end})")

    overlap_counts = {}
    for idx, row in case4.iterrows():
        overlaps = tree[row['Start']:row['End']]
        for overlap in overlaps:
            prot_idx, prot_id = overlap.data
            if prot_idx not in overlap_counts:
                overlap_counts[prot_idx] = {'count': 1, 'Prot Id': prot_id, 'Locus Tags': [row['Locus Tag']]}
            else:
                overlap_counts[prot_idx]['count'] += 1
                overlap_counts[prot_idx]['Locus Tags'].append(row['Locus Tag'])
    return overlap_counts

def get_overlaptreeMOD(tree_df, overlap_df, tree_id_col, overlap_id_col):
    """
    Builds an Interval Tree using intervals from 'tree_df' and counts overlaps with intervals in 'overlap_df'.

    Parameters:
    - tree_df (DataFrame): DataFrame to build the Interval Tree. Must contain 'Start' and 'End' columns.
    - overlap_df (DataFrame): DataFrame to check for overlaps. Must contain 'Start' and 'End' columns.
    - tree_id_col (str): Column name in 'tree_df' for IDs.
    - overlap_id_col (str): Column name in 'overlap_df' for IDs.

    Returns:
    - overlap_counts (dict): A dictionary mapping indices from 'tree_df' to a dict with:
        - 'count': Number of overlaps.
        - 'Prot Id': ID from 'tree_df'.
        - 'Locus Tags': List of overlapping IDs from 'overlap_df'.

    Notes:
    - Invalid intervals (start > end) in 'tree_df' are reported as errors.
    - Assumes numeric 'Start' and 'End' values.
    """

    tree = IntervalTree()
    for idx, row in tree_df.iterrows():
        start = row['Start']
        end = row['End']
        if start <= end:
            tree[start:end] = (idx, row[tree_id_col])
        else:
            print(f"Invalid interval: Start ({start}) is greater than End ({end})")

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

def main():
    parser = argparse.ArgumentParser(description="Classifies Bakta and BLAST records and generates the respective fasta files.")
    parser.add_argument('-b', '--bakta', help="Path to the Bakta TSV file containing genomes to process.")
    parser.add_argument('-p', '--miniprot', help="Path to the miniprot output")
    parser.add_argument('-bb', '--badbakta', help="Path to the miniprot badbakta output")
    parser.add_argument('-c', '--cds', help="Path to the bakta to reference")
    args = parser.parse_args()

    bakta_df = read_to_pandaBakta(args.bakta)
    miniprot_df = read_to_pandaMiniprot(args.miniprot)
    
    base_name = os.path.basename(args.bakta).split('.')[0]

    case_2, case3 = process_case2_and_case3(bakta_df, miniprot_df, base_name)
    case4, case1 = process_cases(bakta_df, miniprot_df, base_name)
    no_solapan_con_bakta = get_non_overlapping_rows(miniprot_df, bakta_df, base_name)
    bakta_df_bad = read_to_pandaBadbakta(args.badbakta)

    # Assuming 'case1' and 'bakta_df_bad' are your DataFrames
    # Calculate overlap length for all rows in 'case1'
    case1_keep2, final_case1_remove = get_case1(case1, bakta_df_bad, args.cds)
    
    create_and_save_summary(base_name, [case_2, case3, case4, case1, no_solapan_con_bakta, bakta_df,final_case1_remove, case1_keep2], ["case2", "case3", "case4", "case1", "missingRvs", "CDS","case1_remove","case1_keep"])
# ____________________________________________

    overlap_counts = get_overlaptree(miniprot_df, case4)
    overlap_counts = get_overlaptreeMOD(miniprot_df, case4, 'Prot Id', 'Locus Tag')
    
    #b_overlap_counts = get_overlaptreeMOD(case4, miniprot_df)
    b_overlap_counts = get_overlaptreeMOD(case4, miniprot_df, 'Locus Tag', 'Prot Id')

    # Filtrar para obtener solo aquellos que solapan más de una vez
    filtered_overlaps = {k: v for k, v in overlap_counts.items() if v['count'] > 1}
    single_overlaps = {k: v for k, v in overlap_counts.items() if v['count'] == 1}
    grupos_solapamiento = grupos_solapamientos(filtered_overlaps, case4, miniprot_df)
    sgrupos_solapamiento = grupos_solapamientos(single_overlaps, case4, miniprot_df)
    
    grupos = grupos_solapamientosMOD(b_overlap_counts, miniprot_df, case4)


    # Creando un nuevo diccionario con 'Prot Id' como clave
    new_groups = {}
    for clave, valor in grupos.items():
        prot_id = valor['Prot Id']
        new_groups[prot_id] = valor

    generate_plot(sgrupos_solapamiento, base_name + "_single_overlap.html")
    generate_plot(grupos_solapamiento, base_name + "_multiple_overlap.html")

#____________________________________________________________________________________
    #coordinates_case2 = set(zip(case_2['Start'], case_2['End']))
    #coordinates_case3 = set(zip(case3['Start'], case3['End']))
    #coordinates_rv = set(zip(no_solapan_con_bakta['Start'], no_solapan_con_bakta['End']))
    #coordinates_case1 = set(zip(case1_keep2['Start'], case1_keep2['End']))
    #
    #coordenadas_ordenadas = sorted (coordinates_case2 | coordinates_case3 | coordinates_rv | coordinates_case1)

    def create_coordinates_dict(df):
        return {(row['Start'], row['End']): row['Locus Tag'] for index, row in df.iterrows()}

    def create_coordinates_dict2(df):
        return {(row['Start'], row['End']): row['Prot Id'] for index, row in df.iterrows()}    

    # Creación de los diccionarios individuales
    coordinates_dict_case2 = create_coordinates_dict(case_2)
    coordinates_dict_case3 = create_coordinates_dict(case3)
    coordinates_dict_rv = create_coordinates_dict2(no_solapan_con_bakta)
    coordinates_dict_case1 = create_coordinates_dict(case1_keep2)

    # Combinar los cuatro diccionarios en uno
    combined_coordinates_dict = {**coordinates_dict_case2, **coordinates_dict_case3, **coordinates_dict_rv, **coordinates_dict_case1}

    def calcular_porcentaje_solapamiento(rango1, rango2):
        inicio_max = max(rango1[0], rango2[0])
        fin_min = min(rango1[1], rango2[1])
        interseccion = max(0, fin_min - inicio_max)
        longitud_rango1 = rango1[1] - rango1[0]
        longitud_rango2 = rango2[1] - rango2[0]
        longitud_mayor = max(longitud_rango1, longitud_rango2)  # Usar la longitud mayor
        if longitud_mayor == 0:
            return 0
        return (interseccion / longitud_mayor) * 100

    for key, value in sgrupos_solapamiento.items():
        bakta_match_nt = value['case4'][0][1] - value['case4'][0][0]
        rv_match_nt = miniprot_df.loc[miniprot_df['Prot Id'] == value['Prot Id']].iloc[0]['Match Nt']
        #print(rv_match_nt)
        #nueva_coord = value['miniprot_df'][0] if rv_match_nt > bakta_match_nt else value['case4'][0]
        #nuevo_id = value['Prot Id'] if rv_match_nt > bakta_match_nt else value['Locus Tags'][0]
        overlaps = len(new_groups[value['Locus Tags'][0]]['Locus Tags'])
        if overlaps > 1:
            nueva_coord = value['miniprot_df'][0]
            nuevo_id = value['Prot Id']
        elif rv_match_nt > bakta_match_nt :
            nueva_coord = value['miniprot_df'][0]
            nuevo_id = value['Prot Id']
        elif rv_match_nt < bakta_match_nt and abs(bakta_match_nt - rv_match_nt) > 550 :
            nueva_coord = value['case4'][0]
            nuevo_id = value['Locus Tags'][0]
        else:
            nueva_coord = value['miniprot_df'][0]
            nuevo_id = value['Prot Id']      

        coordenadas_para_agregar = {nueva_coord: nuevo_id}
        coordenadas_para_eliminar = set()

        for exist_coord, exist_id in combined_coordinates_dict.items():
            
            if nueva_coord != exist_coord and nueva_coord[1] > exist_coord[0] and nueva_coord[0] < exist_coord[1]:
                porcentaje = calcular_porcentaje_solapamiento(nueva_coord, exist_coord)

                if porcentaje >= 70:
                    len_nueva_coord = nueva_coord[1] - nueva_coord[0]
                    len_exist_coord = exist_coord[1] - exist_coord[0]
                    diff = len_exist_coord - len_nueva_coord
                    rev_diff = len_nueva_coord - len_exist_coord

                    if (nueva_coord[1] - nueva_coord[0]) > (exist_coord[1] - exist_coord[0]) and rev_diff > 550:

                        coordenadas_para_eliminar.add(exist_coord)
                    elif (exist_coord[1] - exist_coord[0]) - (nueva_coord[1] - nueva_coord[0]) < 550 and "Rv" not in exist_id: #new

                        coordenadas_para_eliminar.add(exist_coord) #new

                    else:
                        coordenadas_para_agregar.clear()

                        break
          
        for coord in coordenadas_para_eliminar:
            del combined_coordinates_dict[coord]
        

        combined_coordinates_dict.update(coordenadas_para_agregar)

    def longitud_rango(rango):
        return rango[1] - rango[0]

    for key, value in grupos_solapamiento.items():
        rangos_case4 = value['case4']
        locus_tags = value['Locus Tags']
        rango_miniprot_df = value['miniprot_df'][0]
        prot_id = value['Prot Id']


        mejor_solapamiento = (None, None, -1)  # (rango, locus tag, porcentaje)

        for i, rango in enumerate(rangos_case4):
            porcentaje = calcular_porcentaje_solapamiento(rango, rango_miniprot_df)
            if porcentaje > mejor_solapamiento[2]:
                locus_tag = locus_tags[i] if i < len(locus_tags) else "Desconocido"
                mejor_solapamiento = (rango, locus_tag, porcentaje)

        if mejor_solapamiento[0]:
            longitud_mejor_solapamiento = longitud_rango(mejor_solapamiento[0])
            longitud_miniprot = longitud_rango(rango_miniprot_df)
            diff_long = longitud_mejor_solapamiento - longitud_miniprot

                
            if longitud_miniprot >= longitud_mejor_solapamiento:
                rango_mas_grande = rango_miniprot_df
                identificador_rango_mas_grande = prot_id
            elif  diff_long < 500 :
                rango_mas_grande = rango_miniprot_df
                identificador_rango_mas_grande = prot_id
            else:
                rango_mas_grande = mejor_solapamiento[0]
                identificador_rango_mas_grande = mejor_solapamiento[1]  # Solo el locus tag

            #print(f"Rango más grande: {rango_mas_grande}, Identificador: {identificador_rango_mas_grande}")
            
            coordenadas_para_agregar = {rango_mas_grande: identificador_rango_mas_grande}
            coordenadas_para_eliminar = set()
            nueva_coord = rango_mas_grande

            for exist_coord, exist_id in combined_coordinates_dict.items():
                if nueva_coord != exist_coord and nueva_coord[1] > exist_coord[0] and nueva_coord[0] < exist_coord[1]:
                    porcentaje = calcular_porcentaje_solapamiento(nueva_coord, exist_coord)
                    if porcentaje >= 70:
                        len_nueva_coord = nueva_coord[1] - nueva_coord[0]
                        len_exist_coord = exist_coord[1] - exist_coord[0]
                        diff = len_exist_coord - len_nueva_coord
                        rev_diff = len_nueva_coord - len_exist_coord
                        if (nueva_coord[1] - nueva_coord[0]) > (exist_coord[1] - exist_coord[0]) and rev_diff > 550:
                            coordenadas_para_eliminar.add(exist_coord)
                        else:
                            coordenadas_para_agregar.clear()
                            break
            for coord in coordenadas_para_eliminar:
                del combined_coordinates_dict[coord]

            combined_coordinates_dict.update(coordenadas_para_agregar)
    
    def custom_filter2(group):
        return group.loc[group['Match Nt'].idxmax()]
    
    always_missing = ['Rv0829','Rv0609A', 'Rv0724A', 'Rv1089A', 'Rv1507A', 'Rv2023A', 'Rv2085', 'Rv2160A', 'Rv2274A', 'Rv2369c', 'Rv2408', 'Rv2432c', 'Rv2530A', 'Rv2614A', 'Rv2755c', 'Rv2767c', 'Rv2810c', 'Rv2862A', 'Rv2871', 'Rv2929', 'Rv2970A', 'Rv3566A', 'Rv3724A', 'Rv3922c', 'Rv3374']
    df_filtrado = miniprot_df[miniprot_df['Prot Id'].apply(lambda x: any(id in x for id in always_missing))]
    df_filtrado = df_filtrado.groupby('Prot Id').apply(custom_filter2).reset_index(drop=True)
    for index, row in df_filtrado.iterrows():
        start_end_tuple = (row['Start'], row['End'])
        prot_id = row['Prot Id']
        # Verificar si la tupla no está en el diccionario y añadirla
        if start_end_tuple not in combined_coordinates_dict:
            combined_coordinates_dict[start_end_tuple] = prot_id

    #for key, value in combined_coordinates_dict.items():
    #    if "Rv0493c" in value:
    #        print(key, value)
    #print(combined_coordinates_dict)
    valores = list(combined_coordinates_dict.values())
    duplicados = set([x for x in valores if valores.count(x) > 1])
    from collections import defaultdict
    # Agrupar coordenadas por valores duplicados
    coord_por_valor = defaultdict(list)
    for coord, valor in combined_coordinates_dict.items():
        if valor in duplicados:
            coord_por_valor[valor].append(coord)
    #print(coord_por_valor)
    def calcular_porcentaje_solapamiento2(rango1, rango2):
        inicio_max = max(rango1[0], rango2[0])
        fin_min = min(rango1[1], rango2[1])
        interseccion = max(0, fin_min - inicio_max)

        longitud_rango1 = rango1[1] - rango1[0]
        longitud_rango2 = rango2[1] - rango2[0]
        longitud_menor = min(longitud_rango1, longitud_rango2)

        if longitud_menor == 0:
            return 0
        return (interseccion / longitud_menor) * 100
    lista_dup = []
    for valor_duplicado, coords_duplicadas in coord_por_valor.items():
        for coord_dup in coords_duplicadas:
            for coord, valor in combined_coordinates_dict.items():
            # Evitar comparar la coordenada duplicada con otras del mismo ID duplicado
                if valor not in duplicados:
                    porcentaje = calcular_porcentaje_solapamiento2(coord_dup, coord)
                    if porcentaje == 100:  # Si hay solapamiento del 100%
                        lista_dup.append(coord_dup)
    for coord in lista_dup:
        del combined_coordinates_dict[coord]
    #print(combined_coordinates_dict)
    sorted_dict = dict(sorted(combined_coordinates_dict.items(), key=lambda item: item[0]))
    print(sorted_dict)

if __name__ == "__main__":
    main()