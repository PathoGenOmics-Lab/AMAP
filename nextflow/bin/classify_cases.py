#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import os 

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

    filtered_miniprot = last_miniprot.groupby('Prot Id').apply(custom_filter).reset_index(drop=True)
    
    return filtered_miniprot

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

def main():
    parser = argparse.ArgumentParser(description="Classifies Bakta and BLAST records and generates the respective fasta files.")
    parser.add_argument('-b', '--bakta', help="Path to the Bakta TSV file containing genomes to process.")
    parser.add_argument('-p', '--miniprot', help="Path to the miniprot output")
    args = parser.parse_args()

    bakta_df = read_to_pandaBakta(args.bakta)
    miniprot_df = read_to_pandaMiniprot(args.miniprot)

    base_name = os.path.basename(args.bakta).split('.')[0]

    case_2, case3 = process_case2_and_case3(bakta_df, miniprot_df, base_name="your_base_name")
    case4, case1 = process_cases(bakta_df, miniprot_df, base_name)
    no_solapan_con_bakta = get_non_overlapping_rows(miniprot_df, bakta_df, base_name)
    create_and_save_summary(base_name, [case_2, case3, case4, case1, no_solapan_con_bakta, bakta_df], ["case2", "case3", "case4", "case1", "missingRvs", "CDS"])

if __name__ == "__main__":
    main()