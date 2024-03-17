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