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