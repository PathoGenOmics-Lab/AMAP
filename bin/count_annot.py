#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script reads a master annotation table, calculates the frequency of
eggNOG COG categories and Mycobrowser Functional Categories, saves these
counts to summary TSV files, and generates a bar plot visualizing the
category distributions.
"""

import argparse
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Headless backend — must be set before pyplot import
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Tuple

__author__ = 'Paula Ruiz-Rodriguez'
__credits__ = ['Paula Ruiz-Rodriguez', 'PathoGenOmics']
__license__ = ['GPL 3']
__version__ = '1.0.1'
__maintainer__ = '@paururo'
__email__ = 'paula.ruiz.rodriguez@csic.es'
__status__ = 'finished'
__lab__ = 'PathoGenOmics (I2SysBio)'

# What a gene without a COG category counts as, once. The report reads 'NaN'
# as its category "NA", which it shows as "None".
NO_COG = 'NaN'


def cog_letters(category) -> list:
    """
    The COG categories of a gene, one letter each: "K", or "JK" for a gene in
    both J and K. A gene with none, which eggNOG writes as '-', is NO_COG.
    """
    if pd.isna(category):
        return [NO_COG]
    return sorted(set(re.findall('[A-Z]', str(category).upper()))) or [NO_COG]


def process_cog_categories(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates the count and percentage of each COG category.

    The 'COG_category' column can contain multiple categories (e.g., "JKS").
    This function expands them so that each category is counted individually.
    A gene without a COG category counts once, as NO_COG.
    The percentage is calculated against the total number of genes (rows)
    in the input dataframe.

    Args:
        df: The input DataFrame containing annotation data.

    Returns:
        A DataFrame with COG category counts and percentages.
    """
    total_genes = len(df)
    
    # Create a temporary DataFrame for COG processing
    cog_df = df[['lineage', 'genome', 'COG_category']].copy()
    
    # Each gene's categories. A gene without one used to be filled in as
    # 'Unassigned' and split into letters like a category string, so it
    # counted once in A, D, E, G, I and U and twice in N and S.
    cog_df['COG_category'] = cog_df['COG_category'].apply(cog_letters)
    
    # Create a new row for each COG category character
    cog_df_exploded = cog_df.explode('COG_category')

    # Group by lineage, genome, and category to count occurrences
    cog_counts = cog_df_exploded.groupby(
        ['lineage', 'genome', 'COG_category']
    ).size().reset_index(name='count')

    # Calculate percentage based on the total number of genes.
    # NOTE: Multi-letter COG categories (e.g. "JKS") are exploded to individual
    # letters, so a gene can contribute to multiple categories. This means the
    # sum of percentages across all categories can exceed 100% — this is expected
    # behaviour and matches eggNOG-mapper convention.
    cog_counts['percentage_of_genes'] = (cog_counts['count'] / total_genes) * 100
    
    return cog_counts

def process_mycobrowser_categories(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates the count and percentage of each Mycobrowser Functional Category.

    The percentage is calculated relative to the total number of functionally
    annotated genes within each respective genome.

    Args:
        df: The input DataFrame containing annotation data.

    Returns:
        A DataFrame with Mycobrowser category counts and percentages.
    """
    # Group by lineage, genome, and category to count occurrences
    myco_counts = df.groupby(
        ['lineage', 'genome', 'Functional_Category'], dropna=False
    ).size().reset_index(name='count')
    
    # Replace potential NaN in category column for clarity
    myco_counts['Functional_Category'] = myco_counts['Functional_Category'].fillna('Unassigned')

    # Calculate percentage relative to the total counts within each genome
    total_counts_per_genome = myco_counts.groupby(['lineage', 'genome'])['count'].transform('sum')
    myco_counts['percentage'] = (myco_counts['count'] / total_counts_per_genome) * 100

    return myco_counts

def create_summary_plots(cog_df: pd.DataFrame, myco_df: pd.DataFrame, output_path: Path):
    """
    Generates and saves a figure with two bar plots for COG and Mycobrowser categories.

    Args:
        cog_df: DataFrame with COG category statistics.
        myco_df: DataFrame with Mycobrowser category statistics.
        output_path: The path to save the output PNG image.
    """
    # Aggregate data for plotting (sum counts across all genomes for simplicity)
    cog_plot_data = cog_df.groupby('COG_category')['count'].sum().reset_index()
    myco_plot_data = myco_df.groupby('Functional_Category')['count'].sum().reset_index()

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 15))
    
    # Plot 1: EggNOG COG Categories
    ax1.bar(cog_plot_data['COG_category'], cog_plot_data['count'], color=plt.cm.viridis(np.linspace(0, 1, len(cog_plot_data))))
    ax1.set_title('eggNOG Functional Category Distribution', fontsize=14)
    ax1.set_ylabel('Total Gene Count', fontsize=12)
    ax1.tick_params(axis='x', rotation=45)
    ax1.grid(axis='y', linestyle='--', alpha=0.7)

    # Plot 2: Mycobrowser Functional Categories
    ax2.bar(myco_plot_data['Functional_Category'], myco_plot_data['count'], color=plt.cm.plasma(np.linspace(0, 1, len(myco_plot_data))))
    ax2.set_title('Mycobrowser Functional Category Distribution', fontsize=14)
    ax2.set_ylabel('Total Gene Count', fontsize=12)
    ax2.tick_params(axis='x', labelsize=8)
    ax2.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Improve layout
    plt.setp(ax2.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    fig.tight_layout(pad=3.0)
    
    # Save the figure
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Successfully generated plot at {output_path}")


def main():
    """
    Main function to parse arguments and orchestrate the analysis.
    """
    parser = argparse.ArgumentParser(
        description="Analyze functional categories from a master annotation table and generate summary files and plots."
    )
    parser.add_argument(
        "master_table_file",
        help="Path to the input master annotation TSV file (e.g., *-AMAPmastertable.tsv)."
    )
    parser.add_argument(
        "--sample",
        help="The sample's ID in the sample sheet, for the genome column "
             "(default: the name of the master table)."
    )
    parser.add_argument(
        "--lineage",
        help="The lineage typing called the sample (default: Untyped)."
    )
    args = parser.parse_args()

    input_path = Path(args.master_table_file)
    base_name = input_path.name.replace("-AMAPmastertable.tsv", "")

    # Read and preprocess the data
    df = pd.read_csv(input_path, sep='\t', low_memory=False)
    df.replace("-", np.nan, inplace=True)

    # A gene is a "sample" gene if the 'rvs' column IS empty (NaN).
    # A gene is a "reference" gene if the 'rvs' column is NOT empty.
    is_sample_gene = df['rvs'].isna()

    # Every gene in a master table is the same sample's, and the pipeline says
    # which sample that is and its lineage. They used to be read out of the
    # gene IDs, <name>_<gene>, but the name is <lineage>_<sampleId> after
    # lineage typing and the sample ID alone without it, and a sample ID can
    # hold underscores of its own: L4_H37Rv_test counted as genome H37Rv,
    # samples Beijing_1 and Beijing_2 as one, and H37Rv_test under
    # --skip_typing as lineage H37Rv.
    df['lineage'] = args.lineage or 'Untyped'
    df['genome'] = args.sample or base_name

    # Process COG categories (always available from eggNOG)
    cog_summary = process_cog_categories(df)

    # Process Mycobrowser categories (only if column exists)
    has_mycobrowser = 'Functional_Category' in df.columns and df['Functional_Category'].notna().any()
    if has_mycobrowser:
        mycobrowser_summary = process_mycobrowser_categories(df)
    else:
        mycobrowser_summary = pd.DataFrame(columns=['lineage','genome','Functional_Category','count','percentage'])

    # --- START: SAVE ALL OUTPUT FILES ---
    print("Saving output files...")

    # Define output paths with a consistent naming convention
    cog_output_path = input_path.parent / f"{base_name}-AMAP.COG.tsv"
    myco_output_path = input_path.parent / f"{base_name}-AMAP.MycoBrowser.tsv"
    plot_output_path = input_path.parent / f"{base_name}-AMAP.annot.png"
    no_ref_output_path = input_path.parent / f"{base_name}-AMAP.noRef.tsv"

    # Save summary tables
    cog_summary.to_csv(cog_output_path, sep='\t', index=False)
    if has_mycobrowser:
        mycobrowser_summary.to_csv(myco_output_path, sep='\t', index=False)
        print(f"Successfully saved summary files: {cog_output_path.name}, {myco_output_path.name}")
    else:
        # Write empty file so downstream doesn't crash
        pd.DataFrame().to_csv(myco_output_path, sep='\t', index=False)
        print(f"Mycobrowser data not available — wrote empty {myco_output_path.name}")

    # Save non-reference gene list (where 'rvs' was empty)
    no_ref_df = df.loc[is_sample_gene, ['genome', 'lineage']]
    no_ref_df.to_csv(no_ref_output_path, sep='\t', index=False, header=True)
    print(f"Successfully generated non-reference file: {no_ref_output_path.name}")

    # Generate and save plots
    if has_mycobrowser:
        create_summary_plots(cog_summary, mycobrowser_summary, plot_output_path)
    else:
        # COG-only plot
        fig, ax1 = plt.subplots(1, 1, figsize=(12, 7))
        cog_plot = cog_summary.groupby('COG_category')['count'].sum().reset_index()
        ax1.bar(cog_plot['COG_category'], cog_plot['count'],
                color=plt.cm.viridis(np.linspace(0, 1, len(cog_plot))))
        ax1.set_title('eggNOG Functional Category Distribution', fontsize=14)
        ax1.set_ylabel('Total Gene Count', fontsize=12)
        ax1.tick_params(axis='x', rotation=45)
        ax1.grid(axis='y', linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.savefig(plot_output_path, dpi=300, bbox_inches='tight')
        print(f"Generated COG-only plot at {plot_output_path}")
    # --- END: SAVE ALL OUTPUT FILES ---

if __name__ == "__main__":
    main()
