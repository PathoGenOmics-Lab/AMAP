#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script reads a master annotation table, calculates the frequency of
eggNOG COG categories and Mycobrowser Functional Categories, saves these
counts to summary TSV files, and generates a bar plot visualizing the
category distributions.
"""

import argparse
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
__version__ = '1.0.0'
__maintainer__ = '@paururo'
__email__ = 'paula.ruiz.rodriguez@csic.es'
__status__ = 'finished'
__lab__ = 'PathoGenOmics (I2SysBio)'


def process_cog_categories(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates the count and percentage of each COG category.

    The 'COG_category' column can contain multiple categories (e.g., "JKS").
    This function expands them so that each category is counted individually.
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
    
    # Handle missing values, convert to uppercase, and then to lists of characters
    cog_df['COG_category'] = cog_df['COG_category'].fillna('Unassigned').str.upper().apply(list)
    
    # Create a new row for each COG category character
    cog_df_exploded = cog_df.explode('COG_category')

    # Group by lineage, genome, and category to count occurrences
    cog_counts = cog_df_exploded.groupby(
        ['lineage', 'genome', 'COG_category']
    ).size().reset_index(name='count')

    # Calculate percentage based on the total number of genes
    cog_counts['percentage'] = (cog_counts['count'] / total_genes) * 100
    
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
    args = parser.parse_args()
    
    input_path = Path(args.master_table_file)
    base_name = input_path.name.replace("-AMAPmastertable.tsv", "")

    # Read and preprocess the data
    df = pd.read_csv(input_path, sep='\t', low_memory=False)
    df.replace("-", np.nan, inplace=True)
    
    # --- START: ID PROCESSING LOGIC ---
    # Use the first column as the source for IDs
    id_source = df.iloc[:, 0].astype(str)

    # A gene is a "sample" gene if the 'rvs' column IS empty (NaN).
    # A gene is a "reference" gene if the 'rvs' column is NOT empty.
    is_sample_gene = df['rvs'].isna()
    is_ref_gene = ~is_sample_gene

    # Determine the primary lineage and genome ID from the sample genes (where 'rvs' is empty)
    sample_ids = id_source[is_sample_gene]
    if not sample_ids.empty:
        # Find the first valid sample ID to extract lineage and genome
        first_sample_id_parts = sample_ids.iloc[0].split('_')
        if len(first_sample_id_parts) >= 2:
            sample_lineage = first_sample_id_parts[0]
            sample_genome_id = first_sample_id_parts[1]
        else:
            sample_lineage = 'Unknown_Lineage'
            sample_genome_id = 'Unknown_Genome'
    else:
        # Fallback if no sample genes are found
        sample_lineage = 'Unknown_Lineage'
        sample_genome_id = 'Unknown_Genome'

    # Assign 'lineage' column
    # For reference genes, use the sample's primary lineage. For others, extract it from their ID.
    df['lineage'] = np.where(
        is_ref_gene,
        sample_lineage,
        id_source.str.split('_').str[0]
    )

    # Assign 'genome' column
    # For reference genes, prepend the sample's genome ID.
    # For sample genes, use everything after the first underscore.
    df['genome'] = np.where(
        is_ref_gene,
        sample_genome_id + '_' + id_source.str.split('_').str[0], # e.g., G26689_Rv0029
        id_source.str.split('_', n=1).str[1] # e.g., 'G26689_2049' from 'L8_G26689_2049'
    )
    # --- END: ID PROCESSING LOGIC ---

    # Process COG and Mycobrowser categories
    cog_summary = process_cog_categories(df)
    mycobrowser_summary = process_mycobrowser_categories(df)

    # --- START: SAVE ALL OUTPUT FILES ---
    print("Saving output files...")

    # Define output paths with a consistent naming convention
    cog_output_path = input_path.parent / f"{base_name}-AMAP.COG.tsv"
    myco_output_path = input_path.parent / f"{base_name}-AMAP.MycoBrowser.tsv"
    plot_output_path = input_path.parent / f"{base_name}-AMAP.annot.png"
    no_ref_output_path = input_path.parent / f"{base_name}-AMAP.noRef.tsv"

    # Save summary tables
    cog_summary.to_csv(cog_output_path, sep='\t', index=False)
    mycobrowser_summary.to_csv(myco_output_path, sep='\t', index=False)
    print(f"Successfully saved summary files: {cog_output_path.name}, {myco_output_path.name}")

    # Save non-reference gene list (where 'rvs' was empty)
    no_ref_df = df.loc[is_sample_gene, ['genome', 'lineage']]
    no_ref_df.to_csv(no_ref_output_path, sep='\t', index=False, header=True)
    print(f"Successfully generated non-reference file: {no_ref_output_path.name}")

    # Generate and save plots
    create_summary_plots(cog_summary, mycobrowser_summary, plot_output_path)
    # --- END: SAVE ALL OUTPUT FILES ---

if __name__ == "__main__":
    main()
