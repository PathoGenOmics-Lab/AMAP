#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
from Bio import SeqIO
import numpy as np

__author__ = 'Paula Ruiz-Rodriguez'
__credits__ = ['Paula Ruiz-Rodriguez', 'PathoGenOmics']
__license__ = ['GPL 3']
__version__ = '0.1.0'
__maintainer__ = '@paururo'
__email__ = 'paula.ruiz.rodriguez@csic.es'
__status__ = 'finished'
__lab__ = 'PathoGenOmics (I2SysBio)'


def main():
    """
    Renames FASTA sequence headers based on EggNOG COG category annotations.
    """
    parser = argparse.ArgumentParser(description="Rename FASTA headers using EggNOG annotation.")
    parser.add_argument("egg_nog", help="Input EggNOG annotation file (from emapper).")
    parser.add_argument("cds", help="Input CDS multifasta file.")
    args = parser.parse_args()

    try:
        # Loading the EggNOG annotation.
        # Use comment='#' to skip all comment/header lines regardless of version,
        # then manually set column names from the emapper format.
        df = pd.read_csv(args.egg_nog, sep='\t', comment='#', header=None)
        emapper_columns = [
            "id", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
            "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
            "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
            "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"
        ]
        # Handle files with more or fewer columns than expected
        if len(df.columns) >= len(emapper_columns):
            df.columns = emapper_columns + [f"extra_{i}" for i in range(len(df.columns) - len(emapper_columns))]
        else:
            df.columns = emapper_columns[:len(df.columns)]

    except FileNotFoundError:
        print(f"Error: Annotation file not found at path: {args.egg_nog}")
        return
    except Exception as e:
        print(f"Error reading the annotation file: {e}")
        return


    # Replace hyphens ('-') with NaN (Not a Number) for easier handling
    df.replace("-", np.nan, inplace=True)
    
    # Create a name for the output file
    output_file = args.cds + ".renamed"

    print(f"Reading sequences from: {args.cds}")
    print(f"Writing renamed sequences to: {output_file}")

    renamed_count = 0
    total_count = 0

    with open(args.cds, 'r') as fasta_file, open(output_file, 'w') as output:
        for record in SeqIO.parse(fasta_file, "fasta"):
            total_count += 1
            # Extract the sequence ID (e.g., 'L8_G26689_Rv0001' from 'L8_G26689_Rv0001;+;0-1524')
            seq_id = record.id.split(';')[0]

            # Match by stripping semicolons from both sides.
            # eggNOG query preserves the full FAA header (e.g., 'L8_G26689_Rv0001;+;0-1524')
            # while FASTA record.id also has the full header. We need to match the
            # base ID (before first ';') from both.
            if 'simple_id' not in df.columns:
                df['simple_id'] = df['id'].astype(str).str.split(';').str[0]
            cog_row = df[df['simple_id'] == seq_id]

            # If we find a match in the annotation file
            if not cog_row.empty:
                # Get the COG category value
                cog_value = cog_row['COG_category'].iloc[0]
                
                # Check that the value is not null (it could be a '-' we replaced with NaN)
                if pd.notna(cog_value):
                    if "dup" in seq_id:
                        new_id = seq_id + ";" + str(cog_value)
                    else:
                        new_id = "_".join(seq_id.split("_")[0:3]) + ";" + str(cog_value)
                    
                    record.id = new_id
                    record.description = ''  # Clear the description for a clean header
                    renamed_count += 1

            # Write the record (either modified or original) to the output file
            SeqIO.write(record, output, "fasta")

    print(f"\nProcess completed.")
    print(f"Total sequences processed: {total_count}")
    print(f"Successfully renamed sequences: {renamed_count}")


if __name__ == '__main__':
    main()
