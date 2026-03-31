#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ast
import pandas as pd
from Bio import SeqIO
import argparse

__author__ = 'Paula Ruiz-Rodriguez'
__credits__ = ['Paula Ruiz-Rodriguez', 'PathoGenOmics']
__license__ = ['GPL 3']
__version__ = '0.1.0'
__maintainer__ = '@paururo'
__email__ = 'paula.ruiz.rodriguez@csic.es'
__status__ = 'finished'
__lab__ = 'PathoGenOmics (I2SysBio)'


def extract_sequences(tsv_filepath, fasta_filepath, output_filepath):
    """
    Extracts sequences from a FASTA file based on a list of IDs from a TSV file.

    Args:
        tsv_filepath (str): Path to the input TSV file.
        fasta_filepath (str): Path to the input multifasta file.
        output_filepath (str): Path for the output FASTA file.
    """
    # Read the TSV file into a pandas DataFrame.
    df = pd.read_csv(tsv_filepath, sep='\t')

    # Extract the sequence IDs from the 'Members' column.
    # The IDs are stored in a string that looks like a list, e.g., "['id1', 'id2']".
    # We clean this string and split it to get individual IDs.
    sequence_ids = []
    for row in df['Members']:
        # Safely parse the Python list literal stored as a string
        try:
            parsed = ast.literal_eval(row)
            if isinstance(parsed, list):
                sequence_ids.extend(parsed)
            else:
                sequence_ids.append(str(parsed))
        except (ValueError, SyntaxError):
            # Fallback: strip and split manually
            cleaned_ids = row.strip("[]'").split("', '")
            sequence_ids.extend(cleaned_ids)

    # Read the FASTA file and create a dictionary of sequences, with IDs as keys.
    # This allows for quick lookups.
    available_sequences = {record.id: record for record in SeqIO.parse(fasta_filepath, 'fasta')}
    
    # Create a list of sequences that match the IDs from the TSV file.
    found_sequences = [available_sequences[seq_id] for seq_id in sequence_ids if seq_id in available_sequences]

    # Write the found sequences to a new FASTA file.
    SeqIO.write(found_sequences, output_filepath, 'fasta')
    print(f"Successfully extracted {len(found_sequences)} sequences to {output_filepath}")


def main():
    """
    Main function to parse command-line arguments and run the script.
    """
    parser = argparse.ArgumentParser(
        description="Extracts specific sequences from a multifasta file based on a list of IDs provided in a TSV file."
    )
    parser.add_argument("tsv", help="Path to the input TSV file containing the sequence IDs in a 'Members' column.")
    parser.add_argument("cds", help="Path to the input multifasta file (e.g., CDS sequences).")
    
    args = parser.parse_args()
    
    # Define the output file name by appending '.dup' to the input CDS file name.
    output_file = args.cds + ".dup"
    
    # Call the function to perform the extraction.
    extract_sequences(args.tsv, args.cds, output_file)

if __name__ == '__main__':
    main()
