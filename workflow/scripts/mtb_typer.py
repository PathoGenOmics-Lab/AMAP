#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import argparse
from Bio import SeqIO
import os

def arguments():
    """
    Define and parse command line arguments.

    Returns:
    - argparse.Namespace: Parsed arguments.
    """

    parser = argparse.ArgumentParser(
        prog='mtbtyper',
        description='A package to use k-mers to predict the _Mycobacterium tuberculosis_ complex species'
    )
    
    parser.add_argument('-i', '--input', help='Input file in fasta format', required=True)
    parser.add_argument('-t', '--snp', help='Marker file', required=True)
    parser.add_argument('-r', '--ref', help='Fasta reference', required=True)
    
    return parser.parse_args()

def generate_kmers(sequence, k):
    """
    Generate all k-mers from a given sequence and their starting positions.

    Parameters:
    - sequence (str): The sequence from which k-mers are to be generated.
    - k (int): The length of k-mers to be generated.

    Returns:
    - dict: A dictionary with k-mers as keys and their starting positions in the sequence as values.

    Example:
    >>> generate_kmers("ACGTA", 3)
    {'ACG': 0, 'CGT': 1, 'GTA': 2}
    """
    kmers = {}
    for i in range(len(sequence) - k + 1):
        kmers[sequence[i:i+k]] = i
    return kmers

def find_markers(genome_sequence, markers_kmers, k):
    """
    Find positions of marker k-mers within a genome sequence.

    Parameters:
    - genome_sequence (str): The genome sequence to search within.
    - markers_kmers (dict): Dictionary of marker k-mers with reference position and lineage.
    - k (int): Length of the k-mers.

    Returns:
    - dict: Marker k-mers found in the genome and their respective positions and lineage.
    """
    
    genome_kmers = generate_kmers(genome_sequence, k)
    matched_markers = {}
    
    for marker_kmer, (ref_position, lineage) in markers_kmers.items():
        if marker_kmer in genome_kmers:
            matched_markers[marker_kmer] = genome_kmers[marker_kmer], ref_position, lineage
            
    return matched_markers

def get_positions(in_filepos):
    """
    Extract reference positions, alternative bases, and lineage information from a TSV file.

    Parameters:
    - in_filepos (str): Path to the input TSV file. The file should have columns: position, alt_base, and lineage. 
                        The first row is expected to be a header and is skipped.

    Returns:
    - tuple: Two dictionaries:
        1. reference_positions: Mapping from positions to alternative bases.
        2. markers_lineage: Mapping from positions to lineage information.

    Example file format:
    position    alt_base    lineage
    123         A           lineageA
    456         C           lineageB
    ...         ...         ...
    """
    
    reference_positions = {}
    markers_lineage = {}
    
    try:
        with open(in_filepos, 'r') as tsvfile:
            tsvreader = csv.DictReader(tsvfile, delimiter="\t")
            for row in tsvreader:
                pos = int(row['position'])
                reference_positions[pos] = row['alt_base']
                markers_lineage[pos] = row['lineage']
    except (IOError, KeyError) as e:
        print(f"Error reading the file or processing its contents: {e}")
        return {}, {}
            
    return reference_positions, markers_lineage

def get_ref(fasta_file):
    """
    Extracts the sequence from a given FASTA file and returns it in uppercase.

    Parameters:
    - fasta_file (str): Path to the input FASTA file.

    Returns:
    - str: The sequence from the FASTA file, converted to uppercase.

    Raises:
    - FileNotFoundError: If the provided FASTA file cannot be found.
    - ValueError: If the file cannot be read as a FASTA file or contains multiple records.

    Example:
    >>> seq = get_ref("path_to_file.fasta")
    >>> print(seq)
    "ATCGAAGT..."

    Note:
    This function expects the FASTA file to contain a single sequence record.
    """
    
    try:
        ref_seq_record = SeqIO.read(fasta_file, "fasta")
        ref_seq = str(ref_seq_record.seq).upper()
        return ref_seq
    except FileNotFoundError:
        raise FileNotFoundError(f"The file {fasta_file} was not found.")
    except ValueError:
        raise ValueError(f"The file {fasta_file} either contains multiple sequence records or cannot be read as a FASTA file.")

def generate_markerkmer(reference_positions, ref_seq, markers_lineage, k):
    """
    Generate k-mers centered on specified positions with alternative bases.

    Parameters:
    - reference_positions (dict): Positions mapped to alternative bases.
    - ref_seq (str): The reference sequence.
    - markers_lineage (dict): Positions mapped to lineage information.
    - k (int): Length of k-mers (should be odd).

    Returns:
    - dict: k-mers mapped to their position and lineage.

    Note:
    Positions near the sequence edges are ignored.
    
    Example:
    >>> reference_positions = {2: 'A', 25: 'G'}
    >>> ref_seq = "ATCGTACGTAGCTAGCTAGCTGACTG"
    >>> markers_lineage = {2: 'lineageA', 25: 'lineageB'}
    >>> k = 5
    >>> generate_markerkmer(reference_positions, ref_seq, markers_lineage, k)
    {'ATAGT': (2, 'lineageA'), 'TAGGT': (25, 'lineageB')}
    """
    
    markers_kmers = {}
    for pos, alt_base in reference_positions.items():
        if pos > k//2 and pos < len(ref_seq) - k//2:
            kmer = ref_seq[pos - k//2 - 1: pos + k//2].upper()
            kmer = kmer[:k//2] + alt_base + kmer[k//2 + 1:]
            markers_kmers[kmer] = pos, markers_lineage[pos]
            
    return markers_kmers

def analyze_genome(args):
    """
    Analyzes a genome sequence and matches it against known marker k-mers.

    Parameters:
    - args (tuple):
        - genome_name (str): Name of the genome being analyzed.
        - fasta_path (str): Path to the FASTA file containing the genome sequence.
        - markers_kmers (dict): Dictionary containing marker k-mers and their associated information.
        - k (int): Length of the k-mers.

    Returns:
    - tuple: 
        1. genome_name (str): The name of the genome.
        2. identified_lineages (set): A set of identified unique lineages. 
           If the FASTA file does not exist, this set will be empty.

    Note:
    - If the FASTA file is not found or no markers are matched, the result will reflect this accordingly.
    """
    genome_name, fasta_path, markers_kmers, k = args

    genome_seq = get_ref(fasta_path)
    matched_markers = find_markers(genome_seq, markers_kmers, k)

    identified_lineages = set()
    for _, (_, _, lineage) in matched_markers.items():
        identified_lineages.add(lineage)

    return genome_name, identified_lineages

def rename_fasta_based_on_lineages(fasta_path, genome_name, identified_lineages):
    """
    Renames the given FASTA file and its sequence ID based on the identified lineages.

    Parameters:
    - fasta_path (str): The path to the input FASTA file.
    - genome_name (str): The name of the genome.
    - identified_lineages (set): A set of identified lineages for the genome.
    """

    # Convert identified lineages set to a string separated by pipes
    lineages_str = ";".join(sorted(identified_lineages))

    # Construct the new ID and filename
    new_id = f"{lineages_str}+_+{genome_name}"
    new_file_path = os.path.join(os.path.dirname(fasta_path), new_id + ".fasta")

    # Read the original FASTA, change the ID, and write to the new file
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if records:
        records[0].id = new_id
        records[0].description = ""  # Clear the description to ensure only the ID is written
        SeqIO.write(records, new_file_path, "fasta")