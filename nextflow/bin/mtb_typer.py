#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
from Bio import SeqIO
import os

def generate_kmers(sequence, k):
    '''
    Generate all k-mers from a given sequence and their starting positions.

    Parameters:
    - sequence (str): The sequence from which k-mers are to be generated.
    - k (int): The length of k-mers to be generated.

    Returns:
    - dict: A dictionary with k-mers as keys and their starting positions in the sequence as values.

    Example:
    >>> generate_kmers('ACGTA', 3)
    {'ACG': 0, 'CGT': 1, 'GTA': 2}
    '''
    kmers = {}
    for i in range(len(sequence) - k + 1):
        kmers[sequence[i:i+k]] = i
    return kmers

def find_markers(genome_sequence, markers_kmers, k):
    '''
    Find positions of marker k-mers within a genome sequence.

    Parameters:
    - genome_sequence (str): The genome sequence to search within.
    - markers_kmers (dict): Dictionary of marker k-mers with reference position and lineage.
    - k (int): Length of the k-mers.

    Returns:
    - dict: Marker k-mers found in the genome and their respective positions and lineage.
    '''
    
    genome_kmers = generate_kmers(genome_sequence, k)
    matched_markers = {}
    
    for marker_kmer, (ref_position, lineage) in markers_kmers.items():
        if marker_kmer in genome_kmers:
            matched_markers[marker_kmer] = genome_kmers[marker_kmer], ref_position, lineage
            
    return matched_markers

def get_positions(in_filepos):
    '''
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
    '''
    
    reference_positions = {}
    markers_lineage = {}
    
    try:
        with open(in_filepos, 'r') as tsvfile:
            tsvreader = csv.DictReader(tsvfile, delimiter='\t')
            for row in tsvreader:
                pos = int(row['position'])
                reference_positions[pos] = row['alt_base']
                markers_lineage[pos] = row['lineage']
    except (IOError, KeyError) as e:
        print(f'Error reading the file or processing its contents: {e}')
        return {}, {}
            
    return reference_positions, markers_lineage