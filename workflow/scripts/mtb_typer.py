#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import argparse
from Bio import SeqIO

def arguments():
    """
    Define and parse command line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
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