#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

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