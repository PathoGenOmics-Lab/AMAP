#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import pandas as pd
import os
from io import StringIO
from typing import List, Any, Tuple

def get_seq(chck_rds: list, fasta_file: str, fasta_name: str) -> List[Tuple[str, int]]:
    chck_rds_set = set(chck_rds)  # Convert the list to a set for faster searching
    selected_seqs = [seq_record for seq_record in SeqIO.parse(fasta_file, "fasta") if seq_record.id in chck_rds_set] # Use list comprehension for faster, more readable code
    SeqIO.write(selected_seqs, fasta_name, "fasta") # Write the selected sequences to a fasta file
    seq_lengths = [(seq_record.id, len(seq_record.seq)) for seq_record in selected_seqs] # Use list comprehension to create the list of sequence lengths
    return seq_lengths

def create_database(fasta: str, basename: str) -> None:
    os.system(f'makeblastdb -in {fasta} -dbtype nucl -out {basename}')

def get_tsv(tsvfile: str) -> List[List[str]]:
    with open(tsvfile,'r') as intsv:
        return [line.strip('\n').split('\t') for line in intsv]
