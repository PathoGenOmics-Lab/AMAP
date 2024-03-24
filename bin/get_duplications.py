#! /usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np
import networkx as nx

def read_and_create(p_file):
    df = pd.read_csv(p_file, sep='\t')
    names_columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gaps", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs", "qcovhsp", "qlen", "slen", "nident"]
    df.columns = names_columns
    df['same_seq'] = df.iloc[:, 0] == df.iloc[:, 1]
    df['p_query'] = (df.iloc[:, 16] / df.iloc[:, 14]).replace([np.inf, -np.inf], np.nan) * 100
    df['p_subject'] = (df.iloc[:, 16] / df.iloc[:, 15]).replace([np.inf, -np.inf], np.nan) * 100
    df_filtered = df.loc[df['same_seq'] == False]
    df_final_filtered = df_filtered.loc[(df_filtered['p_query'] >= 80) | (df_filtered['p_subject'] >= 80)]
    return df_final_filtered