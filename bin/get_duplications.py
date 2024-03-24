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

def group_related_pairs(df):
    G = nx.Graph()

    for idx, row in df.iterrows():
        G.add_edge(row['qseqid'], row['sseqid'])

    components = list(nx.connected_components(G))
    df_groups = pd.DataFrame(columns=['Group', 'Members'])

    for i, component in enumerate(components):
        df_groups = df_groups.append({'Group': i, 'Members': list(component)}, ignore_index=True)
    df_groups['Length'] = df_groups['Members'].apply(lambda x: len(x))
    return df_groups