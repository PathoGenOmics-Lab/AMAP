#! /usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network

__author__ = 'Paula Ruiz-Rodriguez'
__credits__ = ['Paula Ruiz-Rodriguez', 'PathoGenOmics']
__license__ = 'GPL 3'
__version__ = '0.1.0'
__maintainer__ = '@paururo'
__email__ = 'paula.ruiz-rodriguez@uv.es'
__status__ = 'developing'
__lab__ = 'PathoGenOmics (I2SysBio)'

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
    group_data = []

    for i, component in enumerate(components):
        group_data.append({'Group': i, 'Members': list(component)})

    df_groups = pd.DataFrame(group_data)
    df_groups['Length'] = df_groups['Members'].apply(len)
    return df_groups

def plot_network(df, output_filename):
    G = nx.Graph()

    for idx, row in df.iterrows():
        G.add_edge(row['qseqid'], row['sseqid'])

    # Calculate connected components (clusters)
    clusters = list(nx.connected_components(G))
    
    # Calculate the number of rows and columns for the grid of subplots
    num_clusters = len(clusters)
    num_rows = int(np.ceil(np.sqrt(num_clusters)))
    num_cols = num_rows if num_clusters > num_rows * (num_rows - 1) else num_rows - 1
    
    # Create subplots in a grid
    fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(20, 20))
    axes = axes.flatten()  # Flatten the axes array for easy indexing
    
    # Adjust the spacing between the subplots
    plt.subplots_adjust(wspace=0.3, hspace=0.3)  # Modify these values as needed for your data
    
    for ax in axes[num_clusters:]:  # Hide unused subplots
        ax.set_visible(False)
    
    for i, nodes in enumerate(clusters):
        subgraph = G.subgraph(nodes)
        pos = nx.spring_layout(subgraph, k=0.15, iterations=50)
        ax = axes[i]
        nx.draw_networkx(subgraph, pos, ax=ax, with_labels=True, node_size=100, node_color='lightblue', edge_color='gray', font_size=8)
        ax.set_title(f'Cluster {i+1}', size=15)
        ax.axis('off')
    
    plt.tight_layout(pad=2.0)  # Add padding around the figure
    plt.savefig(output_filename, format='png', dpi=300)  # Adjust DPI as needed


def plot_network2(df, output_filename):
    G = nx.Graph()

    # Add edges from your dataframe to the graph
    for idx, row in df.iterrows():
        G.add_edge(row['qseqid'], row['sseqid'])

    # Create a pyvis network object from the NetworkX graph
    net = Network(height='1080px', width='100%', bgcolor='#222222', font_color='white')
    net.from_nx(G)

    # Set options (you can adjust these for better visualization)
    net.set_options("""
    var options = {
      "nodes": {
        "font": {
          "size": 12,
          "color": "white"
        }
      },
      "edges": {
        "color": {
          "inherit": true
        },
        "smooth": false
      },
      "physics": {
        "barnesHut": {
          "springLength": 200
        }
      }
    }
    """)

    net.show(output_filename, notebook=False)


def main():
    parser = argparse.ArgumentParser(description="Read blast and generate cluster of duplications")
    parser.add_argument("file", help="Blastp output")
    args = parser.parse_args()

    df = read_and_create(args.file)
    df_with_groups = group_related_pairs(df)

    # Save the groups to a file
    df_with_groups.to_csv(args.file + "_groups.tsv", sep='\t', index=False)

    # Plot and save the network
    plot_network(df, args.file + "_network.png")
    plot_network2(df, args.file + "_network_interactive.html")

if __name__ == "__main__":
    main()
    #python3 generate_plot.py resultados_blastp_gaps.txt
    #blastp -query nt_correctedH37Rv.translated.fas -db base_de_h37rv -out resultados_blastp_gaps.txt -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen nident"
    #makeblastdb -in nt_correctedH37Rv.fa -dbtype prot -out base_de_h37rv
