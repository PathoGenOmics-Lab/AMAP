#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script analyzes BLASTp output to identify, characterize, and visualize
networks of homologous gene communities within and between genomes. It performs
the following steps:
1.  Filters BLAST results to retain significant hits.
2.  Builds a protein-protein interaction network for each genome.
3.  Detects gene communities (modules) within each network using the Louvain method.
4.  Identifies homologous relationships between communities across different genomes.
5.  Calculates various network topology metrics for each community.
6.  Generates a detailed TSV report of the communities and their metrics.
7.  Creates an interactive HTML visualization of the homologous community network.
"""

import argparse
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pyvis.network import Network
from pathlib import Path
from typing import List, Dict, Set, Tuple, Any

__author__ = 'Paula Ruiz-Rodriguez'
__credits__ = ['Paula Ruiz-Rodriguez', 'PathoGenOmics']
__license__ = ['GPL 3']
__version__ = '0.1.0'
__maintainer__ = '@paururo'
__email__ = 'paula.ruiz.rodriguez@csic.es'
__status__ = 'finished'
__lab__ = 'PathoGenOmics (I2SysBio)'


# --- 1. Data Loading and Preparation ---

def load_and_filter_blast_results(blast_filepath: Path) -> pd.DataFrame:
    """
    Reads a BLASTp output file (outfmt 6) and filters for significant alignments.

    Filtering criteria:
    - Removes self-hits.
    - Keeps hits where the alignment covers at least 80% of the query or the subject.

    Args:
        blast_filepath: Path to the BLASTp output file.

    Returns:
        A pandas DataFrame with the filtered BLAST results.
    """
    column_names = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gaps", "qstart",
        "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "nident"
    ]
    try:
        df = pd.read_csv(blast_filepath, sep='\t', comment='#', header=None, names=column_names)
    except FileNotFoundError:
        print(f"Error: Input file not found at {blast_filepath}")
        exit(1)

    # Calculate query and subject coverage percentages
    df['query_coverage'] = (df['nident'] / df['qlen']) * 100
    df['subject_coverage'] = (df['nident'] / df['slen']) * 100

    # Filter out self-hits and low-coverage hits
    df_filtered = df[df['qseqid'] != df['sseqid']]
    df_final = df_filtered[
        (df_filtered['query_coverage'] >= 80) | (df_filtered['subject_coverage'] >= 80)
    ].copy()

    return df_final


def get_genome_ids(df: pd.DataFrame) -> Set[str]:
    """
    Extracts unique genome identifiers from the sequence IDs in the BLAST results.
    Assumes a sequence ID format of 'GenomeID_GeneID...'.

    Args:
        df: The filtered BLAST results DataFrame.

    Returns:
        A set of unique genome identifiers.
    """
    all_seq_ids = set(df['qseqid'].unique()) | set(df['sseqid'].unique())
    genome_ids = {"_".join(seq_id.split('_')[:2]) for seq_id in all_seq_ids}
    return genome_ids


# --- 2. Network Construction and Community Detection ---

def build_intra_genome_networks(df: pd.DataFrame, genome_ids: Set[str]) -> Dict[str, nx.Graph]:
    """
    Builds a protein-protein interaction network for each individual genome.

    Args:
        df: The filtered BLAST results DataFrame.
        genome_ids: A set of unique genome identifiers.

    Returns:
        A dictionary mapping each genome ID to its corresponding NetworkX graph.
    """
    networks = {}
    for genome_id in genome_ids:
        genome_df = df[
            df['qseqid'].str.startswith(genome_id) & df['sseqid'].str.startswith(genome_id)
        ]
        G = nx.from_pandas_edgelist(genome_df, 'qseqid', 'sseqid')
        networks[genome_id] = G
    return networks


def detect_and_label_communities(networks: Dict[str, nx.Graph]) -> Dict[str, Dict[str, Set[str]]]:
    """
    Detects communities in each network using the Louvain method and labels them.

    Args:
        networks: A dictionary of genome networks.

    Returns:
        A nested dictionary: {genome_id: {community_name: {node1, node2, ...}}}.
    """
    labeled_communities = {}
    for genome_id, G in networks.items():
        # Use networkx's implementation of the Louvain method
        communities = nx.community.louvain_communities(G)
        labeled_communities[genome_id] = {
            f"Community_{i}": comm for i, comm in enumerate(communities)
        }
    return labeled_communities


# --- 3. Homologous Community Analysis ---

def find_homologous_communities(
    labeled_communities: Dict[str, Dict[str, Set[str]]],
    blast_df: pd.DataFrame
) -> Dict[Tuple[str, str], List[Tuple[str, str]]]:
    """
    Identifies homologous relationships between communities across different genomes.

    Two communities are considered homologous if at least one protein from the first
    community has a significant BLAST hit to a protein in the second community.

    Args:
        labeled_communities: The dictionary of labeled communities for each genome.
        blast_df: The original filtered BLAST results DataFrame.

    Returns:
        A dictionary mapping a community (genome_id, community_name) to a list of
        its homologous communities [(homolog_genome, homolog_community), ...].
    """
    # Create a set of all homologous pairs for fast lookups
    homologous_pairs = set(zip(blast_df['qseqid'], blast_df['sseqid']))

    homologous_map = {}
    genome_ids = list(labeled_communities.keys())

    # Compare every community from every genome against all others
    for i, g1_id in enumerate(genome_ids):
        for j in range(i, len(genome_ids)):
            g2_id = genome_ids[j]

            for c1_name, c1_nodes in labeled_communities[g1_id].items():
                for c2_name, c2_nodes in labeled_communities[g2_id].items():
                    # Skip self-comparison if it's the same community in the same genome
                    if g1_id == g2_id and c1_name == c2_name:
                        continue

                    # Check for any homologous link between the two communities
                    is_homologous = any(
                        (n1, n2) in homologous_pairs or (n2, n1) in homologous_pairs
                        for n1 in c1_nodes for n2 in c2_nodes
                    )

                    if is_homologous:
                        homologous_map.setdefault((g1_id, c1_name), []).append((g2_id, c2_name))
                        # Add the reverse relationship if genomes are different
                        if g1_id != g2_id:
                            homologous_map.setdefault((g2_id, c2_name), []).append((g1_id, c1_name))
    return homologous_map


# --- 4. Network Metrics Calculation ---

def calculate_subgraph_metrics(subgraph: nx.Graph) -> Dict[str, Any]:
    """Calculates various topological metrics for a given subgraph."""
    metrics = {}
    if not subgraph or subgraph.number_of_nodes() == 0:
        return {
            'avg_degree': 0, 'density': 0, 'avg_clustering': 0,
            'avg_path_length': np.nan, 'modularity': np.nan
        }

    metrics['avg_degree'] = np.mean(list(dict(subgraph.degree()).values()))
    metrics['density'] = nx.density(subgraph)
    metrics['avg_clustering'] = nx.average_clustering(subgraph)
    metrics['modularity'] = nx.community.modularity(subgraph, [set(subgraph.nodes())])
    
    if nx.is_connected(subgraph):
        metrics['avg_path_length'] = nx.average_shortest_path_length(subgraph)
    else:
        # Calculate for the largest connected component
        largest_cc = max(nx.connected_components(subgraph), key=len)
        largest_subgraph = subgraph.subgraph(largest_cc)
        if largest_subgraph.number_of_nodes() > 1:
            metrics['avg_path_length'] = nx.average_shortest_path_length(largest_subgraph)
        else:
            metrics['avg_path_length'] = 0

    return metrics


def get_community_report_row(
    genome_id: str,
    community_name: str,
    labeled_communities: Dict,
    network: nx.Graph
) -> List[Any]:
    """Generates a single row of data for the final report DataFrame."""
    nodes = labeled_communities[genome_id][community_name]
    subgraph = network.subgraph(nodes)
    metrics = calculate_subgraph_metrics(subgraph)
    
    # Extract gene names (assuming format '...;GeneName;...')
    elements = [node.split(';')[1] if len(node.split(';')) > 1 else 'Unknown' for node in nodes]

    return [
        genome_id,
        community_name,
        len(nodes),
        dict(pd.Series(elements).value_counts()),
        metrics.get('avg_degree', 0),
        metrics.get('density', 0),
        metrics.get('avg_clustering', 0),
        metrics.get('avg_path_length', np.nan),
        metrics.get('modularity', np.nan),
        ", ".join(sorted(list(nodes)))
    ]


def create_summary_report(
    labeled_communities: Dict,
    homologous_map: Dict,
    networks: Dict
) -> pd.DataFrame:
    """Creates a DataFrame summarizing all communities and their homologous relationships."""
    report_rows = []
    processed_pairs = set()

    for (g1_id, c1_name), homologs in homologous_map.items():
        for (g2_id, c2_name) in homologs:
            # Create a sorted tuple to represent the pair uniquely
            pair = tuple(sorted(((g1_id, c1_name), (g2_id, c2_name))))
            if pair in processed_pairs:
                continue
            
            row1 = get_community_report_row(g1_id, c1_name, labeled_communities, networks[g1_id])
            row2 = get_community_report_row(g2_id, c2_name, labeled_communities, networks[g2_id])
            report_rows.append(row1 + row2)
            processed_pairs.add(pair)

    columns = [
        "Genome1", "Community1", "Size1", "Elements1", "AvgDegree1", "Density1",
        "AvgClustering1", "AvgPathLength1", "Modularity1", "Nodes1",
        "Genome2", "Community2", "Size2", "Elements2", "AvgDegree2", "Density2",
        "AvgClustering2", "AvgPathLength2", "Modularity2", "Nodes2"
    ]
    return pd.DataFrame(report_rows, columns=columns)


# --- 5. Visualization ---

def visualize_network(
    labeled_communities: Dict,
    homologous_map: Dict,
    output_path: Path
):
    """Creates an interactive HTML network visualization using pyvis."""
    net = Network(height='100vh', width='100%', bgcolor='#FFFFFF', font_color='black', notebook=False)
    net.force_atlas_2based(gravity=-50, central_gravity=0.01, spring_length=100)

    # Assign a unique color to each genome
    genome_ids = list(labeled_communities.keys())
    colors = plt.colormaps['tab20'].resampled(len(genome_ids))
    genome_colors = {gid: mcolors.to_hex(colors(i)) for i, gid in enumerate(genome_ids)}

    # Add nodes (communities)
    for genome_id, communities in labeled_communities.items():
        for comm_name, nodes in communities.items():
            node_id = f"{genome_id}:{comm_name}"
            node_size = int(10 + 40 * (len(nodes) / 100)) # Scale size
            net.add_node(
                node_id,
                label=comm_name,
                size=max(10, node_size),
                color=genome_colors[genome_id],
                title=f"Genome: {genome_id}\nCommunity: {comm_name}\nSize: {len(nodes)}"
            )

    # Add edges (homologous relationships)
    for (g1_id, c1_name), homologs in homologous_map.items():
        for (g2_id, c2_name) in homologs:
            # Avoid adding self-loops for intra-genome homology
            if g1_id == g2_id and c1_name == c2_name:
                continue
            
            source_node = f"{g1_id}:{c1_name}"
            target_node = f"{g2_id}:{c2_name}"
            net.add_edge(source_node, target_node)

    net.show(str(output_path), notebook=False)
    print(f"Successfully generated interactive network at {output_path}")


# --- Main Execution ---

def main():
    """Main function to parse arguments and run the network analysis workflow."""
    parser = argparse.ArgumentParser(
        description="Analyze BLAST results to find and visualize homologous gene communities."
    )
    parser.add_argument(
        "blast_file",
        help="Path to the input BLASTp output file (format 6)."
    )
    args = parser.parse_args()
    
    input_path = Path(args.blast_file)
    output_tsv = input_path.parent / "homologous_communities.tsv"
    output_html = input_path.parent / "homologous_communities.html"

    # --- Workflow ---
    # 1. Load and filter data
    blast_df = load_and_filter_blast_results(input_path)
    genome_ids = get_genome_ids(blast_df)

    # 2. Build networks and detect communities
    networks = build_intra_genome_networks(blast_df, genome_ids)
    labeled_communities = detect_and_label_communities(networks)

    # 3. Find homologous relationships
    homologous_map = find_homologous_communities(labeled_communities, blast_df)

    # 4. Generate summary report and save
    report_df = create_summary_report(labeled_communities, homologous_map, networks)
    report_df.to_csv(output_tsv, sep="\t", index=False)
    print(f"Successfully generated report at {output_tsv}")

    # 5. Generate interactive visualization
    visualize_network(labeled_communities, homologous_map, output_html)

if __name__ == "__main__":
    main()
