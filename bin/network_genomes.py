#! /usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
import community as community_louvain
from networkx.algorithms import community
import matplotlib.colors

__author__ = 'Paula Ruiz-Rodriguez'
__credits__ = ['Paula Ruiz-Rodriguez', 'PathoGenOmics']
__license__ = 'GPL 3'
__version__ = '0.1.0'
__maintainer__ = '@paururo'
__email__ = 'paula.ruiz-rodriguez@uv.es'
__status__ = 'developing'
__lab__ = 'PathoGenOmics (I2SysBio)'

def read_and_create(p_file):
    df = pd.read_csv(p_file, sep='\t', comment='#')
    names_columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gaps", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "nident"]
    df.columns = names_columns
    df['same_seq'] = df.iloc[:, 0] == df.iloc[:, 1]
    df['p_query'] = (df.iloc[:, 14] / df.iloc[:, 12]).replace([np.inf, -np.inf], np.nan) * 100
    df['p_subject'] = (df.iloc[:, 14] / df.iloc[:, 13]).replace([np.inf, -np.inf], np.nan) * 100
    df_filtered = df.loc[df['same_seq'] == False]
    df_final_filtered = df_filtered.loc[(df_filtered['p_query'] >= 80) | (df_filtered['p_subject'] >= 80)]
    return df_final_filtered

def extract_unique_sequences(df):
    """ Extrae todas las secuencias únicas del dataframe, tomando solo los primeros dos elementos de las secuencias. """
    sequences = set()

    for seq in df['qseqid'].tolist() + df['sseqid'].tolist():
        # Dividir la secuencia y tomar los primeros dos elementos
        split_seq = seq.split('_')[:2]
        # Unir los elementos y añadir al conjunto de secuencias
        sequences.add('_'.join(split_seq))

    return sequences
def create_network(df, sequence):
    """ Crea una red para una secuencia (genoma) específica. """
    filtered_df = df[(df['qseqid'].str.startswith(sequence)) & (df['sseqid'].str.startswith(sequence))]
    G = nx.Graph()
    for idx, row in filtered_df.iterrows():
        G.add_edge(row['qseqid'], row['sseqid'])
    return G

def create_pairs_set(df):
    """ Crea un conjunto con todos los pares de secuencias en ambas direcciones. """
    pairs = set()
    for _, row in df.iterrows():
        pairs.add((row['qseqid'], row['sseqid']))
        pairs.add((row['sseqid'], row['qseqid']))
    return pairs

def find_homologous_communities(labeled_communities, pairs):
    homologous_communities = {}
    for genome1, communities1 in labeled_communities.items():
        for genome2, communities2 in labeled_communities.items():
            if genome1 != genome2:
                for com1_name, com1_members in communities1.items():
                    for com2_name, com2_members in communities2.items():
                        if any((m1, m2) in pairs or (m2, m1) in pairs for m1 in com1_members for m2 in com2_members):
                            homologous_communities.setdefault((genome1, com1_name), []).append((genome2, com2_name))
    return homologous_communities

def detect_communities(network):
    communities = community.louvain_communities(network)
    return communities

def label_communities(communities_by_genome):
    labeled_communities = {}
    for genome, communities in communities_by_genome.items():
        labeled_communities[genome] = {f"Community_{i}": comm for i, comm in enumerate(communities)}
    return labeled_communities



def get_unique_colors(n):
    """Genera n colores únicos en formato RGB."""
    colors = plt.cm.get_cmap('tab20', n)
    return [colors(i) for i in range(n)]

def assign_colors_to_genomes(genomes):
    """Asigna un color único a cada genoma."""
    unique_colors = get_unique_colors(len(genomes))
    return {genome: matplotlib.colors.to_hex(color) for genome, color in zip(genomes, unique_colors)}


def visualize_homologous_communities_with_pyvis(homologous_communities, labeled_communities):
    net = Network(height='1080px', width='100%', bgcolor='#FFFFFF', font_color='black')

    # Preparar los datos: identificar los tamaños de las comunidades y las relaciones homólogas
    community_sizes = {}
    for genome, communities in labeled_communities.items():
        for community_name, nodes in communities.items():
            community_sizes[(genome, community_name)] = len(nodes)

    # Normalizar los tamaños para que sean más manejables en la visualización
    max_size = max(community_sizes.values())
    min_size, max_vis_size = 10, 50  # Tamaños mínimo y máximo para la visualización
    normalized_sizes = {key: int(min_size + (max_vis_size - min_size) * (size / max_size)) for key, size in community_sizes.items()}

    # Asignar colores a los genomas
    genomes = labeled_communities.keys()
    genome_colors = assign_colors_to_genomes(genomes)
    # Añadir nodos representando las comunidades
    for (genome, community), size in normalized_sizes.items():
        rgba_color = genome_colors[genome]  # Obtener color RGBA para el genoma
        hex_color = matplotlib.colors.to_hex(rgba_color)  # Convertir a formato hexadecimal
        net.add_node(f"{genome}:{community}", label=f"{genome}:{community}", size=size, color=hex_color)

    # Añadir aristas representando relaciones homólogas
    for (genome, community), homologs in homologous_communities.items():
        for homolog_genome, homolog_community in homologs:
            net.add_edge(f"{genome}:{community}", f"{homolog_genome}:{homolog_community}")



    # Mostrar el gráfico
    net.show("homologous_communities.html", notebook=False)

def calculate_centrality_measures(network):
    degree_centrality = nx.degree_centrality(network)
    closeness_centrality = nx.closeness_centrality(network)
    betweenness_centrality = nx.betweenness_centrality(network)
    return {
        'degree': degree_centrality,
        'closeness': closeness_centrality,
        'betweenness': betweenness_centrality
    }

def calculate_additional_metrics(network, nodes):
    try:
        eigen_centralities = nx.eigenvector_centrality(network)
        avg_eigen = np.mean(list(eigen_centralities.values()))
    except nx.PowerIterationFailedConvergence:
        avg_eigen = np.nan

    try:
        katz_centralities = nx.katz_centrality(network, max_iter=5000, alpha=0.005)
        avg_katz = np.mean(list(katz_centralities.values()))
    except nx.PowerIterationFailedConvergence:
        avg_katz = np.nan

    density = nx.density(network)
    clustering_coefficient = nx.average_clustering(network)
    avg_path_length = nx.average_shortest_path_length(network) if nx.is_connected(network) else float('inf')
    transitivity = nx.transitivity(network)
    communities = [set(nodes)]
    modularity = nx.algorithms.community.modularity(network, communities)

    return {
        "density": density,
        "clustering_coefficient": clustering_coefficient,
        "average_path_length": avg_path_length,
        "average_eigen_centralities": avg_eigen,
        "average_katz_centralities": avg_katz,
        "transitivity": transitivity,
        "modularity": modularity
    }

def get_community_info(genome, community, labeled_communities, network):
    nodes = labeled_communities[genome][community]
    element_counts = {element.split(";")[1]: 0 for element in nodes}
    for node in nodes:
        element = node.split(";")[1]
        element_counts[element] += 1

    subgraph = network.subgraph(nodes)
    centralities = calculate_centrality_measures(subgraph)
    additional_metrics = calculate_additional_metrics(subgraph, nodes)

    avg_degree = np.mean(list(centralities['degree'].values()))
    avg_closeness = np.mean(list(centralities['closeness'].values()))
    avg_betweenness = np.mean(list(centralities['betweenness'].values()))

    return [genome, community, len(nodes), element_counts, avg_degree, avg_closeness, avg_betweenness,
            additional_metrics.get('density', np.nan), additional_metrics.get('clustering_coefficient', np.nan),
            additional_metrics.get('average_path_length', np.nan), additional_metrics.get('average_eigen_centralities', np.nan),
            additional_metrics.get('average_katz_centralities', np.nan), additional_metrics.get('transitivity', np.nan),
            additional_metrics.get('modularity', np.nan), ', '.join(nodes)]

def create_homologous_communities_df(labeled_communities, homologous_communities, networks):
    rows = []
    for genome, network in networks.items():
        for community in labeled_communities[genome]:
            row = get_community_info(genome, community, labeled_communities, network)
            homologs = homologous_communities.get((genome, community), [])
            
            if homologs:
                for homolog_genome, homolog_community in homologs:
                    homolog_row = get_community_info(homolog_genome, homolog_community, labeled_communities, networks[homolog_genome])
                    rows.append(row + homolog_row)
            else:
                empty_row = [np.nan] * len(row)
                rows.append(row + empty_row)

    columns = ["Genome1", "Community1", "Size1", "Elements1", "AverageDegree1", "AverageCloseness1", "AverageBetweenness1",
               "Density1", "ClusteringCoefficient1", "AveragePathLength1", "Eigen1", "Katz1","transitivity1","modularity1", "Nodes1",
               "Genome2", "Community2", "Size2", "Elements2", "AverageDegree2", "AverageCloseness2", "AverageBetweenness2",
               "Density2", "ClusteringCoefficient2", "AveragePathLength2", "Eigen2", "Katz2","transitivity2","modularity2", "Nodes2"]
    df = pd.DataFrame(rows, columns=columns)
    return df


def main():
    parser = argparse.ArgumentParser(description="Read blast and generate cluster of duplications")
    parser.add_argument("file", help="Blastp output")
    args = parser.parse_args()

    df = read_and_create(args.file)
    unique_sequences = extract_unique_sequences(df)
    networks = {seq: create_network(df, seq) for seq in unique_sequences}
    pairs = create_pairs_set(df)
    communities_by_genome = {genome: detect_communities(net) for genome, net in networks.items()}
    labeled_communities = label_communities(communities_by_genome)
    homologous_communities = find_homologous_communities(labeled_communities, pairs)
    for genome, communities in labeled_communities.items():
        for community_id, nodes in communities.items():
            print(f"Genoma: {genome}, Comunidad: {community_id}, Nodos: {list(nodes)}")

    # 4. Presentar Agrupaciones Homólogas
    for (genome, community), homologs in homologous_communities.items():
        print(f"En el genoma {genome}, la comunidad {community} es homóloga a las siguientes comunidades en otros genomas:")
        for homolog_genome, homolog_community in homologs:
            print(f"  - {homolog_community} en el genoma {homolog_genome}")



    visualize_homologous_communities_with_pyvis(homologous_communities, labeled_communities)
    #reporte = generate_report(labeled_communities, homologous_communities, networks)
    #print(reporte)
    df = create_homologous_communities_df(labeled_communities, homologous_communities, networks)
    df.to_csv("homologous_communities.tsv", sep="\t", index=False)
if __name__ == "__main__":
    main()