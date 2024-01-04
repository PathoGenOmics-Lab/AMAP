import sys
from Bio import SeqIO
import sourmash
from sourmash.minhash import MinHash
from sklearn.cluster import DBSCAN
import numpy as np
import argparse
import os

def calculate_jaccard_similarity_matrix(bakta_file, h37rv_file):
    """
    Compute the Jaccard similarity matrix for sequences from two FASTA files.

    Parameters:
    - bakta_file (str): Path to the first FASTA file.
    - h37rv_file (str): Path to the second FASTA file.

    Returns:
    - np.ndarray: Jaccard similarity matrix for the combined sequences.
    - list: All sequences from both files.

    Note: Requires SeqIO (BioPython), MinHash (sourmash), and numpy.
    """
    k = 21
    n = 500

    bakta_records = list(SeqIO.parse(bakta_file, "fasta"))
    h37rv_records = list(SeqIO.parse(h37rv_file, "fasta"))

    all_records = bakta_records + h37rv_records

    minhashes = []

    for record in all_records:
        mh = MinHash(ksize=k, n=n, is_protein=False)
        mh.add_sequence(str(record.seq), force=True)
        minhashes.append(mh)

    n_records = len(all_records)
    similarity_matrix = np.zeros((n_records, n_records))

    for i, mh1 in enumerate(minhashes):
        for j, mh2 in enumerate(minhashes):
            if i < j:
                jaccard = mh1.jaccard(mh2)
                similarity_matrix[i, j] = jaccard
                similarity_matrix[j, i] = jaccard

    return similarity_matrix, all_records

def cluster_records(bakta_file, h37rv_file, eps=0.95, min_samples=1):
    """
    Cluster sequences from two FASTA files based on Jaccard similarity.

    Parameters:
    - bakta_file (str): Path to first FASTA file.
    - h37rv_file (str): Path to second FASTA file.
    - eps (float): The maximum distance between two samples for them to be considered as in the same cluster.
    - min_samples (int): The number of samples in a neighborhood for a point to be considered as a core point.

    Returns:
    - ndarray: Cluster labels for each sequence.
    - list: All sequences from both files.

    Note: Uses DBSCAN clustering and requires scikit-learn.
    """
    similarity_matrix, all_records = calculate_jaccard_similarity_matrix(bakta_file, h37rv_file)
    distance_matrix = 1 - similarity_matrix

    clustering = DBSCAN(eps=eps, min_samples=min_samples, metric='precomputed').fit(distance_matrix)
    labels = clustering.labels_

    return labels, all_records

def main(prefix, h37rv_file):
    """
    Generate clusters for sequences and save them to a file based on their Jaccard similarity.

    Parameters:
    - prefix (str): Prefix to generate the bakta file path (`{prefix}-complete.fa`).
    - h37rv_file (str): Path to the FASTA file of h37rv sequences.

    Output:
    - Writes the clusters to `{prefix}.clusters` file.

    Note: Uses SeqIO for FASTA file operations.
    """
    bakta_file = f"{prefix}-complete.fa"
    labels = cluster_records(bakta_file, h37rv_file)

    bakta_records = list(SeqIO.parse(bakta_file, "fasta"))
    h37rv_records = list(SeqIO.parse(h37rv_file, "fasta"))

    all_records = bakta_records + h37rv_records

    with open(f"{prefix}.clusters", "w") as resultado:
        for label in set(labels[0]):
            lista = []
            for i in range(len(labels[0])):
                if labels[0][i] == label:
                    lista.append(all_records[i].id)
            s = f"Cluster {label}:" + ".".join(lista) + "\n"
            resultado.write(s)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Procesar archivos con prefijo variable.")
    parser.add_argument('-n', dest="input_prefix", help="Fasta name")
    parser.add_argument('-p', dest="path", help="Ruta al archivo Bakta multi-fasta.")
    parser.add_argument('-f', dest="h37rv_file", help="Ruta al archivo Bakta multi-fasta.")
    args = parser.parse_args()
    current = os.getcwd()
    os.chdir(args.path)
    h37rv_file = args.h37rv_file
    main(args.input_prefix, current + "/" + h37rv_file)