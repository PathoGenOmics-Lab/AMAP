import sys
from Bio import SeqIO
import sourmash
from sourmash.minhash import MinHash
from sklearn.cluster import DBSCAN
import numpy as np
import argparse

def calculate_jaccard_similarity_matrix(bakta_file, h37rv_file):
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
