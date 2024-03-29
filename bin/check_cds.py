#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import os

def find_next_stop_codon(fasta_seq, start_index, strand):
    """Find the next stop codon in the DNA sequence using Biopython."""
    sequence = fasta_seq.seq

    if strand == "+":
        # Para la hebra positiva, recorrer la secuencia hacia adelante
        for i in range(start_index, len(sequence), 3):
            codon = sequence[i:i + 3]
            if len(codon) == 3 and codon.translate(table=1) == "*":
                return i
    else:
        # Para la hebra negativa, recorrer la secuencia hacia atrás
        for i in range(start_index, -1, -3):
            codon = sequence[i:i + 3].reverse_complement()
            if len(codon) == 3 and codon.translate(table=1) == "*":
                return i - 2

    return None

def get_goodseq(sequence, start, end, strand):
    """Get the good sequence for the gene with the stop codon"""
    if strand == "+":
        end = find_next_stop_codon(sequence, end, strand) + 3 
    else:
        start = find_next_stop_codon(sequence, start, strand) + 2
    return start, end, strand

def read_to_pandas(file):
    df = pd.read_csv(file, sep='\t', header=None)
    df.columns = ["id", "strand", "start", "end", "rvs"]
    return df

def extract_and_translate_cds(fasta_file, annotation_df):
    # Carga la secuencia única del archivo FASTA
    fasta_seq = next(SeqIO.parse(fasta_file, "fasta"))
    fix_cds_stop = []
    # Extracción y traducción de CDS
    cds_amino_acids = {}
    cds_nts = {}
    for index, row in annotation_df.iterrows():
        start = row['start'] - 1  # Python es 0-indexado
        end = row['end']
        strand = row['strand']
        ide = row['id']
        rvs = row['rvs']
        
        sequence = fasta_seq.seq[start:end]
        if strand == '-':
            sequence = sequence.reverse_complement()
            
        translated_seq = sequence.translate()
        if translated_seq[-1] == "*":
            cds_id = f"{ide};{strand};{start}-{end}"
            cds_amino_acids[cds_id] = translated_seq
            cds_nts[cds_id] = sequence
            
        else:
            if strand == "-":
                start, end, strand = get_goodseq(fasta_seq, start, end, strand)
                sequence_nt = fasta_seq.seq[start:end].reverse_complement()
            else:
                end = find_next_stop_codon(fasta_seq, end, strand)
                start, end, strand = get_goodseq(fasta_seq, start, end, strand)
                sequence_nt = fasta_seq.seq[start:end]
            cds_amino_acids[cds_id] = sequence_nt.translate()
            cds_nts[cds_id] = sequence_nt
            # Actualiza el DataFrame
            annotation_df.at[index, 'start'] = start + 1  # Convertir de nuevo a 1-indexado
            annotation_df.at[index, 'end'] = end   
                     
    return cds_amino_acids, cds_nts, annotation_df

def main():
    parser = argparse.ArgumentParser(description="Extraer y traducir CDS de un archivo FASTA a aminoácidos")
    parser.add_argument("fasta_file", help="Archivo FASTA que contiene la secuencia")
    parser.add_argument("annotation_file", help="Archivo de anotaciones en formato tabulado")
    args = parser.parse_args()
    base_name = args.annotation_file.replace('_corrected.tsv','')
    annotation_df = read_to_pandas(args.annotation_file)
    cds_amino_acids, cds_nts, annotation_df = extract_and_translate_cds(args.fasta_file, annotation_df)
    
    annotation_df.to_csv(base_name+'.fixedcds.tsv', sep='\t', index=False)

    # Imprimir las secuencias de aminoácidos
    with open(base_name+'.faa', 'w') as file:
        for id, amino_acid_seq in cds_amino_acids.items():
            file.write(f">{id}\n{amino_acid_seq}\n")
    with open(base_name+'.fna', 'w') as file:
        for id, amino_acid_seq in cds_nts.items():
            file.write(f">{id}\n{amino_acid_seq}\n")
if __name__ == "__main__":
    main()
    #python3 check_cds.py H37Rv.fasta _ancestor_corrected.tsv
