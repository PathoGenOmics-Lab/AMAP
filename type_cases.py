#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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

def blastn(in_query: str, in_subject: str) -> pd.DataFrame:
    # Create a BLASTN command line with the specified parameters and output format 6.
    blastn_cline_tab = NcbiblastnCommandline(query=in_query, db=in_subject, outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen', dust="no")
    blastn_out_tab, _ = blastn_cline_tab()
    # Convert the tabular BLASTN output to a list of dictionaries.
    tab_data = [
        dict(zip(['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs', 'qcovhsp', 'qlen', 'slen'], line.split("\t")))
        for line in blastn_out_tab.strip().split("\n")
    ]
    df = pd.DataFrame(tab_data) # Convert the tabular BLASTN results to a DataFrame.
    return df    # Return the tabular BLASTN results.

def extract_sequence(row, unique_seq):
    start = row['Start'] - 1
    stop = row['Stop']
    strand = row['Strand']
    locus_gene = row['Locus-Gene']
    seq = unique_seq[start:stop]
    if strand == '-':
        seq = seq.reverse_complement()
    return SeqIO.SeqRecord(seq, id=locus_gene, description="")

def get_gff(name: str):
    with open(name + '.tsv', 'r') as f:
        lines = f.readlines()
        header_row = next((i for i, line in enumerate(lines) if line.startswith('#Sequence')), None)
    df = pd.read_csv(StringIO(''.join(lines)), sep='\t', header=None, skiprows=header_row, comment='#')
    df.columns = lines[header_row].strip('#').strip().split('\t')

    filtered_df = df.loc[df['Type'] == "cds", ['Locus Tag', 'Gene', 'Strand', 'Start', 'Stop']]
    filtered_df['Gene'] = filtered_df['Gene'].fillna('NA')
    filtered_df['Locus-Gene'] = filtered_df['Locus Tag'] + '-' + filtered_df['Gene']
    filtered_df.to_csv(name + '_bakta.coord', sep='\t', index=False)

    sequences = SeqIO.to_dict(SeqIO.parse(name + '.fna', "fasta"))
    unique_seq = list(sequences.values())[0].seq

    extracted_sequences = filtered_df.apply(extract_sequence, axis=1, args=(unique_seq,))
    SeqIO.write(extracted_sequences.tolist(), name + '_bakta.multifa', "fasta")

def classify_bakta_records(bakta_records, blast_records, output_file='classified_cases.tsv'):
    case1 = set()
    case2 = set()
    case3 = set()
    case4 = set()

    # Convertir bakta_records en un DataFrame de pandas
    bakta_df = pd.DataFrame([{'id': record.id, 'sequence_nt': str(record.seq)} for record in bakta_records])

    # Combinar bakta_df y blast_records usando 'id' y 'qseqid' como claves
    merged_df = pd.merge(bakta_df, blast_records, left_on='id', right_on='qseqid', how='left')

    # Reemplazar NaN por -1
    merged_df.fillna(-1, inplace=True)
    merged_df.to_csv(output_file + '_blastres.tsv', sep='\t', index=False)

    # Case 1: Los registros que no tienen coincidencia en los registros de blast
    case1 = set(merged_df[merged_df['qseqid'] == -1]['id'])

    # Case 2: Los registros con una coincidencia perfecta en los registros de blast
    perfect_match_condition = (merged_df['qlen'] == merged_df['slen']) & \
                              (merged_df['mismatch'] == 0) & \
                              (merged_df['gapopen'] == 0) & \
                              (merged_df['qlen'] == merged_df['length']) & \
                              (merged_df['qcovhsp'] == 100)
    case2 = set(merged_df[perfect_match_condition]['id'])

    # Case 3: Los registros que tienen la misma longitud que su coincidencia en los registros de blast, pero no son una coincidencia perfecta
    same_length_condition = (merged_df['qlen'] == merged_df['slen']) & \
                            ~perfect_match_condition
    case3 = set(merged_df[same_length_condition]['id'])

    # Case 4: Los registros cuya longitud difiere de su coincidencia en los registros de blast y que no estan en el caso 3
    length_diff_condition = (merged_df['qlen'] != merged_df['slen']) & \
                            ~merged_df['id'].isin(case3)
    case4 = set(merged_df[length_diff_condition]['id'])

    # Crear un DataFrame para almacenar los resultados
    results = []
    for case, ids in zip(['case1', 'case2', 'case3', 'case4'], [case1, case2, case3, case4]):
        for id in ids:
            results.append({'id': id, 'case': case})

    results_df = pd.DataFrame(results)

    # Guardar los resultados en un archivo TSV
    results_df.to_csv(output_file + '_classifiedcases.tsv', sep='\t', index=False)
    return case1, case2, case3, case4
 
   
