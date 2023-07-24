#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
import pandas as pd
from io import StringIO
from pathlib import Path

def get_seq(chck_rds: list, fasta_file: str, fasta_name: str):
    chck_rds_set = set(chck_rds)  # Faster searching
    selected_seqs = [seq_record for seq_record in SeqIO.parse(fasta_file, "fasta") if seq_record.id in chck_rds_set]
    SeqIO.write(selected_seqs, fasta_name, "fasta") # Write the selected sequences to a fasta file

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

def get_gff(name: str, current):
    
    with open(current+"/out/bakta_output/" + name + ".bakta/" + name + '.tsv', 'r') as f:
        lines = f.readlines()
        header_row = next((i for i, line in enumerate(lines) if line.startswith('#Sequence')), None)
    df = pd.read_csv(StringIO(''.join(lines)), sep='\t', header=None, skiprows=header_row, comment='#')
    df.columns = lines[header_row].strip('#').strip().split('\t')

    filtered_df = df.loc[df['Type'] == "cds", ['Locus Tag', 'Gene', 'Strand', 'Start', 'Stop']]
    filtered_df['Gene'] = filtered_df['Gene'].fillna('NA')
    filtered_df['Locus-Gene'] = filtered_df['Locus Tag'] + '-' + filtered_df['Gene']
    filtered_df.to_csv(name + '_bakta.coord', sep='\t', index=False)

    sequences = SeqIO.to_dict(SeqIO.parse(current+"/out/bakta_output/" + name + ".bakta/" + name + '.fna', "fasta"))
    unique_seq = list(sequences.values())[0].seq

    extracted_sequences = filtered_df.apply(extract_sequence, axis=1, args=(unique_seq,))
    SeqIO.write(extracted_sequences.tolist(), name + '_bakta.multifa', "fasta")

def is_case1(row):
    return row['qseqid'] == -1

def is_case2(row):
    return (row['qlen'] == row['slen'] and 
            row['mismatch'] == 0 and 
            row['gapopen'] == 0 and 
            row['qlen'] == row['length'] and 
            row['qcovhsp'] == 100)

def is_case3(row):
    return row['qlen'] == row['slen']

def is_case4(row, in_case3):
    return not in_case3 and (row['qlen'] > row['slen'] or row['qlen'] < row['slen'])

def read_bakta_records(bakta_records):
    return pd.DataFrame([{'id': record.id, 'sequence_nt': str(record.seq)} for record in bakta_records])

def preprocess_merged_df(merged_df):
    # Replace NaN with -1
    merged_df.fillna(-1, inplace=True)
    
    # Convert relevant columns to int
    for column in ['qlen', 'slen', 'mismatch', 'gapopen', 'length', 'qcovhsp']:
        merged_df[column] = merged_df[column].astype(int)
    
    return merged_df

def write_to_csv(df, output_file):
    df.to_csv(output_file, sep='\t', index=False)

def get_cases_for_id(merged_df, id):
    case1, case2, case3, case4 = set(), set(), set(), set()
    in_case3 = False
    
    for _, row in merged_df[merged_df['id'] == id].iterrows():
        if is_case1(row):
            case1.add(id)
        elif is_case2(row):
            case2.add(id)
            break
        elif is_case3(row):
            case3.add(id)
            in_case3 = True
        elif is_case4(row, in_case3):
            case4.add(id)
    
    return case1, case2, case3 - case4, case4 # Subtract case4 from case3 directly here

def classify_bakta_records(bakta_records, blast_records, output_file='classified_cases.tsv'):
    case1, case2, case3, case4 = set(), set(), set(), set()
    
    bakta_df = read_bakta_records(bakta_records)

    # Combine bakta_df and blast_records using 'id' and 'qseqid' as keys
    merged_df = pd.merge(bakta_df, blast_records, left_on='id', right_on='qseqid', how='left')
    merged_df = preprocess_merged_df(merged_df)
    
    write_to_csv(merged_df, output_file + '_blastres.tsv')
    
    # Classify each id
    for id in bakta_df['id']:
        case1_id, case2_id, case3_id, case4_id = get_cases_for_id(merged_df, id)
        case1 |= case1_id
        case2 |= case2_id
        case3 |= case3_id
        case4 |= case4_id

    # Create a DataFrame to store the results
    results = []
    for case, ids in zip(['case1', 'case2', 'case3', 'case4'], [case1, case2, case3, case4]):
        for id in ids:
            results.append({'id': id, 'case': case})

    results_df = pd.DataFrame(results)

    # Save the results to a TSV file
    write_to_csv(results_df, output_file + '_classifiedcases.tsv')

    return case1, case2, case3, case4

def process_genome(name, blast_db_multi: str, current) -> None:
    get_gff(name, current)
    bakta_records = SeqIO.parse(f"{name}_bakta.multifa", "fasta")
    blast_records = blastn(f"{name}_bakta.multifa", blast_db_multi)
    cases = classify_bakta_records(bakta_records, blast_records, name)
    for case_num, case in enumerate(cases, start=1):
        get_seq(case, f"{name}_bakta.multifa", f"{name}_case{case_num}.bakta.fa")

def main():
    parser = argparse.ArgumentParser(description="Classifies Bakta and BLAST records and generates the respective fasta files.")
    parser.add_argument('-f', '--fasta', help="Path to the Bakta TSV file containing genomes to process.")
    parser.add_argument('-p', '--path', help="Path to working directory to genome of interest.")
    parser.add_argument('-b', '--blast', help="Blast path database.")

    args = parser.parse_args()
    current = os.getcwd()
    path = args.path
    Path(path).mkdir(parents=True, exist_ok=True)
    os.chdir(path)
    blast_db_multi = args.blast

    process_genome(args.fasta, current + "/" + blast_db_multi, current)

if __name__ == "__main__":
    main()