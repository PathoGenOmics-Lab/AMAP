#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import multiprocessing
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import pandas as pd
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

    bakta_df = pd.DataFrame([{'id': record.id, 'sequence_nt': str(record.seq)} for record in bakta_records])

    # Combinar bakta_df y blast_records usando 'id' y 'qseqid' como claves
    merged_df = pd.merge(bakta_df, blast_records, left_on='id', right_on='qseqid', how='left')

    # Reemplazar NaN por -1
    merged_df.fillna(-1, inplace=True)
    
    # Convertir las columnas relevantes a int
    merged_df['qlen'] = merged_df['qlen'].astype(int)
    merged_df['slen'] = merged_df['slen'].astype(int)
    merged_df['mismatch'] = merged_df['mismatch'].astype(int)
    merged_df['gapopen'] = merged_df['gapopen'].astype(int)
    merged_df['length'] = merged_df['length'].astype(int)
    merged_df['qcovhsp'] = merged_df['qcovhsp'].astype(int)

    merged_df.to_csv(output_file + '_blastres.tsv', sep='\t', index=False)
    
    # Clasificar cada id
    for id in bakta_df['id']:
        in_case3 = False
        for _, row in merged_df[merged_df['id'] == id].iterrows():
            if row['qseqid'] == -1:
                case1.add(id)
            elif row['qlen'] == row['slen'] and row['mismatch'] == 0 and row['gapopen'] == 0 and row['qlen'] == row['length'] and row['qcovhsp'] == 100:
                case2.add(id)
                break
            elif row['qlen'] == row['slen']:
                case3.add(id)
                in_case3 = True
            elif not in_case3 and (row['qlen'] > row['slen'] or row['qlen'] < row['slen']):
                case4.add(id)

    # Eliminar los elementos de case4 que también están en case3
    case3 = case3 - case4

    # Crear un DataFrame para almacenar los resultados
    results = []
    for case, ids in zip(['case1', 'case2', 'case3', 'case4'], [case1, case2, case3, case4]):
        for id in ids:
            results.append({'id': id, 'case': case})

    results_df = pd.DataFrame(results)

    # Guardar los resultados en un archivo TSV
    results_df.to_csv(output_file + '_classifiedcases.tsv', sep='\t', index=False)

    return case1, case2, case3, case4
 
def process_genome(genome: Tuple[str, str], current_working_directory: str, blast_db_multi: str) -> None:
    name, path = genome[0], genome[1]
    os.chdir(path)
    get_gff(name)
    bakta_records = SeqIO.parse(f"{name}_bakta.multifa", "fasta")
    blast_db_path = os.path.join(current_working_directory, blast_db_multi)
    blast_records = blastn(f"{name}_bakta.multifa", blast_db_path)
    cases = classify_bakta_records(bakta_records, blast_records, name)
    for case_num, case in enumerate(cases, start=1):
        get_seq(case, f"{name}_bakta.multifa", f"{name}_case{case_num}.bakta.fa")

def main():
    parser = argparse.ArgumentParser(description="Clasifica los registros de Bakta y BLAST y genera los archivos fasta respectivos.")
    parser.add_argument('-t', '--tsv_file', help="Ruta al archivo Bakta TSV que contiene los genomas a procesar.")
    parser.add_argument('-r', '--reference_multi_fa', help="Ruta al archivo multi-fasta de referencia para crear la base de datos de BLAST.")
    parser.add_argument('-n', '--blast_db_multi', help="Nombre de la base de datos de BLAST a crear.")
    parser.add_argument('-p', '--num_processes', default=multiprocessing.cpu_count(), type=int, help="Número de procesos a utilizar para el análisis en paralelo. Por defecto, se usa el número de núcleos de la CPU.")
    args = parser.parse_args()

    create_database(args.reference_multi_fa, args.blast_db_multi)
    tsv_list = get_tsv(args.tsv_file)
    current_working_directory = os.getcwd()

    with multiprocessing.Pool(args.num_processes) as pool:
        pool.starmap(process_genome, [(genome, current_working_directory, args.blast_db_multi) for genome in tsv_list])

if __name__ == "__main__":
    main()
