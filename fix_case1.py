#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
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
    """
    Function to perform a BLASTN search and return the results in a DataFrame with both tabular (format 6) and XML (format 5) outputs.

    Arguments:
        in_query        -- string with the path or filename of the query.
        in_subject      -- string with the path or filename of the subject database.

    Returns:
        df              -- DataFrame containing the combined BLASTN results in tabular and XML formats.
    """
    # Create a BLASTN command line with the specified parameters and output format 6.
    blastn_cline_tab = NcbiblastnCommandline(query=in_query, db=in_subject, outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen', dust="no")
    blastn_out_tab, _ = blastn_cline_tab()

    # Create a BLASTN command line with the specified parameters and output format 5.
    blastn_cline_xml = NcbiblastnCommandline(query=in_query, db=in_subject, outfmt='5', dust="no")
    blastn_out_xml, _ = blastn_cline_xml()

    # Convert the tabular BLASTN output to a list of dictionaries.
    tab_data = [
        dict(zip(['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs', 'qcovhsp', 'qlen', 'slen'], line.split("\t")))
        for line in blastn_out_tab.strip().split("\n")
    ] if len(blastn_out_tab.strip()) > 0 else []

    # Convert the XML BLASTN output to a list of BlastRecords.
    blastn_out_xml = StringIO(blastn_out_xml)

    if blastn_out_xml.getvalue().strip():
        blast_records = list(NCBIXML.parse(blastn_out_xml))
    else:
        blast_records = []

    # Combine the tabular and XML BLASTN outputs.
    combined_results = []

    for tab in tab_data:
        # Find the XML data that corresponds to the current tabular data.
        query_id = tab['qseqid']
        subject_id = tab['sseqid']

        # Find the BlastRecord corresponding to the current tabular data.
        xml_record = next((record for record in blast_records if record.query == query_id and any(subject_id in alignment.hit_def for alignment in record.alignments)), None)

        if xml_record is not None:
            combined_result = {
                "tab_data": tab,
                "xml_data": xml_record,
                "subject_name": xml_record.descriptions[0].title.split()[1] if len(xml_record.descriptions) > 0 else None
            }
            combined_results.append(combined_result)

    # Convert the combined BLASTN results to a DataFrame.
    df = pd.DataFrame.from_records(
        [
            {
                **{"" + k: v for k, v in result["tab_data"].items()},
                "xml_data": result["xml_data"],
                "subject_name": result["subject_name"]
            }
            for result in combined_results
        ]
    )

    # Return the combined BLAST
    return df

def process_blast1_results(blast1: pd.DataFrame, case1bakta_records) -> Tuple[List[str], List[str], List[str]]:
    lista_resultados = []
    lista_111 = set()
    lista_122 = set()

    for _, row in blast1.iterrows():
        qseqid = row['qseqid']
        sseqid = row['sseqid']
        slen = int(row['slen'])
        qlen = int(row['qlen'])
        xml_data = row['xml_data']
        for alignment in xml_data.alignments:
            for hsp in alignment.hsps:
                matches = hsp.match.count('|')
                numb = (matches/slen)*100
                if numb >= 80:
                    lista_111.add(qseqid) 
                    res = f"{qseqid}\t{sseqid}\t{matches}\t{slen}\t{qlen}\t{numb}\t1.1.1. remove\n"
                    logging.info(f"Resultados:\n{qseqid}\t1.1.1. remove\n")
                    lista_resultados.append(res)
                elif qseqid not in lista_111:
                    lista_122.add(qseqid)
                    res = f"{qseqid}\t{sseqid}\t{matches}\t{slen}\t{qlen}\t{numb}\t1.1.2. stay\n"
                    lista_resultados.append(res)

    for record in case1bakta_records:
        if record.id not in lista_111 and record.id not in lista_122:
            lista_122.add(record.id)

    return lista_resultados, list(lista_111), list(lista_122)

def process_blast112_results(blast112: pd.DataFrame, lista_1121: List[str], lista_122: List[str]) -> Tuple[List[str], List[str]]:
    lista_resultados = []
    lista_1122 = set()
    
    processed_qseqids = set()

    for _, row in blast112.iterrows():
        qseqid = row['qseqid']
        sseqid = row['sseqid']
        slen = int(row['slen'])
        qlen = int(row['qlen'])
        xml_data = row['xml_data']
        for alignment in xml_data.alignments:
            for hsp in alignment.hsps:
                matches = hsp.match.count('|')
                numb = (matches/qlen)*100
                if numb >= 80:
                    res = f"{qseqid}\t{sseqid}\t{matches}\t{slen}\t{qlen}\t{numb}\t1.1.2.1. remove\n"
                    lista_resultados.append(res)
                    lista_1121.append(qseqid) 
                elif qseqid not in lista_1121:
                    res = f"{qseqid}\t{sseqid}\t{matches}\t{slen}\t{qlen}\t{numb}\t1.1.2.2. stay\n"
                    lista_resultados.append(res)
                    lista_1122.add(qseqid)

                processed_qseqids.add(qseqid)

    # Clasificar como 1.1.2.2 si no hay resultados de BLAST
    for qseqid in lista_122:
        if qseqid not in processed_qseqids:
            res = f"{qseqid}\tNo hit\t0\t0\t0\t0\t1.1.2.2. stay\n"
            lista_resultados.append(res)
            lista_1122.add(qseqid)

    return lista_resultados, list(lista_1122)

def main():
    """
    Esta es la función principal que ejecuta el análisis.
    """
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-t', dest="tsv_file", help="Path to the Bakta multi-fasta file.")
    args = parser.parse_args()

    caso1h37rv = "case1_H37Rv.fasta"
    referencia = "reference.fa"
    create_database(referencia, "ref_fasta")
    create_database(caso1h37rv, "ref_case1")

    tsv_list = get_tsv(args.tsv_file)
    current_working_directory = os.path.abspath('.')

    with open("resultados.txt","w") as salida:
        for genome in tsv_list:
            name, path = genome
            os.chdir(path)
            caso1anot = f"{name}_case1.bakta.fa"
            case1bakta_records = SeqIO.parse(caso1anot, "fasta")
            blast1 = blastn(caso1anot, os.path.join(current_working_directory, "ref_case1"))

            lista_resultados, lista_111, lista_122 = process_blast1_results(blast1, case1bakta_records)

            get_seq(lista_122, caso1anot, "case112.fasta")
            blast112 = blastn("case112.fasta", os.path.join(current_working_directory, "ref_fasta"))

            lista_1121 = []
            lista_resultados_112, lista_1122 = process_blast112_results(blast112, lista_1121, lista_122)

            lista_resultados.extend(lista_resultados_112)
            
            highest_percentage = process_results(lista_resultados)
        
            updated_list = [format_entry(qseqid, entry) for qseqid, entry in highest_percentage.items()]
            
            with open(f"{name}_case1.final.tsv","w") as salida1:
                for i in updated_list:
                    sentencia = f"{name}\t{i}"
                    data = i.split("\t")
                    salida.write(sentencia)
                    sentencia = f"{data[0]}\t{data[6].strip()}\t0\t0\tNA\n"
                    salida1.write(sentencia)

def process_results(lista_resultados):
    """
    Procesa la lista de resultados y devuelve el porcentaje más alto para cada qseqid.
    """
    highest_percentage = {}
    for item in lista_resultados:
        data = item.split("\t")
        qseqid = data[0]
        sseqid = data[1]
        matches = int(data[2])
        slen = int(data[3])
        qlen = int(data[4])
        numb = float(data[5])
        category = data[6]
        
        if qseqid in highest_percentage:
            current_category = highest_percentage[qseqid]["category"]
            if category.startswith("remove") and not current_category.startswith("remove"):
                highest_percentage[qseqid] = {"sseqid": sseqid, "matches": matches, "slen": slen, "qlen": qlen, "numb": numb, "category": category}
            elif numb > highest_percentage[qseqid]["numb"]:
                highest_percentage[qseqid] = {"sseqid": sseqid, "matches": matches, "slen": slen, "qlen": qlen, "numb": numb, "category": category}
        else:
            highest_percentage[qseqid] = {"sseqid": sseqid, "matches": matches, "slen": slen, "qlen": qlen, "numb": numb, "category": category}
    return highest_percentage

def format_entry(qseqid, entry):
    """
    Formatea una entrada de la lista de resultados.
    """
    return "\t".join([qseqid, entry["sseqid"], str(entry["matches"]), str(entry["slen"]), str(entry["qlen"]), str(entry["numb"]), entry["category"]])

if __name__ == "__main__":
    main()

