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
    """
    Reads a FASTA file and extracts only the sequences whose names are in chck_rds.
    Writes the selected sequences to a new FASTA file.

    Parameters:
        chck_rds    (list)  -- A list of sequence IDs to select from the FASTA file.
        fasta_file  (str)   -- The path to the input FASTA file.
        fasta_name  (str)   -- The name of the output FASTA file.

    Returns:
        A list of tuples containing the sequence names and their lengths.
    """
    # Open the FASTA file and select only the sequences with the specified names
    selected_seqs = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        if seq_record.id in chck_rds:
            selected_seqs.append(seq_record)
    # Write the selected sequences to a new FASTA file
    SeqIO.write(selected_seqs, fasta_name, "fasta")

    # Create a list of tuples with the sequence names and their lengths
    seq_lengths = []
    for seq_record in selected_seqs:
        seq_lengths.append((seq_record.id, len(seq_record.seq)))

    return seq_lengths

def create_database(fasta: str, basename: str) -> None:
    """
    Create a BLAST database from a FASTA file using the makeblastdb command line tool.

    Parameters:
        fasta    (str)  -- The path to the input FASTA file.
        basename (str)  -- The basename to use for the BLAST database files.

    Returns:
        None: This function does not return anything, it simply creates the BLAST database files.
    """
    os.system(f'makeblastdb -in {fasta} -dbtype nucl -out {basename}')

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
    lista_111 = []
    lista_122 = []

    for index, row in blast1.iterrows():
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
                    lista_111.append(qseqid) 
                    res = "\t".join([qseqid,sseqid,str(matches),str(slen),str(qlen),str(numb),"1.1.1. remove\n"])
                    logging.info("Resultados:\n%s", "\t".join([qseqid,"1.1.1. remove\n"]))
                    lista_resultados.append(res)
                else:
                    if qseqid not in lista_111:
                        lista_122.append(qseqid)
                        res = "\t".join([qseqid,sseqid,str(matches),str(slen),str(qlen),str(numb),"1.1.2. stay\n"])
                        lista_resultados.append(res)

    for record in case1bakta_records:
        if record.id not in lista_111:
            if record.id not in lista_122:
                lista_122.append(record.id)

    return lista_resultados, lista_111, lista_122

def process_blast112_results(blast112: pd.DataFrame, lista_1121: List[str], lista_122: List[str]) -> Tuple[List[str], List[str]]:
    lista_resultados = []
    lista_1122 = []
    
    processed_qseqids = []

    for index, row in blast112.iterrows():
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
                    res = "\t".join([qseqid, sseqid, str(matches), str(slen), str(qlen), str(numb), "1.1.2.1. remove\n"])
                    lista_resultados.append(res)
                    lista_1121.append(qseqid) 
                else:
                    if qseqid not in lista_1121:
                        res = "\t".join([qseqid, sseqid, str(matches), str(slen), str(qlen), str(numb), "1.1.2.2. stay\n"])
                        lista_resultados.append(res)
                        lista_1122.append(qseqid)

                processed_qseqids.append(qseqid)

    # Clasificar como 1.1.2.2 si no hay resultados de BLAST
    for qseqid in lista_122:
        if qseqid not in processed_qseqids:
            res = "\t".join([qseqid, "No hit", "0", "0", "0", "0", "1.1.2.2. stay\n"])
            lista_resultados.append(res)
            lista_1122.append(qseqid)

    return lista_resultados, lista_1122


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-n', dest="name", help="Fasta name")
    parser.add_argument('-p', dest="path", help="Ruta al archivo Bakta multi-fasta.")
    parser.add_argument('-b', dest="badbakta", help="Ruta al archivo Bakta multi-fasta.")
    parser.add_argument('-r', dest="reference", help="Ruta al archivo Bakta multi-fasta.")
    args = parser.parse_args()
    os.chdir(args.path)
    caso1h37rv = args.badbakta
    referencia = args.reference
    create_database(referencia, "ref_fasta")
    create_database(caso1h37rv, "ref_case1")

    name = args.name
    os.chdir(args.path)
    caso1anot = name + "_case1.bakta.fa"
    case1bakta_records = SeqIO.parse(caso1anot, "fasta")
    blast1 = blastn(caso1anot, "ref_case1")

    lista_resultados, lista_111, lista_122 = process_blast1_results(blast1, case1bakta_records)

    get_seq(lista_122, caso1anot, "case112.fasta")
    blast112 = blastn("case112.fasta", "ref_fasta")

    lista_1121 = []
    lista_resultados_112, lista_1122 = process_blast112_results(blast112, lista_1121, lista_122)

    lista_resultados.extend(lista_resultados_112)
    # Actualizar la lista de resultados
    highest_percentage = {}
    updated_list = []
    
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
    
    for qseqid, entry in highest_percentage.items():
        updated_entry = "\t".join([qseqid, entry["sseqid"], str(entry["matches"]), str(entry["slen"]), str(entry["qlen"]), str(entry["numb"]), entry["category"]])
        updated_list.append(updated_entry)
    print(updated_list)
    with open(name+"_case1.final.tsv","w") as salida1:
        for i in updated_list:
            sentencia = name+"\t"+i
            print(sentencia)
            data = i.split("\t")
            sentencia = data[0]+"\t"+data[6].strip("\n")+"\t0\t0\tNA\n"
            salida1.write(sentencia)
            

if __name__ == "__main__":
    main()
