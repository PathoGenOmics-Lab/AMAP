#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbitblastnCommandline
import subprocess
from tempfile import NamedTemporaryFile
from io import StringIO
import re
import argparse
import glob
from functools import partial
import os
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from typing import Tuple, Dict
from concurrent.futures import ThreadPoolExecutor, as_completed

def run_blast(in_query: str, in_subject: str, outfmt: str) -> str:
    '''Runs blastn with the given parameters and returns the output.'''
    cline = NcbitblastnCommandline(query=in_query, db=in_subject, outfmt=outfmt, seg="no")
    out, _ = cline()
    return out

def run_blast_tab(in_query: str, in_subject: str): 
    '''Runs blastn with the given parameters and returns the output in tabular format.'''
    return run_blast(in_query, in_subject, '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen')

def run_blast_xml(in_query: str, in_subject: str): 
    '''Runs blastn with the given parameters and returns the output in xml format.'''
    return run_blast(in_query, in_subject, '5')

def parse_tab_data(blastn_out_tab: str) -> list:
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs', 'qcovhsp', 'qlen', 'slen']
    return [dict(zip(columns, line.split("\t"))) for line in blastn_out_tab.strip().split("\n") if line]

def parse_xml_data(blastn_out_xml: str) -> dict:
    return {record.query: record for record in NCBIXML.parse(StringIO(blastn_out_xml)) if blastn_out_xml.strip()}

def blastn(in_query: str, in_subject: str) -> pd.DataFrame:
    '''Runs blastn with the given parameters and returns the output in a pandas DataFrame.'''
    
    with ThreadPoolExecutor(max_workers=2) as executor:
        future_tab = executor.submit(run_blast_tab, in_query, in_subject)
        future_xml = executor.submit(run_blast_xml, in_query, in_subject)
        
        blastn_out_tab = future_tab.result()
        blastn_out_xml = future_xml.result()

    tab_data = parse_tab_data(blastn_out_tab)
    blast_records = parse_xml_data(blastn_out_xml)

    combined_data = []
    for tab in tab_data:
        xml_record = blast_records.get(tab['qseqid'])
        if xml_record and any(tab['sseqid'] in alignment.hit_def for alignment in xml_record.alignments):
            combined_data.append({**tab,
                                  "xml_data": xml_record,
                                  "subject_name": xml_record.descriptions[0].title.split()[1] if xml_record.descriptions else None})
    
    return pd.DataFrame(combined_data)

def get_seq(chck_rds, fasta_file: str, fasta_name: str):
    '''Extracts a sequence from a fasta file, translates it, and writes it to a new fasta file.'''
    def selected_seqs_gen():
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            if chck_rds in seq_record.id:
                protein_seq = seq_record.seq.translate()
                protein_record = SeqRecord(protein_seq, id = seq_record.id, description = "")
                yield protein_record
                break
    SeqIO.write(selected_seqs_gen(), fasta_name, "fasta")

def parse_argparse() -> argparse.Namespace:
    '''Parse the command line arguments and return the parsed arguments.'''
    parser = argparse.ArgumentParser(description = 'A description of your script.')
    parser.add_argument('-p', '--prefix', required = True, help = 'Path to the prefix of the files')
    args = parser.parse_args()
    return args

def clean_up_files(prefix: str):
    '''Remove all files matching the pattern '*database*' under the given prefix.'''
    file_pattern = f"{prefix}*database*"
    files = glob.glob(file_pattern)
    for f in files:
        try:
            os.remove(f)
            print(f'Successfully removed file {f}')
        except OSError as e:
            print(f'Failed to remove file {f}. Reason: {e}')

def create_database(fasta: str, basename: str) -> None:
    '''Creates a blast database from a fasta file.'''
    subprocess.run(['makeblastdb', '-in', fasta, '-dbtype', 'nucl', '-out', basename])

def process_case_genes(gene, coordinates_set, output_file):
    gene_info = gene.split(";")
    start, end = map(int, gene_info[1].split(","))
    to_write = f"{gene_info[0]}\t{start}\t{end}\t{gene_info[2]}"
    output_file.write(to_write + "\n")
    coordinates_set.add((start, end))
    return coordinates_set

def process_non_case_genes(gene, coordinates_set, output_file, prefix):
    gen_name, gen_range = gene.split("_")
    start, end = map(int, gen_range.split("-"))
    gene_details = [gen_name, abs(end - start), gene]
    return process_non_matching_genes([gene_details], output_file, coordinates_set, prefix)

def cluster_1_genes(line, coordinates_set, output_file, prefix):
    _, gene_data = line.strip().split(":")
    genes = gene_data.split(".")
    for gene in genes:
        if "case" in gene:
            if "case1" in gene:
                coordinates_set = process_case_genes(gene, coordinates_set, output_file)
        else:
            coordinates_set = process_non_case_genes(gene, coordinates_set, output_file, prefix)
    coordinates_set = set(sorted(coordinates_set))
    return coordinates_set

def process_coordinates(coords, bakta_coords, file_writer, gene_id):
    ordered_coords = (min(int(coords[0]), int(coords[1])), max(int(coords[0]), int(coords[1])))
    noorder_len = int(coords[1]) - int(coords[0])
    len_difference = ordered_coords[1] - ordered_coords[0]
    length_indicator = "+" if noorder_len >= 0 else "-"
    if len_difference > 100:
        if ordered_coords not in bakta_coords:
            file_writer.write(f"{gene_id}\t{ordered_coords[0]}\t{ordered_coords[1]}\t{length_indicator}\n")
            bakta_coords.add(ordered_coords)
            return True
    return False

def process_gene(gene_data, bakta_coords, file_writer):
    gene_id, coordinates_list = gene_data
    for coords in coordinates_list:
        if process_coordinates(coords, bakta_coords, file_writer, gene_id):
            break

def process_non_matching_genes(non_matching_gen: list, file_writer, bakta_coords: set, prefix) -> set:
    fixed_process_element_with_blastn = partial(process_element_with_blastn, prefix=prefix)
    num_workers = 4 

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(fixed_process_element_with_blastn, gen): gen for gen in non_matching_gen}
        
        for future in as_completed(futures):
            try:
                result = future.result()
                for gene_data in result.items():
                    process_gene(gene_data, bakta_coords, file_writer)
            except Exception as e:
                print(f"Exception occurred: {e}")

    return bakta_coords

def process_element_with_blastn(element: Tuple, prefix: str) -> Dict:
    '''Processes a single element using blastn and returns results'''
    name = os.path.basename(prefix) + "database"
    pid = os.getpid()

    with NamedTemporaryFile(suffix=f".{pid}.fasta") as query_file:
        get_seq(element[0], "/storage/PGO/results/mtb/proj_anc2023/bakta/last_step/nt_correctedH37Rv.fa", query_file.name)        
        blast_out = blastn(query_file.name, name)

    return process_blast_results(blast_out)

def calculate_matches(alignment, qlen):
    exact_matches = sum(a == b for a, b in zip(alignment.query, alignment.match))
    return (exact_matches / float(qlen)) * 100

def is_valid_alignment(alignment, start, end, qlen):
    orientation = "+" if end - start >= 0 else "-"

    #if alignment.match[-1] == '*' or alignment.sbjct[-1] != '*':
    if alignment.match[-1] == '*':
        return True, (start, end, qlen)
    else:
        longitud, coordenadas = extract_amino_acid_sequence_with_stop_codon(start, end, orientation)
        calculo = (min(int(longitud), int(qlen)) / max(int(longitud), int(qlen))) * 100
        if calculo >= 80:
            return True, (start, coordenadas[0], coordenadas[0])
    return False, None

def process_blast_results(blast_out):
    blast_results = defaultdict(list)
    for blast_result in blast_out.itertuples(index=False):
        pident = float(blast_result.pident)
        if pident <= 70:
            continue
        alignment = blast_result.xml_data.alignments[0].hsps[0]
        matches = calculate_matches(alignment, blast_result.qlen)
        if matches < 80:
            continue
        start, end = map(int, (blast_result.sstart, blast_result.send))
        is_valid, result = is_valid_alignment(alignment, start, end, blast_result.qlen)
        if is_valid:
            blast_results[blast_result.qseqid].append(result)
    return blast_results

def translate_sequence(sequence, start, end, step):
    amino_acid_sequence = ""
    for i in range(start, end, step):
        next_codon = sequence[i:i+3]
        if step < 0:
            next_codon = next_codon.reverse_complement()
        next_amino_acid = next_codon.translate(table="Bacterial", to_stop=False)
        amino_acid_sequence += next_amino_acid
        if next_amino_acid == "*":
            return amino_acid_sequence, i
    return amino_acid_sequence, i

def extract_amino_acid_sequence_with_stop_codon(start, end, orientation):
    start = start - 1
    global SEQ
    fasta_file = SEQ
    with open(SEQ, "r") as fasta_file:
        sequence = [seq_record.seq for seq_record in SeqIO.parse(fasta_file, "fasta")][0]
    
    seq = sequence[start:end]
    if orientation == "-":
        seq = seq.reverse_complement()
    
    amino_acid_sequence = seq.translate(table="Bacterial", to_stop=False)

    if orientation == "+":
        step = 3
        aa_seq, stop_codon_index = translate_sequence(sequence, end, len(sequence), step)
    else: # orientation == "-"
        step = -3
        aa_seq, stop_codon_index = translate_sequence(sequence, start, -1, step)
        
    return len(amino_acid_sequence + aa_seq), (start, stop_codon_index + step)

def count_cases(genes):
    '''Counts the number of genes in each case and the number of Rv genes.'''
    g_genes = genes.split(".")
    cases = ["case1", "case2", "case3", "case4"]
    counts = {case: genes.count(case) for case in cases}
    rvs = genes.count("Rv")
    diff = rvs - sum(counts.values())
    return rvs, g_genes, diff, counts["case1"], counts["case2"], counts["case3"], counts["case4"]

def process_element(element):
    '''Processes an element from the cluster file and returns the element name, length and the element itself.'''
    if "case" in element:
        length = int(element.split(";")[1].split(",")[1]) - int(element.split(";")[1].split(",")[0])
        return [element.split(";")[0], length, element]
    else:
        return [element.split("_")[0], int(element.split("_")[1].split("-")[1]) - int(element.split("_")[1].split("-")[0]), element]
    
def process_coordinates_and_write(element_info, rv_element, coordinates_list, output):
    """Process coordinates and write to output if the coordinates are not in the list."""
    element_info_split = element_info[2].split(";")
    real_coordinates = tuple(map(int, element_info_split[1].split(",")))

    if real_coordinates not in coordinates_list:

        coordinates_list.add(real_coordinates)
        output.write(f"{element_info[0]}\t{real_coordinates[0]}\t{real_coordinates[1]}\t{element_info_split[2]}\t{rv_element}\n")
    return coordinates_list

def match_elements(l_bak, l_gen):
    '''Matches the elements from the cluster file with the elements from the blastn output.'''
    non_matching_bak = [bak for bak in l_bak if not any(bak[1] == gen[1] for gen in l_gen)]
    non_matching_gen = [gen for gen in l_gen if not any(gen[1] == bak[1] for bak in l_bak)]
    matching_elements = [(bak, gen) for bak in l_bak for gen in l_gen if bak[1] == gen[1]]
    return non_matching_bak, non_matching_gen, matching_elements