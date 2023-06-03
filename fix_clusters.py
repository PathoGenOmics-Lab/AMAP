#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from concurrent.futures import ThreadPoolExecutor
import subprocess
from tempfile import NamedTemporaryFile
from io import StringIO
import re
import argparse
import glob

def get_seq(chck_rds, fasta_file: str, fasta_name: str):
    def selected_seqs_gen():
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            if chck_rds in seq_record.id:
                yield seq_record
                break
    SeqIO.write(selected_seqs_gen(), fasta_name, "fasta")


def create_database(fasta: str, basename: str) -> None:
    subprocess.run(['makeblastdb', '-in', fasta, '-dbtype', 'nucl', '-out', basename])

def run_blast(in_query: str, in_subject: str, outfmt: str) -> str:
    cline = NcbiblastnCommandline(query=in_query, db=in_subject, outfmt=outfmt, dust="no")
    out, _ = cline()
    return out

from concurrent.futures import ThreadPoolExecutor
from io import StringIO
from Bio.Blast import NCBIXML
import pandas as pd

def run_blast_tab(in_query: str, in_subject: str): 
    return run_blast(in_query, in_subject, '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen')

def run_blast_xml(in_query: str, in_subject: str): 
    return run_blast(in_query, in_subject, '5')

def blastn(in_query: str, in_subject: str) -> pd.DataFrame:
    with ThreadPoolExecutor(max_workers=2) as executor:
        future_tab = executor.submit(run_blast_tab, in_query, in_subject)
        future_xml = executor.submit(run_blast_xml, in_query, in_subject)
        
        blastn_out_tab = future_tab.result()
        blastn_out_xml = future_xml.result()

    tab_data = [
        dict(zip(['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs', 'qcovhsp', 'qlen', 'slen'], line.split("\t")))
        for line in blastn_out_tab.strip().split("\n")
    ] if len(blastn_out_tab.strip()) > 0 else []

    blast_records = {record.query: record for record in NCBIXML.parse(StringIO(blastn_out_xml)) if blastn_out_xml.strip()} or {}

    combined_results = [
        {**tab,
            "xml_data": xml_record,
            "subject_name": xml_record.descriptions[0].title.split()[1] if xml_record.descriptions else None}
        for tab in tab_data
        if (xml_record := blast_records.get(tab['qseqid'])) and any(tab['sseqid'] in alignment.hit_def for alignment in xml_record.alignments)
    ]

    return pd.DataFrame(combined_results)


def count_cases(genes):
    g_genes = genes.split(".")
    cases = ["case1", "case2", "case3", "case4"]
    counts = {case: genes.count(case) for case in cases}
    rvs = genes.count("Rv")
    diff = rvs - sum(counts.values())
    return rvs, g_genes, diff, counts["case1"], counts["case2"], counts["case3"], counts["case4"]

def process_element(element):
    if "case" in element:
        length = int(element.split(";")[1].split(",")[1]) - int(element.split(";")[1].split(",")[0])
        return [element.split(";")[0], length, element]
    else:
        return [element.split("_")[0], int(element.split("_")[1].split("-")[1]) - int(element.split("_")[1].split("-")[0]), element]

def match_elements(l_bak, l_gen):
    non_matching_bak = [bak for bak in l_bak if not any(bak[1] == gen[1] for gen in l_gen)]
    non_matching_gen = [gen for gen in l_gen if not any(gen[1] == bak[1] for bak in l_bak)]
    matching_elements = [(bak, gen) for bak in l_bak for gen in l_gen if bak[1] == gen[1]]
    return non_matching_bak, non_matching_gen, matching_elements

def process_element_with_blastn(element: tuple) -> dict:
    pid = os.getpid()  # Get current process ID
    with NamedTemporaryFile(suffix=f".{pid}.fasta", delete=False) as query_file:  # Create a temporary file with the process ID
        get_seq(element[0], "nt_h37rv.fasta", query_file.name)  # Assuming get_seq is a defined function in your code
        blast_out = blastn(query_file.name, "G48915_L4database")  # Run blastn
        os.remove(query_file.name)  # Remove the temporary file

    blast_results = {}
    for blast_result in blast_out.itertuples(index=False):
        if float(blast_result.pident) > 90:  # Filter results with pident > 90
            matches = blast_result.xml_data.alignments[0].hsps[0].match.count('|') / float(blast_result.qlen) * 100
            if matches >= 80:  # Filter results with matches >= 80
                if blast_result.qseqid not in blast_results:
                    blast_results[blast_result.qseqid] = []
                blast_results[blast_result.qseqid].append((blast_result.sstart, blast_result.send, blast_result.qlen))

    return blast_results

def process_non_matching_genes(non_matching_gen: list, file_writer, bakta_coords: set) -> set:
    with ThreadPoolExecutor() as executor:
        results = list(executor.map(process_element_with_blastn, non_matching_gen))

    for result in results:
        for gene_id, coordinates_list in result.items():
            for coords in coordinates_list:
                ordered_coords = (min(int(coords[0]), int(coords[1])), max(int(coords[0]), int(coords[1])))
                coord_range = ordered_coords[0], ordered_coords[1]
                len_difference = int(coords[1]) - int(coords[0])
                length_indicator = "+" if len_difference >= 0 else "-"
                if coord_range not in bakta_coords:
                    file_writer.write(f"{gene_id}\t{ordered_coords[0]}\t{ordered_coords[1]}\t{length_indicator}\n")
                    bakta_coords.add(coord_range)
                    break  # Break once you find a match that hasn't been used yet

    return bakta_coords

def cluster1_genes(line: str, coord_list: set, file_writer) -> set:
    line_elements = line.strip("\n").split(":")
    genes = line_elements[1].split(".")
    for element in genes:
        if "case" not in element:
            gene = element.split("_")
            gene_details = [gene[0], abs(int(gene[1].split("-")[1]) - int(gene[1].split("-")[0])), element]
            coord_list = process_non_matching_genes([gene_details], file_writer, coord_list)
        elif "case1" in element:
            file_writer.write(str(element) + "\n") 
    return coord_list

def process_non_matching_elements_modified(non_matching_gen, bakta_coords):
    lista = []
    with ThreadPoolExecutor() as executor:
        results = list(executor.map(process_element_with_blastn, non_matching_gen))
    
    for result in results:        
        for qseqid, coords_list in result.items():
            for coords in coords_list:
                ordered_coords = (min(int(coords[0]), int(coords[1])), max(int(coords[0]), int(coords[1])))
                f_coord = (ordered_coords[0], ordered_coords[1])
                lengg = int(coords[1]) - int(coords[0])
                lenn = "+"
                if lengg < 0:
                    lenn = "-"
                if f_coord not in bakta_coords:
                    lista.append([qseqid, ordered_coords[0], ordered_coords[1],lenn])
                    bakta_coords.add(f_coord)
                    break  # Sal del bucle una vez que encuentres una coincidencia que no se haya utilizado todav�a
    return lista


def process_line(line, prefix):
    name = prefix.split('/')[-1] + "_"
    # Divide la linea en partes
    line_parts = line.split(".")
    # Filtra las partes para 'g_elements' y 'rv_elements'
    g_elements = (part for part in line_parts if name in part)
    rv_elements = (part for part in line_parts if "Rv" in part)

    # Analiza 'rv_elements' una vez antes del bucle usando expresiones regulares para extraer n�meros
    rv_elements_parsed = [(element, abs(int(matches[0]) - int(matches[1]))) 
                          for element in rv_elements 
                          for matches in re.findall(r'(\d+)-(\d+)', element)]

    matching_dict = {}  # Diccionario para almacenar elementos coincidentes
    used_rv_elements = set()  # Conjunto para almacenar 'rv_elements' usados
    # Procesa cada 'g_element', calcula su longitud y busca coincidencias
    for g_element in g_elements:
        el = process_element(g_element)  # Asume que process_element es una funci�n definida en tu c�digo
        g_start, g_end = map(int, el[2].split(";")[1].split(","))
        g_length = abs(g_end - g_start)

        # Empareja 'g_elements' con 'rv_elements' bas�ndose en la longitud
        for rv_element, rv_length in rv_elements_parsed:
            if rv_element not in used_rv_elements and g_length == rv_length:
                matching_dict[g_element] = rv_element
                used_rv_elements.add(rv_element)
                break  # Sal del bucle una vez que encuentres una coincidencia que no se haya utilizado todav�a
    return matching_dict

def parse_argparse():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-p', '--prefix', help='')
    args = parser.parse_args()
    return args

def main(prefix):
    
    create_database(f"{prefix}.fna", f"{prefix}database")
    
    with open(f"{prefix}.cleancoord.tsv", "w") as enter:
        with open(f"{prefix}.clusters", "r") as in_file:
            lista_coordenadas = set()
            #lista_coordenadas.add((1164572,1165123))
            #lista_coordenadas.add((3378001,3378243))

            for line in in_file:

                if "Cluster -1" in line:
                    lista_coordenadas = cluster1_genes(line, lista_coordenadas, enter)

                else:
                    l = line.strip("\n").split(":")
                    cluster = l[0]
                    rvs, g_genes, diff, case1, case2, case3, case4 = count_cases(l[1])  # print(cluster,case1,case2,case3,case4,diff,rvs)

                    elements = l[1].split(".")
                    l_bak = [process_element(element) for element in elements if "case" in element]
                    l_gen = [process_element(element) for element in elements if "case" not in element]
                    non_matching_bak, non_matching_gen, matching_elements = match_elements(l_bak, l_gen)

                    if "case4" not in line:
                        if diff == 0:
                            matching_dict = process_line(l[1], prefix)
                            for g_element, rv_element in matching_dict.items():
                                el = process_element(g_element)
                                el_split = el[2].split(";")

                                real_coord = tuple(map(int, el_split[1].split(",")))
                                if real_coord not in lista_coordenadas:
                                    lista_coordenadas.add(real_coord)
                                    enter.write(f"{el[0]}\t{real_coord[0]}\t{real_coord[1]}\t{el_split[2]}\t{rv_element}.{cluster}\n")
                        else:                            
                            for bak, gen in matching_elements:
                                if bak[1] == gen[1]:
                                    bak_split = bak[2].split(";")
                                    real_coord = tuple(map(int, bak_split[1].split(",")))
                                    if real_coord not in lista_coordenadas:
                                        lista_coordenadas.add(real_coord)
                                        enter.write(f"{bak[0]}\t{real_coord[0]}\t{real_coord[1]}\t{bak_split[2]}\t{gen[0]}.{cluster}\n")

                            lista_coordenadas = process_non_matching_genes(non_matching_gen, enter, lista_coordenadas)

                    if "case4" in line:
                        if diff < 0 or diff == 0:
                            non_matching_bak, non_matching_gen, matching_elements = match_elements(l_bak, l_gen)
                            if len(non_matching_gen) != 0:                                              
                                lista_set = set()
                                for matching_element in matching_elements:
                                    bak, gen = matching_element
                                    if bak[1] == gen[1]:                                    
                                        real_coord = tuple(map(int, bak[2].split(";")[1].split(",")))
                                        if tuple(real_coord) not in lista_coordenadas:
                                            aaa = bak[2].split(";")[2]
                                            aa = bak[2].split(";")[1].split(",")
                                            enter.write(f"{bak[0]}\t{aa[0]}\t{aa[1]}\t{aaa}\t{gen[0]}.{cluster}\n")
                                            lista_coordenadas.add(tuple(real_coord))
                                            lista_set.add(tuple(real_coord))
                                        else:
                                            if gen not in non_matching_gen:
                                                non_matching_gen.append(gen)
                                lista_res = process_non_matching_elements_modified(non_matching_gen, lista_coordenadas)  
                                no_bak = []
                                for element in lista_res:
                                    write_ele = False
                                    original_coord = element[0].split("_")[1].split("-")
                                    original_len = int(original_coord[1])-int(original_coord[0])
                                    long = int(element[2]) - int(element[1])

                                    diff_rv = long - original_len

                                    for el in non_matching_bak:
                                        diff_bk = el[1] - original_len
                                        if diff_rv == 0 or long > el[1] or (diff_bk < 450 and diff_rv < diff_bk):
                                            write_ele = True
                                            if el not in no_bak:
                                                no_bak.append(el)
                                    if write_ele:
                                        coord = (int(element[1]),int(element[2]))
                                        if tuple(coord) not in lista_set:
                                            towrite =  [str(item) for item in element]
                                            enter.write("\t".join(towrite)+"\n")
                                            lista_coordenadas.add(tuple(coord))
                                            lista_set.add(tuple(coord))
                                for bakta in non_matching_bak:
                                    if bakta not in no_bak:
                                        print(long)
                                        print(diff_rv)
                                        print("aqui1",bakta)
                                        enter.write(str(bakta) + "\n")
                            else:
                                no_bak = []
                                non_matching_gen = []
                                for sublistas in matching_elements:
                                    for elemento in sublistas:
                                        if "case" not in elemento[2]:
                                            if elemento not in non_matching_gen:
                                                non_matching_gen.append(elemento)
                                lista_set = set()
                                lista_res = process_non_matching_elements_modified(non_matching_gen, lista_coordenadas)

                                for element in lista_res:
                                    write_ele = False
                                    original_coord = element[0].split("_")[1].split("-")
                                    original_len = int(original_coord[1])-int(original_coord[0])
                                    long = int(element[2]) - int(element[1])
                                    diff_rv = long - original_len
                                    for el in non_matching_bak:
                                        diff_bk = el[1] - original_len
                                        if (long > el[1]) or (diff_bk < 450 and diff_rv < diff_bk):
                                            write_ele = True
                                            no_bak.append(el)
                                    
                                    if write_ele: 
                                        coord = (int(element[1]),int(element[2]))

                                        if tuple(coord) not in lista_set:
                                            towrite =  [str(item) for item in element]
                                            enter.write("\t".join(towrite)+"\n")

                                            lista_set.add(tuple(coord))
                                for bakta in non_matching_bak:
                                    if bakta not in no_bak:
                                        print(long)
                                        print(diff_rv)
                                        print("aqui2",bakta)
                                        enter.write(str(bakta) + "\n")

                        elif diff > 0:
                            non_matching_bak, non_matching_gen, matching_elements = match_elements(l_bak, l_gen)
                            lista_coordenadas = process_non_matching_genes(l_gen, enter, lista_coordenadas) 

if __name__ == "__main__":
    args = parse_argparse() 
    main(args.prefix)               
    #directory_path = os.path.dirname(args.prefix)
    ## All files under '/mypath/' that match '*database*'
    #file_pattern = f"{args.prefix}*database*" 
    #files = glob.glob(file_pattern) # Use glob to match the pattern
    #for f in files: # Loop over the matched files
    #    try:
    #        os.remove(f)  # Try to remove the file
    #        print(f'Successfully removed file {f}')
    #    except Exception as e:  # Catch any exceptions (like the file not existing)
    #        print(f'Failed to remove file {f}. Reason: {e}')
