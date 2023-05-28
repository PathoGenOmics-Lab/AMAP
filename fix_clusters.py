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
import os
from collections import defaultdict

def get_seq(chck_rds: list, fasta_file: str, fasta_name: str) -> List[Tuple[str, int]]:
    chck_rds_set = set(chck_rds)  # Convert the list to a set for faster searching
    selected_seqs = [seq_record for seq_record in SeqIO.parse(fasta_file, "fasta") if seq_record.id in chck_rds_set] # Use list comprehension for faster, more readable code
    SeqIO.write(selected_seqs, fasta_name, "fasta") # Write the selected sequences to a fasta file
    seq_lengths = [(seq_record.id, len(seq_record.seq)) for seq_record in selected_seqs] # Use list comprehension to create the list of sequence lengths
    return seq_lengths

def create_database(fasta: str, basename: str) -> None:
    subprocess.run(['makeblastdb', '-in', fasta, '-dbtype', 'nucl', '-out', basename])

def blastn(in_query: str, in_subject: str) -> pd.DataFrame:
    blastn_cline_tab = NcbiblastnCommandline(query=in_query, db=in_subject, outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen', dust="no")
    blastn_out_tab, _ = blastn_cline_tab()
    blastn_cline_xml = NcbiblastnCommandline(query=in_query, db=in_subject, outfmt='5', dust="no")
    blastn_out_xml, _ = blastn_cline_xml()

    tab_data = [
        dict(zip(['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs', 'qcovhsp', 'qlen', 'slen'], line.split("\t")))
        for line in blastn_out_tab.strip().split("\n")
    ] if len(blastn_out_tab.strip()) > 0 else []

    blastn_out_xml = StringIO(blastn_out_xml)
    if blastn_out_xml.getvalue().strip():
        blast_records = list(NCBIXML.parse(blastn_out_xml))
    else:
        blast_records = []

    combined_results = []

    for tab in tab_data:
        query_id = tab['qseqid']
        subject_id = tab['sseqid']
        xml_record = next((record for record in blast_records if record.query == query_id and any(subject_id in alignment.hit_def for alignment in record.alignments)), None)

        if xml_record is not None:
            combined_result = {
                "tab_data": tab,
                "xml_data": xml_record,
                "subject_name": xml_record.descriptions[0].title.split()[1] if len(xml_record.descriptions) > 0 else None
            }
            combined_results.append(combined_result)

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

    return df

def count_cases(genes):
    g_genes = genes.split(".")
    cases = ["case1", "case2", "case3", "case4"]
    counts = {case: genes.count(case) for case in cases}
    rvs = genes.count("Rv")
    diff = rvs - sum(counts.values())
    return rvs, g_genes, diff, counts["case1"], counts["case2"], counts["case3"], counts["case4"]

def process_element(element):
    if "case" in element:
        parts = element.split(";")[1].split(",")
        length = int(parts[1]) - int(parts[0])
        return [element.split(";")[0], length, element]
    else:
        parts = element.split("_")[1].split("-")
        return [element.split("_")[0], int(parts[1]) - int(parts[0]), element]

def match_elements(l_bak, l_gen):
    non_matching_bak = [bak for bak in l_bak if not any(bak[1] == gen[1] for gen in l_gen)]
    non_matching_gen = [gen for gen in l_gen if not any(gen[1] == bak[1] for bak in l_bak)]
    matching_elements = [(bak, gen) for bak in l_bak for gen in l_gen if bak[1] == gen[1]]
    return non_matching_bak, non_matching_gen, matching_elements

def process_element_processsafe2(element):
    pid = os.getpid()  # Get process ID
    with NamedTemporaryFile(suffix=f".{pid}.fasta", delete=False) as query_file:  # Use process ID in the file name
        get_seq(element[0], "nt_h37rv.fasta", query_file.name)
        blast_out = blastn(query_file.name, "G48915_L4database")
        os.remove(query_file.name)

    blast_results = {}
    for blast_result in blast_out.itertuples(index=False):
        qseqid, pident, slen, qlen, xml_data, sstart, send = blast_result.qseqid, blast_result.pident, blast_result.slen, blast_result.qlen, blast_result.xml_data, blast_result.sstart, blast_result.send

        if float(pident) > 90:
            matches = xml_data.alignments[0].hsps[0].match.count('|') / float(qlen) * 100
            if matches >= 80:
                if qseqid not in blast_results:
                    blast_results[qseqid] = []
                blast_results[qseqid].append((sstart, send, qlen))
    return blast_results

def process_non_matching_elements(non_matching_gen, enter, bakta_coords):
    with ThreadPoolExecutor() as executor:
        results = list(executor.map(process_element_processsafe2, non_matching_gen))
    
    for result in results:
        for qseqid, coords_list in result.items():
            for coords in coords_list:
                ordered_coords = (min(int(coords[0]), int(coords[1])), max(int(coords[0]), int(coords[1])))
                full_coords = (qseqid, ordered_coords[0], ordered_coords[1])
                f_coord = (ordered_coords[0], ordered_coords[1])
                lengg = int(coords[1]) - int(coords[0])
                lenn = "+"
                if lengg < 0:
                    lenn = "-"
                if f_coord not in bakta_coords:
                    enter.write(f"{qseqid}\t{ordered_coords[0]}\t{ordered_coords[1]}\t{lenn}\n")
                    bakta_coords.add(f_coord)
                    break  # Sal del bucle una vez que encuentres una coincidencia que no se haya utilizado todavía
    
    return bakta_coords

def process_non_matching_elements_modified(non_matching_gen, enter, bakta_coords):
    lista = []
    with ThreadPoolExecutor() as executor:
        results = list(executor.map(process_element_processsafe2, non_matching_gen))
    
    for result in results:        
        for qseqid, coords_list in result.items():
            for coords in coords_list:
                ordered_coords = (min(int(coords[0]), int(coords[1])), max(int(coords[0]), int(coords[1])))
                full_coords = (qseqid, ordered_coords[0], ordered_coords[1])
                f_coord = (ordered_coords[0], ordered_coords[1])
                lengg = int(coords[1]) - int(coords[0])
                lenn = "+"
                if lengg < 0:
                    lenn = "-"
                if f_coord not in bakta_coords:
                    lista.append([qseqid, ordered_coords[0], ordered_coords[1],lenn])
                    bakta_coords.add(f_coord)
                    break  # Sal del bucle una vez que encuentres una coincidencia que no se haya utilizado todavía
    return lista

def process_line(line):
    line_parts = line.split(".")
    g_elements = [part for part in line_parts if "G48915_L4_" in part]
    rv_elements = [part for part in line_parts if "Rv" in part]
    
    # Parse rv_elements once before the loop
    rv_elements_dict = {element: int(element.split("_")[1].split("-")[1]) - int(element.split("_")[1].split("-")[0]) for element in rv_elements}

    matching_dict = {}
    for g_element in g_elements:
        el = process_element(g_element)
        g_start, g_end = map(int, el[2].split(";")[1].split(","))
        g_length = g_end - g_start

        # Find a match in rv_elements_dict
        matched_rv = next((rv for rv, rv_length in rv_elements_dict.items() if g_length == rv_length), None)
        if matched_rv is not None:
            matching_dict[g_element] = matched_rv
            del rv_elements_dict[matched_rv]  # This ensures that a RV element can't be used twice

    return matching_dict


create_database("G48915_L4.fna", "G48915_L4database")

with open("AQUI.G48915_L4.cleancoord.tsv", "w") as enter:
    with open("G48915_L4.clusters", "r") as in_file:
        lista_coordenadas = set()
        lista_coordenadas.add((1164572,1165123))
        lista_coordenadas.add((3378001,3378243))

        for line in in_file:

            if "Cluster -1" in line:
                l = line.strip("\n").split(":")
                cluster = l[0]
                genes = l[1]
                g_genes = genes.split(".")
                for element in g_genes:
                    if "case" not in element:
                        gene = element
                        gen = gene.split("_")
                        lgene =[gen[0], int(gen[1].split("-")[1]) - int(gen[1].split("-")[0]), gene]
                        
                        lista_coordenadas = process_non_matching_elements([lgene], enter, lista_coordenadas)     

                    elif "case1" in element:

                        enter.write(str(element) + "\n") 

            if "Cluster -1" not in line:
                l = line.strip("\n").split(":")
                cluster = l[0]
                rvs, g_genes, diff, case1, case2, case3, case4 = count_cases(l[1])  # print(cluster,case1,case2,case3,case4,diff,rvs)

                if "case4" not in line:
                    if diff == 0:
                        matching_dict = process_line(l[1])
                        for g_element, rv_element in matching_dict.items():
                            el = process_element(g_element)
                            real_coord = tuple(map(int, el[2].split(";")[1].split(",")))
                            if real_coord not in lista_coordenadas:
                                aaa = el[2].split(";")[2]
                                aa = el[2].split(";")[1].split(",")
                                enter.write(f"{el[0]}\t{aa[0]}\t{aa[1]}\t{aaa}\t{rv_element}.{cluster}\n")
                                lista_coordenadas.add(real_coord)
                    else:
                        l_bak = [process_element(element) for element in l[1].split(".") if "case" in element]
                        l_gen = [process_element(element) for element in l[1].split(".") if "case" not in element]
                        non_matching_bak, non_matching_gen, matching_elements = match_elements(l_bak, l_gen)

                        for matching_element in matching_elements:
                            bak, gen = matching_element
                            if bak[1] == gen[1]:
                                aaa = bak[2].split(";")[2]
                                aa = bak[2].split(";")[1].split(",")
                                enter.write(f"{bak[0]}\t{aa[0]}\t{aa[1]}\t{aaa}\t{gen[0]}.{cluster}\n")
                                real_coord = tuple(map(int, bak[2].split(";")[1].split(",")))
                                lista_coordenadas.add(tuple(real_coord))
                        
                        lista_coordenadas = process_non_matching_elements(non_matching_gen, enter, lista_coordenadas)
                
                if "case4" in line:
                    if diff < 0 or diff == 0:
                        l_bak = [process_element(element) for element in l[1].split(".") if "case" in element]
                        l_gen = [process_element(element) for element in l[1].split(".") if "case" not in element]
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
                            lista_res = process_non_matching_elements_modified(non_matching_gen, enter, lista_coordenadas)  
                            no_bak = []
                            for element in lista_res:
                                write_ele = False
                                original_coord = element[0].split("_")[1].split("-")
                                original_len = int(original_coord[1])-int(original_coord[0])
                                long = int(element[2]) - int(element[1])
                                
                                diff_rv = long - original_len

                                for el in non_matching_bak:
                                    diff_bk = el[1] - original_len
                                    if diff_rv == 0:
                                        write_ele = True
                                        if el not in no_bak:
                                            no_bak.append(el)                                       
                                    if long > el[1]:
                                        write_ele = True
                                        if el not in no_bak:
                                            no_bak.append(el)  
                                    if diff_bk < 450:
                                        if diff_rv < diff_bk:
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
                            lista_res = process_non_matching_elements_modified(non_matching_gen, enter, lista_coordenadas)
                            
                            for element in lista_res:
                                write_ele = False
                                original_coord = element[0].split("_")[1].split("-")
                                original_len = int(original_coord[1])-int(original_coord[0])
                                long = int(element[2]) - int(element[1])
                                diff_rv = long - original_len

                                for el in non_matching_bak:
                                    diff_bk = el[1] - original_len
                                    if long > el[1]:
                                        write_ele = True
                                        no_bak.append(el)
                                    elif diff_bk < 450:
                                        if diff_rv < diff_bk:
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
                        l_bak = [process_element(element) for element in l[1].split(".") if "case" in element]
                        l_gen = [process_element(element) for element in l[1].split(".") if "case" not in element]
                        non_matching_bak, non_matching_gen, matching_elements = match_elements(l_bak, l_gen)
                                            
                        lista_coordenadas = process_non_matching_elements(l_gen, enter, lista_coordenadas)                    
                        
