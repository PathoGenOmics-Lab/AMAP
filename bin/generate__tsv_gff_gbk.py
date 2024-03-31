#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
import csv
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

__author__ = 'Paula Ruiz-Rodriguez'
__credits__ = ['Paula Ruiz-Rodriguez', 'PathoGenOmics']
__license__ = ['GPL 3']
__version__ = '0.1.0'
__maintainer__ = '@paururo'
__email__ = 'paula.ruiz-rodriguez@uv.es'
__status__ = 'developing'
__lab__ = 'PathoGenOmics (I2SysBio)'

def read_pandas(open_file, name_columns):
    df = pd.read_csv(open_file, comment='#',names=name_columns, sep='\t')
    return df

def filter_bakta(bakta_pd):
    bakta_pd_cds = bakta_pd[bakta_pd['Type'] == 'cds']
    bakta_pd_rest = bakta_pd[bakta_pd['Type'] != 'cds']
    return bakta_pd_cds, bakta_pd_rest

def tsv_to_gbk(tsv_file, fasta_file, gbk_file, molecule_type='DNA'):
    secuencias = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    with open(tsv_file, 'r') as tsv:
        tsv_reader = csv.reader(tsv, delimiter='\t')
        print(tsv_reader)
        next(tsv_reader, None) 
        
        for record in secuencias.values():
            record.annotations['molecule_type'] = molecule_type
            for row in tsv_reader:
                if row[0].startswith('#'):
                    continue
                seq_id, feature_type, start, end, strand, locus_tag, gene, product, db_xrefs = row
                if seq_id not in secuencias:
                    continue 
                start, end = int(start) - 1, int(end)
                strand = 1 if strand == '+' else -1 if strand == '-' else 0
                qualifiers = {
                    'locus_tag': [locus_tag],
                    'gene': [gene],
                    'product': [product],
                    'db_xref': db_xrefs.split(', ')
                }
                if feature_type.lower() == 'cds':
                    sequence = secuencias[seq_id].seq[start:end]
                    if strand == -1:
                        sequence = sequence.reverse_complement()
                    qualifiers['translation'] = [str(sequence.translate(to_stop=True))]
                feature = SeqFeature(FeatureLocation(start, end, strand=strand), type=feature_type, qualifiers=qualifiers)
                secuencias[seq_id].features.append(feature)
            tsv.seek(0) 
        with open(gbk_file, 'w') as gbk:
            SeqIO.write(secuencias.values(), gbk, 'genbank')

def columns():
    name_eggnog = [
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
    "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
    "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
    "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"
    ]
    name_correctedtsv = [
    "id", "strand", "start", "end", "rvs"
    ]
    name_bakta = [
        "Sequence Id","Type", "Start","Stop","Strand","Locus Tag", "Gene","Product","DbXrefs"
    ]
    name_myco = [
        "Refseq_ID", "Mycobrowser_Version", "Feature", "Start", "Stop", "Score", 
        "Strand", "Frame", "Locus", "Name", "Function", "Product", "Comments", "UniProt_AC", 
        "Is_Pseudogene", "Functional_Category", "Protein Data Bank", "PFAM", "UniProt", 
        "Gene Ontology", "Enzyme Classification", "Drug Resistance Mutations", "InterPro", 
        "UniProtKB/TrEMBL", "SWISS-MODEL", "Orthologues M. leprae", "Orthologues M. marinum", 
        "Orthologues M. smegmatis", "Orthologues M. bovis", "Orthologues M. lepromatosis", 
        "Orthologues M. tuberculosis", "Orthologues M. abscessus", "Orthologues M. haemophilum", 
        "Orthologues M. orygis"
    ]
    return name_eggnog, name_correctedtsv, name_bakta, name_myco

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("egg_nog", help="EggnogTsvFile")
    parser.add_argument("correctedtsv", help="CorrectedTsvFile")
    parser.add_argument("bakta", help="BaktaTsvFile")
    parser.add_argument("myco", help="MycobrowserTsvFile")
    parser.add_argument("fasta", help="FastaFile")
    args = parser.parse_args()
    base_name = args.correctedtsv.replace('_corrected.tsv','')

    name_eggnog, name_correctedtsv, name_bakta, name_myco = columns()

    eggnog_pd = read_pandas(args.egg_nog, name_eggnog)
    eggnog_pd['simple_id'] = eggnog_pd['query'].str.split(';').str[0]
    correctedtsv_pd = read_pandas(args.correctedtsv, name_correctedtsv)
    correctedtsv_pd['simple_id'] = correctedtsv_pd['id'].str.split(';').str[0]
    baktatsv_pd = read_pandas(args.bakta, name_bakta)
    bakta_pd_cds, bakta_pd_rest = filter_bakta(baktatsv_pd)
    myco_pd = read_pandas(args.myco, name_myco)

    merged_df = pd.merge(correctedtsv_pd, eggnog_pd, on='simple_id', how='left')
    merged_df = pd.merge(merged_df, bakta_pd_cds, left_on='simple_id', right_on='Locus Tag', how='left')
    merged_df = merged_df.drop(columns=['query'])
    merged_df.drop_duplicates(subset='simple_id', inplace=True)
    merged_df.reset_index(drop=True, inplace=True)

    merged_df['Locus'] = correctedtsv_pd['rvs'].str.split('_').str[0]

    merged_df = pd.merge(merged_df, myco_pd, on='Locus', how='left')
    columns_to_drop = ['Start_x', 'Stop_x', 'Start_y', 'Stop_y', 'Strand_x', 'Strand_y', 
                       'Product_y', 'Sequence Id', 'Type', 'Locus Tag']
    merged_df = merged_df.drop(columns=columns_to_drop)
    merged_df.drop_duplicates(subset='simple_id', inplace=True)
    merged_df.to_csv(base_name + "-mastertable.tsv", sep='\t', index=False)

    nuevo_df = merged_df[['id','strand','start','end','rvs','Preferred_name', 'Gene', 'Product_x','Description']]
    nuevo_df['Preferred_name'] = nuevo_df['Preferred_name'].replace('-', np.nan)
    nuevo_df['Description'] = nuevo_df['Description'].replace('-', np.nan)
    nuevo_df['Combined_Gene'] = nuevo_df['Preferred_name'] + ';' + nuevo_df['Gene']
    nuevo_df['Combined_Description'] = nuevo_df['Product_x'] + ';' + nuevo_df['Description'] 
    nuevo_df = nuevo_df.drop(columns=['Gene','Preferred_name', 'Product_x', 'Description'])
    nuevo_df['Combined_Gene'] = nuevo_df['Combined_Gene'].str.split(';').str[1]
    nuevo_df['Combined_Description'] = nuevo_df['Combined_Description'].str.split(';').str[1]
    nuevo_df['rvs'] = 'AMAP_ids:' + nuevo_df['rvs'].astype(str)
    nuevos_nombres = {'id': 'Locus Tag', 'strand': 'Strand', 'start': 'Start', 'end':'Stop', 'strand':'Strand', 'Combined_Gene':'Gene', 'Combined_Description':'Product', 'rvs':'DbXrefs'}  # Ajusta según tu necesidad

    nuevo_df.rename(columns=nuevos_nombres, inplace=True)
    nuevo_df['Sequence Id'] = "contig_1"
    nuevo_df['Type'] = "CDS"
    orden = ["Sequence Id", "Type", "Start", "Stop", "Strand", "Locus Tag", "Gene", "Product", "DbXrefs"]
    nuevo_df = nuevo_df[orden]
    nuevo_df = pd.concat([bakta_pd_rest, nuevo_df], ignore_index=True)
    nuevo_df = nuevo_df.sort_values(by='Start')
    
    for secuencia in SeqIO.parse(args.fasta, "fasta"):
        nuevo_df['Sequence Id'] = secuencia.id
    
    lineas_adicionales = [
    "# Annotated with AMAP",
    "# Software: v0.1.0",
    "# URL: https://github.com/PathoGenOmics-Lab/AMAP",
    "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs"
    ]

    cabecera = '\n'.join(lineas_adicionales) + '\n'

    ruta_archivo = str(base_name+"-AMAP.tsv")
    with open(ruta_archivo, 'w') as archivo:
        archivo.write(cabecera)

        nuevo_df.to_csv(archivo, sep='\t', index=False, header=False)
    
    nuevo_df['Attributes'] = nuevo_df.apply(lambda row: f"ID={row['Locus Tag']};Name={row['Gene']};product={row['Product']};Dbxref={row['DbXrefs']}", axis=1)
    nuevo_df = nuevo_df[['Sequence Id', 'Type', 'Start', 'Stop', 'Strand', 'Attributes']]

    nuevo_df['Source'] = 'AMAP'
    nuevo_df['Score'] = '.'
    nuevo_df['Phase'] = '.'

    ruta_archivo_gff = str(base_name+"-AMAP.gff3")
    nuevo_df = nuevo_df[['Sequence Id', 'Source', 'Type', 'Start', 'Stop', 'Score', 'Strand', 'Phase', 'Attributes']]
    nuevo_df.to_csv(ruta_archivo_gff, sep='\t', index=False, header=False)
    with open(ruta_archivo_gff, 'a') as archivo_gff:
        archivo_gff.write("\n##FASTA\n")

        for secuencia in SeqIO.parse(args.fasta, "fasta"):
            archivo_gff.write(f">{secuencia.id}\n{secuencia.seq}\n")

    ruta_archivo_fasta = args.fasta

    tsv_to_gbk(ruta_archivo, ruta_archivo_fasta, base_name + '-AMAP.gbk')

if __name__ == "__main__":
    main()
#python3 generate_gff.py test3.emapper.annotations _ancestor_corrected.tsv _ancestor.tsv  Mycobacterium_tuberculosis_H37Rv_txt_v4.txt /home/ruizro/test_nextflow/mis_fastas/MTB_ancestor_reference.fasta