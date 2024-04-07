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

def generate_master_table(eggnog_pd, correctedtsv_pd, bakta_pd_cds, myco_pd, base_name):
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
    merged_df.to_csv(base_name + "-AMAPmastertable.tsv", sep='\t', index=False)
    return merged_df


def generate_tsv(merged_df, bakta_rest_df, fasta_file_path, output_file_path):
    # Create a copy of the merged DataFrame to avoid SettingWithCopyWarning
    new_df = merged_df[['id', 'strand', 'start', 'end', 'rvs', 'Preferred_name', 'Gene', 'Product_x', 'Description']].copy()

    # Use loc to set values and avoid warnings
    new_df.loc[:, 'Preferred_name'] = new_df['Preferred_name'].replace('-', np.nan)
    new_df.loc[:, 'Description'] = new_df['Description'].replace('-', np.nan)
    new_df.loc[:, 'Combined_Gene'] = new_df['Preferred_name'] + ';' + new_df['Gene']
    new_df.loc[:, 'Combined_Description'] = new_df['Product_x'] + ';' + new_df['Description']
    
    # Dropping columns
    new_df.drop(columns=['Gene', 'Preferred_name', 'Product_x', 'Description'], inplace=True)

    # Splitting and renaming columns
    new_df.loc[:, 'Combined_Gene'] = new_df['Combined_Gene'].str.split(';').str[1]
    new_df.loc[:, 'Combined_Description'] = new_df['Combined_Description'].str.split(';').str[1]
    new_df.loc[:, 'rvs'] = 'AMAP_ids:' + new_df['rvs'].astype(str)
    
    # Renaming columns
    new_column_names = {'id': 'Locus Tag', 'strand': 'Strand', 'start': 'Start', 'end': 'Stop', 'Combined_Gene': 'Gene', 'Combined_Description': 'Product', 'rvs': 'DbXrefs'}
    new_df.rename(columns=new_column_names, inplace=True)

    # Setting additional columns
    new_df['Sequence Id'] = "contig_1"
    new_df['Type'] = "CDS"
    
    # Reordering columns
    column_order = ["Sequence Id", "Type", "Start", "Stop", "Strand", "Locus Tag", "Gene", "Product", "DbXrefs"]
    new_df = new_df[column_order]

    # Concatenating with bakta_rest_df and sorting
    new_df = pd.concat([bakta_rest_df, new_df], ignore_index=True).sort_values(by='Start')

    # Update 'Sequence Id' for each sequence in FASTA file (this might need adjustment based on specific needs)
    for sequence in SeqIO.parse(fasta_file_path, "fasta"):
        new_df.loc[new_df['Sequence Id'] == "contig_1", 'Sequence Id'] = sequence.id

    # Header for TSV file
    header_content = "\n".join([
        "# Annotated with AMAP",
        "# Software: v0.1.0",
        "# URL: https://github.com/PathoGenOmics-Lab/AMAP",
        "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs"
    ]) + "\n"
    
    # Writing to file
    with open(output_file_path, 'w') as file:
        file.write(header_content)
        new_df.to_csv(file, sep='\t', index=False, header=False)

    return new_df

def generate_gff(new_df, gff_file_path, fasta_file_path):
    new_df = new_df.copy()
    # Use apply to create the 'Attributes' column
    new_df['Attributes'] = new_df.apply(lambda row: f"ID={row['Locus Tag']};Name={row['Gene']};product={row['Product']};Dbxref={row['DbXrefs']}", axis=1)
    new_df = new_df[['Sequence Id', 'Type', 'Start', 'Stop', 'Strand', 'Attributes']]
    # Use loc to set new values
    new_df.loc[:, 'Source'] = 'AMAP'
    new_df.loc[:, 'Score'] = '.'
    new_df.loc[:, 'Phase'] = '.'

    # Rearrange columns
    new_df = new_df[['Sequence Id', 'Source', 'Type', 'Start', 'Stop', 'Score', 'Strand', 'Phase', 'Attributes']]
   
    # Writing the DataFrame to a GFF file
    new_df.to_csv(gff_file_path, sep='\t', index=False, header=False)

    # Appending FASTA sequences to the GFF file
    with open(gff_file_path, 'a') as file_gff:
        file_gff.write("\n##FASTA\n")

        for sequence in SeqIO.parse(fasta_file_path, "fasta"):
            file_gff.write(f">{sequence.id}\n{sequence.seq}\n")

    return new_df

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("egg_nog", help="EggnogTsvFile")
    parser.add_argument("correctedtsv", help="CorrectedTsvFile")
    parser.add_argument("bakta", help="BaktaTsvFile")
    parser.add_argument("myco", help="MycobrowserTsvFile")
    parser.add_argument("fasta", help="FastaFile")
    args = parser.parse_args()
    base_name = args.correctedtsv.replace('.AMAP.tsv','')
    ruta_archivo_fasta = args.fasta
    ruta_archivo = str(base_name+"-AMAP.tsv")
    ruta_archivo_gff = str(base_name+"-AMAP.gff3")
    name_eggnog, name_correctedtsv, name_bakta, name_myco = columns()

    eggnog_pd = read_pandas(args.egg_nog, name_eggnog)
    eggnog_pd['simple_id'] = eggnog_pd['query'].str.split(';').str[0]
    correctedtsv_pd = read_pandas(args.correctedtsv, name_correctedtsv)
    correctedtsv_pd['simple_id'] = correctedtsv_pd['id'].str.split(';').str[0]
    baktatsv_pd = read_pandas(args.bakta, name_bakta)
    bakta_pd_cds, bakta_pd_rest = filter_bakta(baktatsv_pd)
    myco_pd = read_pandas(args.myco, name_myco)

    merged_df = generate_master_table(eggnog_pd, correctedtsv_pd, bakta_pd_cds, myco_pd, base_name)
    nuevo_df = generate_tsv(merged_df,bakta_pd_rest,ruta_archivo_fasta, ruta_archivo)
    generate_gff(nuevo_df, ruta_archivo_gff,ruta_archivo_fasta)
    tsv_to_gbk(ruta_archivo, ruta_archivo_fasta, base_name + '-AMAP.gbk')

if __name__ == "__main__":
    main()
    #python3 generate_gff.py test3.emapper.annotations _ancestor_corrected.tsv _ancestor.tsv  Mycobacterium_tuberculosis_H37Rv_txt_v4.txt /home/ruizro/test_nextflow/mis_fastas/MTB_ancestor_reference.fasta