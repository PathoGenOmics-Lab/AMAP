#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import csv
from Bio.SeqFeature import SeqFeature, FeatureLocation

__author__ = 'Paula Ruiz-Rodriguez'
__credits__ = ['Paula Ruiz-Rodriguez', 'PathoGenOmics']
__license__ = ['GPL 3']
__version__ = '0.1.0'
__maintainer__ = '@paururo'
__email__ = 'paula.ruiz.rodriguez@csic.es'
__status__ = 'finished'
__lab__ = 'PathoGenOmics (I2SysBio)'

def read_pandas(open_file, name_columns):
    df = pd.read_csv(open_file, comment='#',names=name_columns, sep='\t')
    return df

def filter_bakta(bakta_pd):
    bakta_pd_cds = bakta_pd[bakta_pd['Type'] == 'cds']
    bakta_pd_rest = bakta_pd[bakta_pd['Type'] != 'cds']
    return bakta_pd_cds, bakta_pd_rest

def tsv_to_gbk(tsv_file, fasta_file, gbk_file, molecule_type='DNA'):
    """Convert a TSV annotation file + FASTA to GenBank format."""
    secuencias = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    
    for record in secuencias.values():
        record.annotations['molecule_type'] = molecule_type

    # Read all features from the TSV file
    with open(tsv_file, 'r') as tsv:
        tsv_reader = csv.reader(tsv, delimiter='\t')
        for row in tsv_reader:
            # Skip header and comment lines
            if not row or row[0].startswith('#'):
                continue
            if len(row) < 9:
                continue
            seq_id, feature_type, start, end, strand, locus_tag, gene, product, db_xrefs = row[:9]
            if seq_id not in secuencias:
                continue
            try:
                start, end = int(start) - 1, int(end)
            except ValueError:
                continue
            strand_int = 1 if strand == '+' else -1 if strand == '-' else 0
            qualifiers = {
                'locus_tag': [locus_tag],
                'gene': [gene],
                'product': [product],
                'db_xref': db_xrefs.split(', ')
            }
            if feature_type.lower() == 'cds':
                sequence = secuencias[seq_id].seq[start:end]
                if strand_int == -1:
                    sequence = sequence.reverse_complement()
                qualifiers['translation'] = [str(sequence.translate(to_stop=True))]
            feature = SeqFeature(FeatureLocation(start, end, strand=strand_int), type=feature_type, qualifiers=qualifiers)
            secuencias[seq_id].features.append(feature)

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
    merged_df['Locus'] = merged_df['rvs'].str.split('_').str[0]
    merged_df = pd.merge(merged_df, myco_pd, on='Locus', how='left')
    columns_to_drop = ['Start_x', 'Stop_x', 'Start_y', 'Stop_y', 'Strand_x', 'Strand_y', 
                       'Product_y', 'Sequence Id', 'Type', 'Locus Tag']
    # Only drop columns that actually exist after the merges
    existing_cols_to_drop = [c for c in columns_to_drop if c in merged_df.columns]
    merged_df = merged_df.drop(columns=existing_cols_to_drop)
    merged_df.drop_duplicates(subset='simple_id', inplace=True)
    merged_df.to_csv(base_name + "-AMAPmastertable.tsv", sep='\t', index=False)
    return merged_df


def generate_tsv(merged_df, bakta_rest_df, fasta_file_path, output_file_path):
    # Create a copy of the merged DataFrame to avoid SettingWithCopyWarning
    new_df = merged_df[['id', 'strand', 'start', 'end', 'rvs', 'Preferred_name', 'Gene', 'Product_x', 'Description']].copy()

    # Use loc to set values and avoid warnings
    new_df.loc[:, 'Preferred_name'] = new_df['Preferred_name'].replace('-', np.nan)
    new_df.loc[:, 'Description'] = new_df['Description'].replace('-', np.nan)

    # Combine gene names: prefer Bakta Gene over eggNOG Preferred_name, fallback chain
    new_df.loc[:, 'Combined_Gene'] = new_df['Gene'].fillna(new_df['Preferred_name'])
    # Combine descriptions: prefer eggNOG Description over Bakta Product_x, fallback chain
    new_df.loc[:, 'Combined_Description'] = new_df['Description'].fillna(new_df['Product_x'])

    # Dropping columns
    new_df.drop(columns=['Gene', 'Preferred_name', 'Product_x', 'Description'], inplace=True)
    new_df.loc[:, 'rvs'] = 'AMAP_ids:' + new_df['rvs'].fillna('').astype(str)
    
    # Renaming columns
    new_column_names = {'id': 'Locus Tag', 'strand': 'Strand', 'start': 'Start', 'end': 'Stop', 'Combined_Gene': 'Gene', 'Combined_Description': 'Product', 'rvs': 'DbXrefs'}
    new_df.rename(columns=new_column_names, inplace=True)

    # Setting additional columns — assign Sequence Id per gene
    # For multi-contig assemblies, match each gene to the correct contig by coordinates
    fasta_records = list(SeqIO.parse(fasta_file_path, "fasta"))
    if len(fasta_records) == 1:
        new_df['Sequence Id'] = fasta_records[0].id
    else:
        # Sort contigs by length desc, assign each gene to the smallest fitting contig
        sorted_contigs = sorted(fasta_records, key=lambda r: len(r.seq), reverse=True)
        def find_contig_id(row):
            end = row['Stop'] if not pd.isna(row['Stop']) else 0
            best = sorted_contigs[0].id
            for contig in sorted_contigs:
                if int(end) <= len(contig.seq):
                    best = contig.id
            return best
        new_df['Sequence Id'] = new_df.apply(find_contig_id, axis=1)
    new_df['Type'] = "CDS"
    
    # Reordering columns
    column_order = ["Sequence Id", "Type", "Start", "Stop", "Strand", "Locus Tag", "Gene", "Product", "DbXrefs"]
    new_df = new_df[column_order]

    # Concatenating with bakta_rest_df and sorting
    new_df = pd.concat([bakta_rest_df, new_df], ignore_index=True).sort_values(by='Start')

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
    new_df.loc[:, 'Phase'] = new_df['Type'].apply(lambda t: '0' if t == 'CDS' else '.')

    # Rearrange columns
    new_df = new_df[['Sequence Id', 'Source', 'Type', 'Start', 'Stop', 'Score', 'Strand', 'Phase', 'Attributes']]
   
    # Writing the DataFrame to a GFF3 file with proper header
    with open(gff_file_path, 'w') as file_gff:
        file_gff.write("##gff-version 3\n")
        new_df.to_csv(file_gff, sep='\t', index=False, header=False)

        # Append FASTA sequences
        file_gff.write("##FASTA\n")
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

    # Mycobrowser is optional — skip if "NONE" or file doesn't exist
    use_mycobrowser = args.myco != "NONE" and os.path.isfile(args.myco)
    if use_mycobrowser:
        myco_pd = read_pandas(args.myco, name_myco)
    else:
        myco_pd = pd.DataFrame(columns=['Locus'] + name_myco)

    merged_df = generate_master_table(eggnog_pd, correctedtsv_pd, bakta_pd_cds, myco_pd, base_name)
    nuevo_df = generate_tsv(merged_df,bakta_pd_rest,ruta_archivo_fasta, ruta_archivo)
    generate_gff(nuevo_df, ruta_archivo_gff,ruta_archivo_fasta)
    tsv_to_gbk(ruta_archivo, ruta_archivo_fasta, base_name + '-AMAP.gbk')

if __name__ == "__main__":
    main()
    #python3 generate_gff.py test3.emapper.annotations _ancestor_corrected.tsv _ancestor.tsv  Mycobacterium_tuberculosis_H37Rv_txt_v4.txt /home/ruizro/test_nextflow/mis_fastas/MTB_ancestor_reference.fasta
