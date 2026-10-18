#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
import csv
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Allow importing from bin/classify/ when run as a script, the same way
# classify_cases.py does. classify.ids only needs `re`, so this stays cheap.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classify.ids import clean_id

__author__ = 'Paula Ruiz-Rodriguez'
__credits__ = ['Paula Ruiz-Rodriguez', 'PathoGenOmics']
__license__ = ['GPL 3']
__version__ = '0.1.0'
__maintainer__ = '@paururo'
__email__ = 'paula.ruiz.rodriguez@csic.es'
__status__ = 'finished'
__lab__ = 'PathoGenOmics (I2SysBio)'

def read_pandas(open_file, name_columns):
    df = pd.read_csv(open_file, comment='#', names=name_columns, sep='\t')
    # `names=` disables header inference, so a file that DOES carry a header
    # row yields it as a data row (a phantom gene named "id", and Start/Stop
    # parsed as text). GET_OUTPUTS strips the header before calling this, but
    # the script is also run standalone, so drop it here when it is present.
    # A table with fewer columns than names (one from before the contig column
    # was added) leaves the last ones empty, so match the header on its prefix.
    first = [str(v) for v in df.iloc[0] if pd.notna(v)] if len(df) > 0 else []
    if first and first == list(name_columns)[:len(first)]:
        df = df.iloc[1:].reset_index(drop=True)
        # The stray header row forced every column to object dtype; restore
        # numeric typing so coordinates sort numerically instead of as text.
        for col in df.columns:
            converted = pd.to_numeric(df[col], errors='coerce')
            if converted.notna().equals(df[col].notna()):
                df[col] = converted
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
                # Use full translation (not to_stop) to handle pseudogenes
                # with internal stop codons. Strip trailing * if present.
                full_translation = str(sequence.translate()).rstrip('*')
                # If there are internal stops, mark as pseudogene in qualifiers
                if '*' in full_translation:
                    qualifiers['pseudo'] = ['']
                    qualifiers['note'] = [f'contains {full_translation.count("*")} internal stop codon(s)']
                    # GenBank convention: pseudogenes use to_stop translation
                    qualifiers['translation'] = [str(sequence.translate(to_stop=True))]
                else:
                    qualifiers['translation'] = [full_translation]
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
    # check_cds.py writes the contig last; a table from before contigs were
    # carried through reads with it empty.
    name_correctedtsv = [
    "id", "strand", "start", "end", "rvs", "pseudogene", "contig"
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

def locus_from_rvs(value):
    """Build the Mycobrowser join key from an entry of the `rvs` column.

    `rvs` holds reference protein IDs of the form GENE_START-END, comma
    separated when a gene maps to more than one, e.g.
    'Rv0001_1-1524' or 'Rv0001_1-1524, Rv0002_2052-3260'.

    This used to be rvs.split('_')[0], which truncates any gene whose name
    itself contains an underscore: PE_PGRS1_703830-704057 became 'PE', so every
    PE/PPE gene of a sample collapsed onto the same key and joined against
    whatever Mycobrowser had under 'PE'. Prefixed identifiers in
    organism-agnostic mode broke the same way (ECOLI_thrA_337-2799 -> 'ECOLI').

    Take the first ID and strip only the trailing coordinates, which is what
    clean_id already does for the rest of the pipeline.
    """
    if pd.isna(value):
        return value
    first = str(value).split(',')[0].strip()
    return clean_id(first)


def generate_master_table(eggnog_pd, correctedtsv_pd, bakta_pd_cds, myco_pd, base_name):
    merged_df = pd.merge(correctedtsv_pd, eggnog_pd, on='simple_id', how='left')
    merged_df = pd.merge(merged_df, bakta_pd_cds, left_on='simple_id', right_on='Locus Tag', how='left')
    merged_df = merged_df.drop(columns=['query'])
    merged_df.drop_duplicates(subset='simple_id', inplace=True)

    # Add explicit lineage column for downstream R report.
    # base_name is "L4_G26689" (typed) or "G26689" (untyped).
    # Without this column, the R report extracts lineage from the first column
    # using sub("_.*", ""), which fails for untyped samples (returns genome ID).
    parts = base_name.split("_")
    merged_df['lineage'] = parts[0] if len(parts) >= 2 else "Untyped"
    merged_df.reset_index(drop=True, inplace=True)
    merged_df['Locus'] = merged_df['rvs'].apply(locus_from_rvs)
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
    keep = ['id', 'strand', 'start', 'end', 'rvs', 'Preferred_name', 'Gene', 'Product_x', 'Description']
    new_df = merged_df[keep + (['contig'] if 'contig' in merged_df.columns else [])].copy()

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

    # Setting additional columns: assign Sequence Id per gene, from the contig
    # check_cds.py recorded. Guessing it from the coordinates, the smallest
    # contig they fit in, put 90% of the genes of a 40-contig draft on the
    # wrong contig; it remains only for tables that carry no contig.
    fasta_records = list(SeqIO.parse(fasta_file_path, "fasta"))
    fasta_ids = {r.id for r in fasta_records}
    if len(fasta_records) == 1:
        new_df['Sequence Id'] = fasta_records[0].id
    elif 'contig' in new_df.columns and new_df['contig'].isin(fasta_ids).all():
        new_df['Sequence Id'] = new_df['contig']
    else:
        print("WARNING: no contig recorded for the genes of a multi-contig assembly; "
              "placing each on the smallest contig its coordinates fit in", file=sys.stderr)
        # Assign each gene to the smallest contig where it fits (tightest match)
        sorted_contigs = sorted(fasta_records, key=lambda r: len(r.seq), reverse=True)
        def find_contig_id(row):
            start = int(row['Start']) if not pd.isna(row['Start']) else 0
            end = int(row['Stop']) if not pd.isna(row['Stop']) else 0
            best = sorted_contigs[0].id
            best_len = float('inf')
            for contig in sorted_contigs:
                clen = len(contig.seq)
                if start >= 1 and end <= clen and clen < best_len:
                    best = contig.id
                    best_len = clen
            return best
        new_df['Sequence Id'] = new_df.apply(find_contig_id, axis=1)
    new_df['Type'] = "CDS"
    
    # Bakta renames every sequence to contig_N unless run with
    # --keep-contig-headers, and then its tRNA, rRNA and other features sit on
    # a sequence the FASTA does not have: tsv_to_gbk dropped them and the GFF3
    # listed them on an undeclared contig. One sequence can only be that one.
    bakta_rest_df = bakta_rest_df.copy()
    if len(fasta_records) == 1:
        bakta_rest_df['Sequence Id'] = fasta_records[0].id
    elif len(bakta_rest_df) and not bakta_rest_df['Sequence Id'].isin(fasta_ids).all():
        print("WARNING: Bakta's sequence names do not match the FASTA; run it with "
              "--keep-contig-headers", file=sys.stderr)

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

def _escape_gff3(value):
    """Escape special characters in GFF3 attribute values (spec §2.1)."""
    if pd.isna(value) or value is None:
        return ''
    value = str(value)
    # GFF3 requires URL-encoding of: tab, newline, carriage return, %, ;, =, &, ,
    for char, enc in [('%', '%25'), ('\t', '%09'), ('\n', '%0A'), ('\r', '%0D'),
                      (';', '%3B'), ('=', '%3D'), ('&', '%26'), (',', '%2C')]:
        value = value.replace(char, enc)
    return value


def generate_gff(merged_df, tsv_df, gff_file_path, fasta_file_path):
    """Generate a GFF3 file with proper spec compliance and eggNOG annotation.
    
    Follows GFF3 specification (https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md):
    - ##gff-version 3 header
    - ##sequence-region pragmas for each contig
    - Gene features as parents of CDS features
    - eggNOG functional attributes: Ontology_term, eC_number, COG, KEGG, PFAMs
    - Proper attribute escaping (special chars URL-encoded)
    - ##FASTA section with embedded sequences
    - Sorted by coordinate within each contig
    """
    tsv_df = tsv_df.copy()
    
    # Build a lookup from locus tag to eggNOG annotations
    eggnog_lookup = {}
    if merged_df is not None:
        for _, row in merged_df.iterrows():
            gene_id = str(row.get('id', '')) if not pd.isna(row.get('id', np.nan)) else ''
            if not gene_id:
                continue
            info = {}
            for col in ['GOs', 'EC', 'COG_category', 'Description', 'KEGG_ko',
                        'KEGG_Pathway', 'PFAMs']:
                val = row.get(col, '')
                if pd.isna(val) or str(val).strip() in ('', '-'):
                    continue
                info[col] = str(val).strip()
            if info:
                eggnog_lookup[gene_id] = info

    # Sort by contig and start position for spec compliance
    tsv_df['Start_int'] = pd.to_numeric(tsv_df['Start'], errors='coerce').fillna(0).astype(int)
    tsv_df = tsv_df.sort_values(['Sequence Id', 'Start_int']).drop(columns=['Start_int'])

    # Parse FASTA for sequence-region pragmas
    fasta_records = list(SeqIO.parse(fasta_file_path, "fasta"))
    contig_lengths = {rec.id: len(rec.seq) for rec in fasta_records}

    with open(gff_file_path, 'w') as file_gff:
        file_gff.write("##gff-version 3\n")
        
        # Write ##sequence-region pragmas (GFF3 spec §4.1)
        for contig_id in sorted(contig_lengths.keys()):
            file_gff.write(f"##sequence-region {contig_id} 1 {contig_lengths[contig_id]}\n")

        for _, row in tsv_df.iterrows():
            seq_id = row['Sequence Id']
            feat_type = row['Type']
            start = row['Start']
            stop = row['Stop']
            strand = row['Strand']
            locus_tag = _escape_gff3(row['Locus Tag'])
            gene = _escape_gff3(row['Gene'])
            product = _escape_gff3(row['Product'])
            dbxrefs = _escape_gff3(row['DbXrefs'])
            phase = '0' if feat_type == 'CDS' else '.'

            if feat_type == 'CDS':
                # Emit gene feature first (GFF3 spec compliance)
                gene_id = f"gene-{locus_tag}"
                gene_attrs = f"ID={gene_id};Name={gene};locus_tag={locus_tag}"
                file_gff.write(f"{seq_id}\tAMAP\tgene\t{start}\t{stop}\t.\t{strand}\t.\t{gene_attrs}\n")

                # Build CDS attributes with eggNOG functional annotation
                cds_attrs = f"ID={locus_tag};Parent={gene_id};Name={gene};product={product};Dbxref={dbxrefs}"

                # Add eggNOG attributes if available
                einfo = eggnog_lookup.get(row['Locus Tag'], {})
                if 'GOs' in einfo:
                    cds_attrs += f";Ontology_term={_escape_gff3(einfo['GOs'])}"
                if 'EC' in einfo:
                    cds_attrs += f";eC_number={_escape_gff3(einfo['EC'])}"
                if 'COG_category' in einfo:
                    cds_attrs += f";COG_category={_escape_gff3(einfo['COG_category'])}"
                if 'Description' in einfo:
                    cds_attrs += f";Note={_escape_gff3(einfo['Description'])}"
                if 'KEGG_ko' in einfo:
                    cds_attrs += f";KEGG_ko={_escape_gff3(einfo['KEGG_ko'])}"
                if 'KEGG_Pathway' in einfo:
                    cds_attrs += f";KEGG_Pathway={_escape_gff3(einfo['KEGG_Pathway'])}"
                if 'PFAMs' in einfo:
                    cds_attrs += f";PFAMs={_escape_gff3(einfo['PFAMs'])}"

                file_gff.write(f"{seq_id}\tAMAP\t{feat_type}\t{start}\t{stop}\t.\t{strand}\t{phase}\t{cds_attrs}\n")
            else:
                # Non-CDS features (tRNA, rRNA, etc.)
                attrs = f"ID={locus_tag};Name={gene};product={product};Dbxref={dbxrefs}"
                file_gff.write(f"{seq_id}\tAMAP\t{feat_type}\t{start}\t{stop}\t.\t{strand}\t{phase}\t{attrs}\n")

        # Append FASTA sequences (GFF3 spec §4.2)
        file_gff.write("##FASTA\n")
        for sequence in fasta_records:
            file_gff.write(f">{sequence.id}\n{sequence.seq}\n")

    return tsv_df

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
        # name_myco already contains 'Locus' (the merge key). Prepending it
        # again produced a duplicate column label, and pd.merge(on='Locus')
        # then raised "The column label 'Locus' is not unique", killing every
        # run without Mycobrowser (--skip_mycobrowser / organism-agnostic mode).
        myco_pd = pd.DataFrame(columns=name_myco)

    merged_df = generate_master_table(eggnog_pd, correctedtsv_pd, bakta_pd_cds, myco_pd, base_name)
    nuevo_df = generate_tsv(merged_df,bakta_pd_rest,ruta_archivo_fasta, ruta_archivo)
    generate_gff(merged_df, nuevo_df, ruta_archivo_gff,ruta_archivo_fasta)
    tsv_to_gbk(ruta_archivo, ruta_archivo_fasta, base_name + '-AMAP.gbk')

if __name__ == "__main__":
    main()
    #python3 generate_gff.py test3.emapper.annotations _ancestor_corrected.tsv _ancestor.tsv  Mycobacterium_tuberculosis_H37Rv_txt_v4.txt /home/ruizro/test_nextflow/mis_fastas/MTB_ancestor_reference.fasta
