#! /usr/bin/python3
# -*- coding: utf-8 -*-

import csv
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import shutil

def tsv_to_gbk(tsv_file, fasta_file, gbk_file, molecule_type='DNA'):
    with open(tsv_file, 'r') as tsv:
        tsv_reader = csv.reader(tsv, delimiter='\t')
        
        for record in SeqIO.parse(fasta_file, 'fasta'):
            fasta_name = fasta_file.split('/')[-1].split('.')[0]
            seq_record = SeqRecord(Seq(str(record.seq)), id=fasta_name, description=record.description)

        seq_record.annotations['molecule_type'] = molecule_type
        
        for row in tsv_reader:
            if row[0].startswith('#'):
                continue
            feature_type = row[1]
            start = int(row[2]) - 1  # Convert to 0-based index
            end = int(row[3])
            strand = 1 if row[4] == '+' else -1 if row[4] == '-' else 0
            locus_tag = row[5]
            gene = row[6]
            product = row[7]
            db_xrefs = row[8]
            
            qualifiers = {
                'locus_tag': locus_tag,
                'gene': gene,
                'product': product,
                'db_xref': db_xrefs.split(', '),
            }
            if feature_type == 'cds':
                sequence = seq_record.seq[start:end]
                if strand == -1:
                    sequence = sequence.reverse_complement()
                translation = str(sequence.translate(to_stop=False))
                qualifiers['translation'] = translation
                
                gene_qualifiers = {
                    'locus_tag': locus_tag,
                    'gene': gene
                }
                gene_feature = SeqFeature(FeatureLocation(start, end, strand=strand), type='gene', qualifiers=gene_qualifiers)
                seq_record.features.append(gene_feature)
            
            feature = SeqFeature(FeatureLocation(start, end, strand=strand), type=feature_type, qualifiers=qualifiers)
            seq_record.features.append(feature)

        with open(gbk_file, 'w') as gbk:
            SeqIO.write(seq_record, gbk, 'genbank')

def parse_args():
    parser = argparse.ArgumentParser(description='Convert a TSV file to a GenBank file')
    parser.add_argument('-i', dest='prefix', required=True, help='Input TSV file')
    
    return parser.parse_args()  # Returns an instance of the class ArgumentParser

def main():
    args = parse_args()
    fasta = f"{args.prefix}.fna"
    inputt = f"{args.prefix}.final.tsv"
    s = args.prefix
    # Divide el string por '/' y lo guarda en una lista
    path_list = s.split('/')
    # Cambia el último elemento de la lista
    path_list[-1] = path_list[-2].split('_')[1] + '_' + path_list[-2].split('_')[0]
    # Une la lista nuevamente en un string con '/'
    new_path = '/'.join(path_list)
    print(new_path)  # Output: "/storage/PGO/results/mtb/proj_anc2023/bakta/res/G49171_L4/L4_G49171"

    output = f"{new_path}.gbk"
    directory_path = os.path.dirname(args.prefix)
    os.chdir(directory_path)
    tsv_to_gbk(inputt, fasta, args.gbk_file, output)

    destination_path = "/storage/PGO/results/mtb/proj_anc2023/bakta/fix_gbk/all_gbk"

    # use the move function to move the file
    shutil.move(output, destination_path)
   

if __name__ == '__main__':
    main()