#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import argparse
import pandas as pd
def parse_args():
    parser = argparse.ArgumentParser(description='Explore pangenome')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input JSON file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output JSON file')
    args = parser.parse_args()
    return args

def process_json(json_file):
    with open(json_file, 'r') as json_in:
        json_data = json.load(json_in)

        df_list = []  # Inicializar una lista vacía para almacenar los DataFrames temporales

        for d in json_data:
            gene_id = d['geneId']
            locus = d['locus']
            gene_len = d['geneLen']
            dupli = d['dupli']
            ann = d['ann']
            msa = d['msa']
            divers= d['divers']
            event = d['event']
            GName = d['GName']
            locus = d['locus']

            components = [item.split('_')[0] for item in locus.split()]

            counts = {pattern: components.count(pattern) for pattern in ['L5', 'L6', 'L9', 'A1', 'A2', 'A3', 'A4']}
            counts['gene_id'] = gene_id
            counts['gene_len'] = gene_len
            counts['dupli'] = dupli
            counts['ann'] = ann
            counts['msa'] = msa
            counts['divers'] = divers
            counts['event'] = event
            counts['GName'] = GName
            counts['locus'] = locus

            # Crear un DataFrame temporal y agregarlo a la lista
            temp_df = pd.DataFrame([counts])
            df_list.append(temp_df)

        # Concatenar todos los DataFrames temporales en un DataFrame final
        df = pd.concat(df_list, ignore_index=True)

        return df



def main():
    args = parse_args()
    #process_json(args.input)
    df = process_json(args.input)
    df.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    main()