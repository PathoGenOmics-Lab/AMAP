#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Explore pangenome')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input JSON file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output JSON file')
    args = parser.parse_args()
    return args

def process_json(json_file):
    with open(json_file, 'r') as json_in:
        json_data = json.load(json_in)
        # Iterar sobre la lista de diccionarios
    for d in json_in:
        # Acceder a los campos que necesites
        gene_id = d['geneId']
        gene_len = d['geneLen']
        count = d['count']
        dupli = d['dupli']

def main():
    args = parse_args()


if __name__ == '__main__':
    main()