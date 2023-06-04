#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import argparse
def main():
    args = parse_argparse() 
    file_1 = f"{args.input}.tsv"
    file_2 = f"{args.input}.cleancoord.tsv"  
    # Define column names
    column_names = ["Sequence Id", "Type", "Start", "Stop", "Strand", "Locus Tag", "Gene", "Product", "DbXrefs"]

    # Read first file, ignore lines that start with "#"
    df1 = pd.read_csv(file_1, sep='\t', comment='#', names=column_names)

    df2 = pd.read_csv(file_2, sep='\t', header=None, names=["Locus Tag", "Start", "Stop", "Strand", "Cluster"])
    df2_unique = df2[~df2['Locus Tag'].isin(df1['Locus Tag'])]
    df2_unique['Locus Tag'] = df2_unique['Locus Tag'].str.split('_').str[0]
    df1.set_index('Locus Tag', inplace=True)
    df2.set_index('Locus Tag', inplace=True)
    df2_unique.set_index('Locus Tag', inplace=True)
    df2_unique = df2_unique.drop('Cluster', axis=1)
    df2_unique['Type'] = 'cds'
    df2_unique['Sequence Id'] = 'contig_1'

    df1_not_cds = df1[df1["Type"] != "cds"]

    merged = df1.join(df2, lsuffix='_file1', rsuffix='_file2', how='inner') # Changed 'outer' to 'inner'
    merged['Start'] = np.where(pd.isnull(merged['Start_file2']), merged['Start_file1'], merged['Start_file2'])
    merged['Stop'] = np.where(pd.isnull(merged['Stop_file2']), merged['Stop_file1'], merged['Stop_file2'])
    merged['Strand'] = np.where(pd.isnull(merged['Strand_file2']), merged['Strand_file1'], merged['Strand_file2'])
    merged.drop(columns=['Start_file1', 'Stop_file1', 'Strand_file1', 'Start_file2', 'Stop_file2', 'Strand_file2','Cluster'], inplace=True)
    merged = merged[merged['Type'] == 'cds']
    print(merged)
    # Concatenate df1_not_cds and merged
    final_df = pd.concat([df1_not_cds, merged])


    # Concatenate final_df and df2_unique
    final_df = pd.concat([final_df, df2_unique])
    final_df.sort_values(by='Start', inplace=True)

    final_df = final_df.reset_index()
    final_df = final_df.reindex(columns=["Sequence Id", "Type", "Start", "Stop", "Strand", "Locus Tag", "Gene", "Product", "DbXrefs"])
    # Write to a TSV file
    final_df.to_csv(f"{args.input}.final.tsv"  , sep='\t', index=False)

def parse_argparse():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input', help='Input path', required=True)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    main()