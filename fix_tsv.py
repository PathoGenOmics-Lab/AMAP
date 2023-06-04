#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd

def main():
    df1 = pd.read_csv('G26689_L8.tsv', sep='\t', header=0)
    df2 = pd.read_csv('G26689_L8.cleancoord.tsv', sep='\t', header=None, names=['Locus Tag', 'Start', 'Stop', 'Strand', 'Cluster'])

if '__name__' == '__main__':
    main()