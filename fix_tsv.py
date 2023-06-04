#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd

def main():
    df1 = pd.read_csv('data/2018-01-01.tsv', sep='\t')
    df2 = pd.read_csv('data/2018-01-02.tsv', sep='\t')

if '__name__' == '__main__':
    main()