#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd

# Definimos los nombres de las columnas
column_names = ["Sequence Id", "Type", "Start", "Stop", "Strand", "Locus Tag", "Gene", "Product", "DbXrefs"]

# Leemos el primer archivo, ignorando las líneas que comienzan con "#"
df1 = pd.read_csv('G26689_L8.tsv', sep='\t', comment='#', names=column_names)

df2 = pd.read_csv('G26689_L8.cleancoord.tsv', sep='\t', header=None, names=["Locus Tag", "Start", "Stop", "Strand", "Cluster"])
