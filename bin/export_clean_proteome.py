#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Export a clean proteome: no pseudogenes, no internal stop codons, no empty sequences.

This produces the definitive protein FASTA for downstream analyses that require
valid translations (orthology, phylogenomics, pangenome tools, dN/dS).

Usage:
    export_clean_proteome.py input.faa mastertable.tsv output.clean.faa

Rules:
1. Skip sequences flagged as pseudogenes in the master table (Is_Pseudogene == True)
2. Skip sequences containing internal '*' (stop codons within the ORF)
3. Strip trailing '*' (terminal stop codon — valid but unwanted in most tools)
4. Skip empty sequences (0 amino acids after cleaning)
5. Preserve original FASTA headers
"""

import sys
import os
import pandas as pd


def load_pseudogene_ids(mastertable_path):
    """Load set of gene IDs flagged as pseudogenes from the master table."""
    pseudo_ids = set()
    try:
        df = pd.read_csv(mastertable_path, sep='\t', low_memory=False)
        # Check both column name variants
        for col in ['Is_Pseudogene', 'pseudogene', 'is_pseudogene']:
            if col in df.columns:
                pseudo_mask = df[col].astype(str).str.lower().isin(['true', 'yes', '1'])
                id_col = 'simple_id' if 'simple_id' in df.columns else df.columns[0]
                pseudo_ids.update(df.loc[pseudo_mask, id_col].dropna().tolist())
                break
    except Exception as e:
        print(f"WARNING: Could not read pseudogene info from {mastertable_path}: {e}", file=sys.stderr)
    return pseudo_ids


def clean_proteome(input_faa, mastertable_path, output_faa):
    """Filter FAA to produce clean proteome."""
    pseudo_ids = load_pseudogene_ids(mastertable_path)

    total = 0
    kept = 0
    skipped_pseudo = 0
    skipped_stop = 0
    skipped_empty = 0

    with open(input_faa) as fin, open(output_faa, 'w') as fout:
        header = None
        seq_parts = []

        def process_record():
            nonlocal kept, skipped_pseudo, skipped_stop, skipped_empty
            if header is None:
                return
            
            seq = ''.join(seq_parts)
            gene_id = header.split()[0].lstrip('>')

            # Check pseudogene
            if gene_id in pseudo_ids:
                skipped_pseudo += 1
                return

            # Strip trailing stop
            seq = seq.rstrip('*')

            # Check internal stops
            if '*' in seq:
                skipped_stop += 1
                return

            # Check empty
            if len(seq) == 0:
                skipped_empty += 1
                return

            fout.write(f">{header}\n")
            # Write in 70-char lines (standard FASTA)
            for i in range(0, len(seq), 70):
                fout.write(seq[i:i+70] + '\n')
            kept += 1

        for line in fin:
            line = line.strip()
            if line.startswith('>'):
                process_record()
                total += 1
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        
        # Process last record
        process_record()

    # Summary
    print(f"Clean proteome: {kept}/{total} proteins exported "
          f"({skipped_pseudo} pseudogenes, {skipped_stop} internal stops, "
          f"{skipped_empty} empty)", file=sys.stderr)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} input.faa mastertable.tsv output.clean.faa", file=sys.stderr)
        sys.exit(1)
    clean_proteome(sys.argv[1], sys.argv[2], sys.argv[3])
