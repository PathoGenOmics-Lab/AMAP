#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Export a clean proteome: no pseudogenes, no internal stop codons, no empty sequences.

This produces the definitive protein FASTA for downstream analyses that require
valid translations (orthology, phylogenomics, pangenome tools, dN/dS).

Usage:
    export_clean_proteome.py input.faa mastertable.tsv output.clean.faa

Rules:
1. Skip the genes AMAP calls pseudogenes, those whose `pseudogene` is yes in
   the master table, matched on the ID a FASTA header starts with
2. Skip sequences containing internal '*' (stop codons within the ORF)
3. Strip trailing '*' (terminal stop codon, valid but unwanted in most tools)
4. Skip empty sequences (0 amino acids after cleaning)
5. Preserve original FASTA headers
"""

import sys
import os
import pandas as pd


def load_pseudogene_ids(mastertable_path):
    """The IDs of the genes AMAP calls pseudogenes: check_cds.py's `pseudogene`
    column, the flag behind the GFF3's pseudo=true and the GenBank's /pseudo.

    Mycobrowser's Is_Pseudogene was read first, and the master table always
    has it: empty when Mycobrowser is off, and otherwise marking the sample's
    genes named after one of the reference's pseudogenes, whatever AMAP made
    of them. A table without these columns stops the export, where a warning
    let every pseudogene through.
    """
    df = pd.read_csv(mastertable_path, sep='\t', usecols=['simple_id', 'pseudogene'],
                     dtype=str, keep_default_na=False)
    return set(df.loc[df['pseudogene'] == 'yes', 'simple_id'])


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
            # FIX_CDS writes ID;strand;start-end, and the master table has the
            # ID alone, as simple_id
            gene_id = header.split()[0].split(';')[0]

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
