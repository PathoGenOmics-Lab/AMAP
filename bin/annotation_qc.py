#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-sample annotation quality metrics.

Computes:
- Total CDS count, pseudogene count, hypothetical %
- Coding density (% of genome covered by CDS)
- Mean/median gene length
- GC content
- Genome size, N50 (if multi-contig)
- eggNOG annotation rate (% of genes with functional annotation)
- Internal stop codon count

Usage:
    annotation_qc.py mastertable.tsv genome.fasta output.qc.tsv
"""

import sys
import os
import pandas as pd
from collections import Counter

# Allow importing from bin/classify/ when run as a script, as classify_cases.py does.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classify import circular


def parse_fasta_stats(fasta_path):
    """Compute basic assembly metrics from FASTA."""
    contig_lengths = []
    total_bases = Counter()
    
    with open(fasta_path) as f:
        seq_parts = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq_parts:
                    seq = ''.join(seq_parts)
                    contig_lengths.append(len(seq))
                    total_bases.update(seq.upper())
                seq_parts = []
            else:
                seq_parts.append(line)
        if seq_parts:
            seq = ''.join(seq_parts)
            contig_lengths.append(len(seq))
            total_bases.update(seq.upper())

    genome_size = sum(contig_lengths)
    n_contigs = len(contig_lengths)
    
    # N50
    sorted_lens = sorted(contig_lengths, reverse=True)
    running = 0
    n50 = 0
    for l in sorted_lens:
        running += l
        if running >= genome_size / 2:
            n50 = l
            break

    # GC content
    gc = total_bases.get('G', 0) + total_bases.get('C', 0)
    atgc = gc + total_bases.get('A', 0) + total_bases.get('T', 0)
    gc_pct = (gc / atgc * 100) if atgc > 0 else 0

    return {
        'genome_size': genome_size,
        'n_contigs': n_contigs,
        'n50': n50,
        'gc_pct': round(gc_pct, 2)
    }


def compute_qc(mastertable_path, fasta_path, output_path):
    """Compute annotation QC metrics."""
    # Assembly metrics
    asm = parse_fasta_stats(fasta_path)
    
    # Read master table
    try:
        df = pd.read_csv(mastertable_path, sep='\t', low_memory=False)
    except Exception:
        # Empty master table — write minimal QC
        qc = {**asm, 'total_cds': 0}
        pd.DataFrame([qc]).to_csv(output_path, sep='\t', index=False)
        return

    total_cds = len(df)
    
    # Pseudogenes
    pseudo_count = 0
    for col in ['Is_Pseudogene', 'pseudogene', 'is_pseudogene']:
        if col in df.columns:
            pseudo_count = df[col].astype(str).str.lower().isin(['true', 'yes', '1']).sum()
            break

    # Hypothetical proteins
    hypo_count = 0
    for col in ['Product_x', 'Product', 'product', 'Description']:
        if col in df.columns:
            hypo_count = df[col].astype(str).str.lower().str.contains('hypothetical', na=False).sum()
            break
    hypo_pct = round(hypo_count / total_cds * 100, 1) if total_cds > 0 else 0

    # Gene lengths
    start_col = end_col = None
    for s, e in [('start', 'end'), ('Start', 'End'), ('Start_x', 'Stop_x')]:
        if s in df.columns and e in df.columns:
            start_col, end_col = s, e
            break

    if start_col and end_col:
        starts = pd.to_numeric(df[start_col], errors='coerce')
        ends = pd.to_numeric(df[end_col], errors='coerce')
        lengths = ends - starts + 1
        # A gene across the origin of a circular sequence (start after end)
        # runs on from the sequence's last base to its first.
        across = (starts > ends) & (ends > 0)
        if across.any():
            sizes = circular.read_lengths(fasta_path)
            contigs = (df['contig'].astype(str) if 'contig' in df.columns
                       else pd.Series(next(iter(sizes), ''), index=df.index))
            lengths[across] += contigs[across].map(sizes)
        lengths = lengths.dropna()
        mean_len = round(lengths.mean(), 0) if len(lengths) > 0 else 0
        median_len = round(lengths.median(), 0) if len(lengths) > 0 else 0
        total_coding_bp = lengths.sum()
    else:
        mean_len = median_len = total_coding_bp = 0

    coding_density = round(total_coding_bp / asm['genome_size'] * 100, 1) if asm['genome_size'] > 0 else 0

    # eggNOG annotation rate.
    # emapper writes '-' for every field it could not fill, and those survive
    # into the master table (count_annot.py normalises them on read, this did
    # not). notna() therefore counted an unannotated gene as annotated and the
    # rate sat near 100% regardless of how the run had actually gone, which is
    # the opposite of useful in a QC metric.
    eggnog_cols = ['COG_category', 'seed_ortholog', 'Description']
    annotated = 0
    for col in eggnog_cols:
        if col in df.columns:
            values = df[col]
            missing = values.isna() | values.astype(str).str.strip().isin(['', '-'])
            annotated = int((~missing).sum())
            break
    annot_rate = round(annotated / total_cds * 100, 1) if total_cds > 0 else 0

    # Compile QC
    qc = {
        'genome_size': asm['genome_size'],
        'n_contigs': asm['n_contigs'],
        'n50': asm['n50'],
        'gc_pct': asm['gc_pct'],
        'total_cds': total_cds,
        'pseudogene_count': pseudo_count,
        'hypothetical_pct': hypo_pct,
        'mean_gene_length': int(mean_len),
        'median_gene_length': int(median_len),
        'coding_density_pct': coding_density,
        'eggnog_annotation_pct': annot_rate,
    }

    pd.DataFrame([qc]).to_csv(output_path, sep='\t', index=False)
    
    # Print summary
    print(f"QC: {total_cds} CDS, {pseudo_count} pseudo, {hypo_pct}% hypothetical, "
          f"{coding_density}% coding density, {annot_rate}% annotated", file=sys.stderr)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} mastertable.tsv genome.fasta output.qc.tsv", file=sys.stderr)
        sys.exit(1)
    compute_qc(sys.argv[1], sys.argv[2], sys.argv[3])
