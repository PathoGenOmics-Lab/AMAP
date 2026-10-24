#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The proteins PANGENE searches every genome with, and the matrix it writes.

pangene builds its graph from one set of proteins aligned to every genome, the
way its documentation describes: a gene is then one row and each genome one
column, its own included. AMAP's set is the reference proteome, each protein
named after its gene (Rv0001 for Rv0001_1-1524), and the genes of each sample
that no reference gene matches (an empty `rvs` in <sample>.AMAP.tsv), named
after the sample's gene. A sequence already in the set is not added again.

    pangenome.py proteins --reference REF.faa --out proteins.faa SAMPLE.AMAP.tsv [...]

reads each sample's proteins from the <sample>.AMAP.faa next to its table.

    pangenome.py matrix --out pangenome.Rtab GENES_BY_GENOME.Rtab GENOME.fasta [...]

turns what `pangene.js gfa2matrix` writes into the matrix AMAP keeps: one
column per genome, named after its file (<name>.fasta), where pangene names it
<name>#0, and a column of zeros for a genome none of the proteins was found
in, which pangene leaves out.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

# Rv0001_1-1524: a reference gene and where the reference has it
COORDS = re.compile(r'_\d+-\d+$')
# pangene names a genome after its PAF and the haplotype: <name>#0
HAPLOTYPE = re.compile(r'#\d+$')


def read_fasta(path):
    """(first word of the header, sequence) for each record. The H37Rv
    proteome in the container has Windows line ends."""
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                words = line[1:].split()
                name, seq = (words[0] if words else ''), []
            elif line:
                seq.append(line)
    if name is not None:
        yield name, ''.join(seq)


def gene_name(name):
    """pangene reads a protein's name as <gene>:<protein> and takes the gene
    up to the first colon, so a colon of the gene's own would cut it short."""
    return name.replace(':', '_')


def genes_outside_reference(table):
    """IDs of the genes of a <sample>.AMAP.tsv that no reference gene matches."""
    with open(table, newline='') as fh:
        rows = csv.DictReader(fh, delimiter='\t')
        if 'id' not in (rows.fieldnames or []) or 'rvs' not in rows.fieldnames:
            sys.exit(f"{table}: no id and rvs columns, as FIX_CDS writes them")
        return {row['id'] for row in rows if not (row['rvs'] or '').strip()}


def protein_file(table):
    """<sample>.AMAP.faa, next to <sample>.AMAP.tsv."""
    table = Path(table)
    return table.with_name(table.name[:-len('.tsv')] + '.faa')


def write_proteins(reference, tables, out):
    """Write the set; return (reference proteins, sample proteins, repeats left out)."""
    seen, genes = set(), set()
    n_ref = n_own = n_repeat = 0
    with open(out, 'w') as fh:
        for name, seq in read_fasta(reference):
            seq = seq.rstrip('*')
            if not seq:
                continue
            gene = gene_name(COORDS.sub('', name) or name)
            fh.write(f">{gene}:{name}\n{seq}\n")
            seen.add(seq)
            genes.add(gene)
            n_ref += 1
        for table in sorted(tables):
            outside = genes_outside_reference(table)
            if not outside:
                continue
            for name, seq in read_fasta(protein_file(table)):
                gene_id = name.split(';')[0]
                seq = seq.rstrip('*')
                if gene_id not in outside or not seq:
                    continue
                if seq in seen:
                    n_repeat += 1
                    continue
                gene = gene_name(gene_id)
                if gene in genes:
                    # Never two genes under one name: pangene would take
                    # them for forms of the same gene.
                    k = 2
                    while f"{gene}_{k}" in genes:
                        k += 1
                    gene = f"{gene}_{k}"
                fh.write(f">{gene}:{name}\n{seq}\n")
                seen.add(seq)
                genes.add(gene)
                n_own += 1
    return n_ref, n_own, n_repeat


def genome_name(path):
    """A genome's name, as PANGENE names its alignments: the file's name
    without .fasta."""
    name = Path(path).name
    return name[:-len('.fasta')] if name.endswith('.fasta') else name


def write_matrix(raw, genomes, out):
    """Write the matrix; return (genes, genomes, genes in every genome)."""
    with open(raw) as fh:
        header = fh.readline().rstrip('\r\n').split('\t')
        rows = [line.rstrip('\r\n').split('\t') for line in fh if line.strip()]
    # gfa2matrix with no genome at all writes "Gene" and a tab
    cols = [HAPLOTYPE.sub('', c) for c in header[1:]]
    names = sorted({genome_name(g) for g in genomes} | {c for c in cols if c})
    where = {name: i for i, name in enumerate(names)}
    n_core = 0
    with open(out, 'w') as fh:
        fh.write('\t'.join(['Gene'] + names) + '\n')
        for row in rows:
            values = [0] * len(names)
            for col, value in zip(cols, row[1:]):
                if col:
                    i = where[col]
                    values[i] = max(values[i], 1 if value.strip() not in ('', '0') else 0)
            n_core += all(values)
            fh.write('\t'.join([row[0]] + [str(v) for v in values]) + '\n')
    return len(rows), len(names), n_core


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0].strip())
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('proteins', help="Write the proteins every genome is searched with.")
    p.add_argument('--reference', required=True, help="The reference proteome.")
    p.add_argument('--out', required=True, help="The protein set to write.")
    p.add_argument('tables', nargs='*', help="Each sample's <sample>.AMAP.tsv.")
    m = sub.add_parser('matrix', help="Name the matrix's columns after the genomes.")
    m.add_argument('--out', required=True, help="The matrix to write.")
    m.add_argument('raw', help="What `pangene.js gfa2matrix` wrote.")
    m.add_argument('genomes', nargs='*', help="Every genome aligned to, <name>.fasta.")
    args = parser.parse_args(argv)

    if args.command == 'proteins':
        n_ref, n_own, n_repeat = write_proteins(args.reference, args.tables, args.out)
        print(f"{n_ref} reference proteins and {n_own} of genes outside the reference "
              f"({n_repeat} more already in the set)")
    else:
        n_genes, n_genomes, n_core = write_matrix(args.raw, args.genomes, args.out)
        print(f"{n_genes} genes x {n_genomes} genomes, {n_core} in every genome")
    return 0


if __name__ == '__main__':
    sys.exit(main())
