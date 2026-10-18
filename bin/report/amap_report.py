#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Write amap_report.html: AMAP's interactive report, in one self-contained file.

It is organised by the questions a run answers, not by the kind of data:
did each genome annotate well, where on the reference genes are missing or
broken, what happened to one gene in every genome, what the correction did,
which genes are pseudogenes, what lies outside the reference, and how the
functional profile differs between lineages. Methods and provenance close it.

Everything the page shows is computed here from the pipeline's own outputs and
embedded as one compressed JSON document; the page's JavaScript draws it. No
file the page needs lives anywhere else, so it opens with no network and
survives being emailed.

Inputs, one folder per kind (the file names are the ones the pipeline writes):

    --amap DIR      <name>.AMAP.tsv, the corrected annotation (FIX_CDS)
    --cases DIR     <name>-summary.tsv, <name>-case1..4.tsv, <name>-rvmissing.tsv
                    (TYPE_CASES)
    --bakta DIR     Bakta's folder of each sample: products (<name>.tsv), genome
                    size, contigs, GC and N50 (<name>.txt), contig lengths (.gff3)
    --genomes DIR   <name>.fasta, the genomes, instead of Bakta's statistics
    --master DIR    <name>-AMAPmastertable.tsv, with functional annotation
    --ref-proteome  the reference proteome, headers NAME_START-END

and single files: the Mycobrowser table, the merged pangenome matrix, the
homologous communities, params.yml, software_versions.tsv and the sample IDs
of the sample sheet (see --help). Only --amap and --ref-proteome are required. Whatever else is missing, the
section that needs it says so instead of drawing anything.
"""

import argparse
import base64
import csv
import json
import os
import re
import statistics
import sys
import zlib
from collections import Counter, OrderedDict, defaultdict
from datetime import datetime

__author__ = 'Paula Ruiz-Rodriguez'
__license__ = 'GPL 3'
__lab__ = 'PathoGenOmics (I2SysBio)'

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATE = os.path.join(HERE, 'amap_report.html')
STYLES = os.path.join(HERE, 'amap_report.css')
SCRIPT = os.path.join(HERE, 'amap_report.js')
LOGO = os.path.join(HERE, 'amap-logo.png')

# Some master table cells (GO terms) are long; csv's default limit is 128 kB.
csv.field_size_limit(2 ** 31 - 1)

# A reference protein's header: its name, then where the reference has it.
HEADER = re.compile(r'^(.+)_(\d+)-(\d+)$')

# What one reference gene is in one genome. Intact genes are left out of the
# embedded lists; every other state is one letter.
MISSING, PSEUDO, PARTIAL, COPIES, LENGTH, SHARED = 'M', 'P', 'T', 'C', 'L', 'S'

# Where each final gene came from, as the case tables and the final table
# say: a Bakta call kept, by its case (1 to 4), or a gene placed from the
# reference alignment, over a Bakta call that disagreed (R) or where Bakta
# called nothing (A). Without case tables, only the kind of identifier is
# known: a Bakta locus tag (b) or a reference gene's name (r).

# Lineages in the order MTBC phylogenies draw them; anything else follows,
# alphabetically, so a custom reference's groups still sort.
LINEAGE_ORDER = ['L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8', 'L9', 'L10',
                 'A1', 'A2', 'A3', 'A4', 'La1', 'La2', 'La3']
NO_LINEAGE = ('NA', 'Untyped', '')

# How far from the other genomes a value has to be, in robust standard
# deviations (1.4826 x the median absolute deviation), to be flagged.
OUTLIER_Z = 3.5
# A lineage is its own comparison group from this many genomes on; a smaller
# one is compared with every genome in the run.
GROUP_MIN = 5
# The QC values compared, the direction that is bad, and the smallest gap from
# the median worth a flag, so that identical genomes (deviation 0) do not flag
# a difference of one gene.
# The floor is a fraction of the median ('rel') or a plain difference ('abs').
QC_METRICS = [
    ('size', 'both', ('rel', 0.01), 'genome size'),
    ('contigs', 'high', ('abs', 5), 'contigs'),
    ('n50', 'low', ('rel', 0.05), 'N50'),
    ('found_pct', 'low', ('abs', 0.5), 'reference genes found'),
    ('missing', 'high', ('abs', 5), 'missing reference genes'),
    ('pseudo', 'high', ('abs', 5), 'pseudogenes'),
    ('partial', 'high', ('abs', 5), 'partial genes'),
    ('noref', 'high', ('abs', 5), 'genes not in the reference'),
]


# --------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------

def _rows(path):
    """The rows of a TSV with a header line, as dicts of stripped strings."""
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            yield {(k or '').strip(): (v or '').strip() for k, v in row.items() if k is not None}


def _int(value, default=None):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _files(folder, suffix):
    """{sample name: path} of the files ending with suffix anywhere under
    folder: staged flat in a task, or in a folder per sample as the pipeline
    publishes them (results/AMAP/<name>/, results/BAKTA/<name>/)."""
    found = {}
    if not folder or not os.path.isdir(folder):
        return found
    for root, dirs, names in os.walk(folder, followlinks=True):
        dirs.sort()
        for fname in sorted(names):
            if fname.endswith(suffix) and len(fname) > len(suffix):
                found.setdefault(fname[:-len(suffix)], os.path.join(root, fname))
    return found


def read_proteome(path):
    """The reference genes, from the proteome's headers.

    Returns a list of dicts (name, start, end, aa) in reference order and a
    list of warnings. A header without coordinates still names a gene, which
    then has no place on the map.
    """
    genes, seen, twice, warnings = [], set(), [], []
    name = start = end = None
    length = 0

    def close():
        if name is None:
            return
        if name in seen:
            twice.append(name)
            return
        seen.add(name)
        genes.append({'name': name, 'start': start, 'end': end, 'aa': length})

    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                close()
                header = line[1:].strip().split()[0] if line[1:].strip() else ''
                m = HEADER.match(header)
                if m:
                    name, start, end = m.group(1), int(m.group(2)), int(m.group(3))
                    if start > end:
                        # Across the origin of a circular reference: no one
                        # place on a straight axis
                        start = end = None
                else:
                    name, start, end = header, None, None
                length = 0
            else:
                length += len(line.strip().rstrip('*'))
    close()
    if twice:
        shown = ', '.join(twice[:5]) + (', ...' if len(twice) > 5 else '')
        warnings.append(f'The reference proteome names {len(twice)} gene{"s" if len(twice) > 1 else ""} '
                        f'more than once ({shown}); the report places each where its first protein is.')
    genes.sort(key=lambda g: (g['start'] is None, g['start'] or 0, g['end'] or 0, g['name']))
    return genes, warnings


def read_amap_table(path):
    """A sample's corrected annotation, <name>.AMAP.tsv, in its own order."""
    genes = []
    for r in _rows(path):
        start, end = _int(r.get('start')), _int(r.get('end'))
        if start is None or end is None:
            continue
        rvs = [x.strip() for x in r.get('rvs', '').split(',') if x.strip()]
        # A gene across the origin of a circular sequence keeps its start
        # after its end, as AMAP and Bakta write it (see classify/circular.py)
        genes.append({
            'id': r.get('id', ''),
            'strand': '-' if r.get('strand') == '-' else '+',
            'start': start,
            'end': end,
            'rvs': rvs,
            'pseudo': r.get('pseudogene', '').lower() in ('yes', 'true', '1'),
            'contig': r.get('contig', ''),
            'partial': r.get('partial', ''),
            'note': r.get('note', ''),
        })
    return genes


def read_summary(path):
    """A case summary: the counts, and the sample ID and lineage it names."""
    counts, sample, lineage, sequence = {}, None, None, None
    for r in _rows(path):
        n = _int(r.get('Number'))
        if r.get('Type') and n is not None:
            counts[r['Type']] = n
        sequence = sequence or r.get('Sequence')
        sample = sample or r.get('Sample')
        lineage = lineage if lineage is not None else r.get('Lineage')
    return {'counts': counts, 'sample': sample, 'lineage': lineage, 'sequence': sequence}


def read_case_table(path):
    """The Bakta calls of one case: locus tag, coordinates, strand, contig."""
    calls = []
    for r in _rows(path):
        start, end = _int(r.get('Start')), _int(r.get('End'))
        if not r.get('Locus Tag') or start is None or end is None:
            continue
        calls.append({'tag': r['Locus Tag'], 'start': start, 'end': end,
                      'strand': '-' if r.get('Strand') == '-' else '+',
                      'contig': r.get('Contig', ''), 'prot': r.get('Prot Id', '')})
    return calls


def read_rvmissing(path):
    """The reference proteins miniprot placed where Bakta called nothing."""
    names = set()
    for r in _rows(path):
        prot = HEADER.sub(r'\1', r.get('Prot Id', ''))
        if prot:
            names.add(prot)
    return names


def read_fasta_stats(path):
    """Contigs (name, length) in file order and the GC percentage."""
    contigs, gc, acgt = [], 0, 0
    name, length = None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    contigs.append((name, length))
                header = line[1:].strip()
                name, length = (header.split()[0] if header else ''), 0
            else:
                s = line.strip().upper()
                length += len(s)
                g, c = s.count('G'), s.count('C')
                gc += g + c
                acgt += g + c + s.count('A') + s.count('T')
    if name is not None:
        contigs.append((name, length))
    return contigs, (round(100.0 * gc / acgt, 2) if acgt else None)


def n50(lengths):
    total, running = sum(lengths), 0
    for n in sorted(lengths, reverse=True):
        running += n
        if running * 2 >= total:
            return n
    return 0


def read_bakta_stats(path):
    """Genome size, contigs, GC and N50 from Bakta's summary, <name>.txt."""
    keys = {'Length': 'size', 'Count': 'contigs', 'GC': 'gc', 'N50': 'n50'}
    out = {}
    with open(path) as fh:
        for line in fh:
            if ':' not in line:
                continue
            key, _, value = line.partition(':')
            if key.strip() in keys and keys[key.strip()] not in out:
                try:
                    v = float(value.strip())
                except ValueError:
                    continue
                out[keys[key.strip()]] = v if key.strip() == 'GC' else int(v)
            if line.startswith('Annotation'):
                break
    return out


def read_gff_regions(path):
    """Contigs (name, length) from a GFF3's ##sequence-region lines, which
    Bakta writes before each contig's features rather than all at the top."""
    contigs = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('##sequence-region'):
                f = line.split()
                if len(f) >= 4:
                    contigs.append((f[1], int(f[3]) - int(f[2]) + 1))
            elif line.startswith('##FASTA'):
                break
    return contigs


def read_bakta_products(path):
    """{locus tag: (gene, product)} of Bakta's CDS."""
    products = {}
    with open(path, newline='') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) >= 8 and f[1].lower() == 'cds':
                products[f[5]] = (f[6], f[7])
    return products


MASTER_COLUMNS = ('COG_category', 'Description', 'Preferred_name', 'KEGG_Pathway',
                  'Product_x', 'Product', 'Gene')


def read_master(path):
    """{gene id: the functional columns the report uses} of a master table."""
    table = {}
    for r in _rows(path):
        gid = r.get('id')
        if gid:
            table[gid] = {k: r.get(k, '') for k in MASTER_COLUMNS if r.get(k, '') not in ('', '-')}
    return table


def read_mycobrowser(path):
    """{locus: name, product, category, reference pseudogene, drugs}."""
    table = {}
    for r in _rows(path):
        locus = r.get('Locus')
        if not locus or r.get('Feature', 'CDS') not in ('CDS', ''):
            continue
        table[locus] = {
            'name': r.get('Name', ''),
            'product': r.get('Product', ''),
            'category': r.get('Functional_Category', ''),
            'pseudo': r.get('Is_Pseudogene', ''),
            'drugs': r.get('Drug Resistance Mutations', ''),
        }
    return table


def read_versions(path):
    """[[tool, version], ...] from SOFTWARE_VERSIONS' table."""
    out = []
    for r in _rows(path):
        tool, version = r.get('Tool', ''), r.get('Version', '')
        if tool:
            out.append([tool, version])
    return out


def read_params(path):
    """The parameters and the run's header comments from DUMP_PARAMS.

    Read line by line rather than as YAML: before AMAP wrote each value as a
    JSON string under `params:`, DUMP_PARAMS indented only the first parameter
    there, and those files should still read.
    """
    params, run = OrderedDict(), OrderedDict()
    in_params = False
    with open(path) as fh:
        for line in fh:
            text = line.strip()
            if text.startswith('#'):
                m = re.match(r'#\s*([A-Za-z ]+):\s*(.*)$', text)
                if m:
                    run[m.group(1).strip().lower()] = m.group(2).strip()
                continue
            if text == 'params:':
                in_params = True
                continue
            m = re.match(r'([\w.-]+):\s*(.*)$', text)
            if in_params and m:
                value = m.group(2).strip()
                if len(value) >= 2 and value[0] == value[-1] == '"':
                    try:
                        value = json.loads(value)
                    except ValueError:
                        value = value[1:-1]
                params[m.group(1)] = value
    return params, run


def read_rtab(path):
    """The merged pangene matrix: gene names, column names and 0/1 rows."""
    with open(path, newline='') as fh:
        lines = [l.rstrip('\n') for l in fh if l.strip()]
    if len(lines) < 2:
        return None
    header = lines[0].split('\t')
    columns = header[1:]
    rows = []
    for line in lines[1:]:
        f = line.split('\t')
        rows.append((f[0], [1 if v.strip() not in ('', '0') else 0 for v in f[1:len(header)]]))
    return {'columns': columns, 'rows': rows} if columns else None


def read_communities(path):
    """The homologous communities GET_ALLDUP found across genomes."""
    rows = list(_rows(path))
    return rows


# --------------------------------------------------------------------------
# One genome
# --------------------------------------------------------------------------

def gene_status(hits, ref_len):
    """The state of one reference gene in one genome, from the genes naming it.

    Genes that name the reference gene alone decide it, the best first: not a
    pseudogene, then the longest. A gene that also names other reference genes
    (an IS6110 copy names all seventeen) cannot tell them apart, so a gene seen
    only that way is SHARED. Returns (state, detail, best gene).
    """
    if not hits:
        return MISSING, '', None
    own = [h for h in hits if len(h['rvs']) == 1]
    if not own:
        return SHARED, str(len(hits)), hits[0]
    best = sorted(own, key=lambda h: (h['pseudo'], -(gene_length(h) or 0)))[0]
    if best['pseudo']:
        return PSEUDO, '', best
    if best['partial']:
        return PARTIAL, best['partial'].replace("'", ''), best
    if len(own) > 1:
        return COPIES, str(len(own)), best
    length = gene_length(best)
    if ref_len and length and length != ref_len:
        return LENGTH, str(length - ref_len), best
    return '', '', best


def gene_length(g):
    """Bases of a gene, across the origin of its circular sequence too, which
    needs the sequence's length (None when it is not known)."""
    if g['start'] <= g['end']:
        return g['end'] - g['start'] + 1
    n = g.get('contig_length')
    return n - g['start'] + 1 + g['end'] if n else None


class Strings:
    """A table of repeated strings (notes, products), referred to by index."""

    def __init__(self):
        self.items, self.index = [], {}

    def __call__(self, text):
        if text not in self.index:
            self.index[text] = len(self.items)
            self.items.append(text)
        return self.index[text]


def rle(values):
    """Run-length encode a list: 'v' or 'v*n', comma separated."""
    out, i, n = [], 0, len(values)
    while i < n:
        j = i
        while j + 1 < n and values[j + 1] == values[i]:
            j += 1
        k = j - i + 1
        out.append(f'{values[i]}*{k}' if k > 1 else f'{values[i]}')
        i = j + 1
    return ','.join(out)


def unrle(text):
    """The inverse of rle, as a list of strings."""
    out = []
    if not text:
        return out
    for token in text.split(','):
        value, _, count = token.rpartition('*') if '*' in token else (token, '', '1')
        out.extend([value] * int(count))
    return out


def case_of_calls(cases):
    """{locus tag: case} of every Bakta call. A call can sit in two case tables
    (case 2 and case 3, when it matches one placement exactly and another
    in part); it counts once, under the better case."""
    out = {}
    for case in ('2', '3', '4', '1'):
        for call in cases.get(case, []):
            out.setdefault(call['tag'], case)
    return out


def origins_of(name, genes, case_of, rvmissing):
    """Where each final gene came from: its Bakta call's case, or the reference.

    With no case tables the origin is unknown, and only the kind of identifier
    is kept: 'b' for a Bakta locus tag, 'r' for one named after the reference.
    """
    tag = re.compile('^' + re.escape(name) + r'_\d+$')
    out = []
    for g in genes:
        if g['id'] in case_of:
            out.append(case_of[g['id']])
        elif not case_of:
            out.append('b' if tag.match(g['id']) else ('r' if g['rvs'] else '?'))
        elif any(r in rvmissing for r in g['rvs']):
            out.append('A')
        elif g['rvs']:
            out.append('R')
        else:
            out.append('?')
    return out


BAKTA_ORIGINS = ('1', '2', '3', '4', 'b')


def encode_genes(name, genes, contigs, ref, ref_strand, origins, notes, products):
    """A genome's genes, in their own order, as compact columns.

    Most genes sit where their neighbours predict: a genome is collinear with
    its reference over long stretches, forward or reversed, so the next gene's
    start follows from the previous one's and the reference's own spacing.
    Only the residual from that prediction is stored, and runs of equal
    residuals are run-length encoded: a closed genome's 4000 genes take a few
    kilobytes. amap_report.js decodes it (decodeGenes); the two must agree,
    and decode_genes below is the same decoder, for the tests.
    """
    cidx = {c: i for i, c in enumerate(contigs)}
    col_c, col_r, col_s, col_l, col_t, col_b = [], [], [], [], [], []
    extra, ids, pseudo, partial, prod, wrap = {}, {}, {}, {}, {}, {}
    last = {}            # contig -> (ref index, start) of its last single-reference gene
    prev_start = {}      # contig -> start of its last gene
    # A Bakta call's number follows the last one's, one further on for each
    # call the reference replaced in between; only the difference is stored.
    prev_num, replaced = 0, 0
    tag = re.compile('^' + re.escape(name) + r'_(\d+)$')
    for row, g in enumerate(genes):
        key = str(row)
        c = cidx.setdefault(g['contig'], len(cidx))
        if c == len(contigs):
            contigs.append(g['contig'])
        col_c.append(c)
        single = len(g['rvs']) == 1 and g['rvs'][0] in ref['index']
        ri = ref['index'][g['rvs'][0]] if single else -1
        if len(g['rvs']) > 1 or (g['rvs'] and not single):
            extra[key] = [ref['index'].get(r, r) for r in g['rvs']]
        length = g['end'] - g['start'] + 1
        if g['start'] > g['end']:
            # Across the origin: its end is written out, its length taken as
            # the reference's
            wrap[key] = g['end']
            length = (ref['end'][ri] - ref['start'][ri] + 1) if ri >= 0 else 0
        reverse = False
        if ri >= 0 and c in last:
            pr, ps = last[c]
            col_r.append(str(ri - pr))
            if ri > pr:
                pred = ps + ref['start'][ri] - ref['start'][pr]
            elif ri < pr:
                pred = ps + ref['end'][pr] - ref['end'][ri]
                reverse = True
            else:
                pred = ps
        else:
            col_r.append(f'={ri}' if ri >= 0 else 'x')
            pred = prev_start.get(c, 0)
        col_s.append(g['start'] - pred)
        col_l.append(length - (ref['end'][ri] - ref['start'][ri] + 1) if ri >= 0 else length)
        plus = g['strand'] == '+'
        if ri >= 0:
            col_t.append(0 if plus == (ref_strand[ri] != reverse) else 1)
        else:
            col_t.append(1 if plus else 0)
        # The identifier: a Bakta call's locus tag from the running number, a
        # gene placed from the reference as <name>_<its gene>, anything else
        # written out.
        m = tag.match(g['id'])
        if origins[row] in BAKTA_ORIGINS and m and g['id'] == f'{name}_{int(m.group(1)):05d}':
            num = int(m.group(1))
            col_b.append(num - prev_num - 1 - replaced)
            prev_num, replaced = num, 0
        elif origins[row] not in BAKTA_ORIGINS and g['rvs'] and g['id'] == f"{name}_{g['rvs'][0]}":
            pass
        else:
            ids[key] = g['id'][len(name) + 1:] if g['id'].startswith(name + '_') else '=' + g['id']
        if origins[row] == 'R':
            replaced += 1
        if g['pseudo']:
            pseudo[key] = notes(g['note'])
        if g['partial']:
            partial[key] = g['partial']
        if not g['rvs'] and g.get('product'):
            prod[key] = products(g['product'])
        if ri >= 0:
            last[c] = (ri, g['start'])
        prev_start[c] = g['start']
    return {
        'n': len(genes), 'c': rle(col_c), 'r': rle(col_r), 's': rle(col_s),
        'l': rle(col_l), 't': rle(col_t), 'o': rle(origins), 'b': rle(col_b),
        'x': extra, 'id': ids, 'p': pseudo, 'q': partial, 'pr': prod, 'w': wrap,
    }


def decode_genes(name, enc, ref, ref_strand, contigs):
    """encode_genes backwards, for the tests: the JavaScript does the same."""
    cols = {k: unrle(enc[k]) for k in ('c', 'r', 's', 'l', 't', 'o', 'b')}
    numbers = iter(cols['b'])
    genes, last, prev_start, prev_num, replaced = [], {}, {}, 0, 0
    for row in range(enc['n']):
        key = str(row)
        c = int(cols['c'][row])
        token = cols['r'][row]
        reverse = False
        if token == 'x':
            ri = -1
            pred = prev_start.get(c, 0)
        elif token.startswith('='):
            ri = int(token[1:])
            pred = prev_start.get(c, 0)
        else:
            pr, ps = last[c]
            ri = pr + int(token)
            if ri > pr:
                pred = ps + ref['start'][ri] - ref['start'][pr]
            elif ri < pr:
                pred = ps + ref['end'][pr] - ref['end'][ri]
                reverse = True
            else:
                pred = ps
        start = pred + int(cols['s'][row])
        length = int(cols['l'][row]) + (ref['end'][ri] - ref['start'][ri] + 1 if ri >= 0 else 0)
        bit = int(cols['t'][row])
        plus = (ref_strand[ri] != reverse) == (bit == 0) if ri >= 0 else bit == 1
        rvs = [ref['names'][i] if isinstance(i, int) else i for i in enc['x'][key]] \
            if key in enc['x'] else ([ref['names'][ri]] if ri >= 0 else [])
        origin = cols['o'][row]
        if key in enc['id']:
            literal = enc['id'][key]
            gid = literal[1:] if literal.startswith('=') else f'{name}_{literal}'
        elif origin in BAKTA_ORIGINS:
            prev_num += 1 + replaced + int(next(numbers))
            replaced = 0
            gid = f'{name}_{prev_num:05d}'
        else:
            gid = f'{name}_{rvs[0]}'
        if origin == 'R':
            replaced += 1
        end = enc['w'][key] if key in enc.get('w', {}) else start + length - 1
        genes.append({'id': gid, 'contig': contigs[c], 'start': start, 'end': end,
                      'strand': '+' if plus else '-', 'rvs': rvs, 'origin': origin,
                      'pseudo': key in enc['p'], 'partial': enc['q'].get(key, '')})
        if ri >= 0:
            last[c] = (ri, start)
        prev_start[c] = start
    return genes


def encode_gone(name, genes, cases, case_of, final_ids):
    """The Bakta calls the final annotation does not keep, by locus tag
    number, each next to the final gene that took its place, if one overlaps
    it on its strand: then only how far its ends are from that gene's is
    stored. amap_report.js decodes it (decodeGone).
    """
    by_place = defaultdict(list)
    for row, g in enumerate(genes):
        by_place[(g['contig'], g['strand'])].append((g['start'], g['end'], row))
    for v in by_place.values():
        v.sort()
    tag = re.compile('^' + re.escape(name) + r'_(\d+)$')
    calls = []
    for case in ('2', '3', '4', '1'):
        for call in cases.get(case, []):
            if call['tag'] in final_ids or case_of.get(call['tag']) != case:
                continue
            m = tag.match(call['tag'])
            calls.append((int(m.group(1)) if m else None, case, call))
    calls.sort(key=lambda x: (x[0] is None, x[0] or 0, x[2]['tag']))
    col_case, col_num, col_row, col_s, col_e, loose, literal = [], [], [], [], [], {}, {}
    prev_num = prev_row = 0
    for k, (num, case, call) in enumerate(calls):
        col_case.append(case)
        if num is None:
            literal[str(k)] = call['tag']
        else:
            col_num.append(num - prev_num)
            prev_num = num
        best, best_ov = -1, 0
        for s, e, row in by_place.get((call['contig'], call['strand']), []):
            if s > call['end']:
                break
            ov = min(e, call['end']) - max(s, call['start']) + 1
            if ov > best_ov:
                best, best_ov = row, ov
        if best >= 0:
            col_row.append(best - prev_row)
            prev_row = best
            col_s.append(call['start'] - genes[best]['start'])
            col_e.append(call['end'] - genes[best]['end'])
        else:
            col_row.append('x')
            col_s.append(call['start'])
            col_e.append(call['end'])
            loose[str(k)] = [call['strand'], call['contig']]
    return {'n': len(calls), 'k': rle(col_case), 'num': rle(col_num), 'row': rle(col_row),
            's': rle(col_s), 'e': rle(col_e), 'loose': loose, 'tag': literal}


def decode_gone(name, enc, genes):
    """encode_gone backwards, for the tests."""
    cols = {k: unrle(enc[k]) for k in ('k', 'num', 'row', 's', 'e')}
    numbers = iter(cols['num'])
    out, prev_num, prev_row = [], 0, 0
    for k in range(enc['n']):
        key = str(k)
        if key in enc['tag']:
            tag = enc['tag'][key]
        else:
            prev_num += int(next(numbers))
            tag = f'{name}_{prev_num:05d}'
        if cols['row'][k] == 'x':
            strand, contig = enc['loose'][key]
            out.append({'tag': tag, 'case': cols['k'][k], 'row': -1, 'start': int(cols['s'][k]),
                        'end': int(cols['e'][k]), 'strand': strand, 'contig': contig})
        else:
            prev_row += int(cols['row'][k])
            g = genes[prev_row]
            out.append({'tag': tag, 'case': cols['k'][k], 'row': prev_row,
                        'start': g['start'] + int(cols['s'][k]), 'end': g['end'] + int(cols['e'][k]),
                        'strand': g['strand'], 'contig': g['contig']})
    return out


# --------------------------------------------------------------------------
# The run
# --------------------------------------------------------------------------

def lineage_key(lineage):
    if lineage in LINEAGE_ORDER:
        return (0, LINEAGE_ORDER.index(lineage), lineage)
    if lineage in NO_LINEAGE:
        return (2, 0, lineage)
    return (1, 0, lineage)


def robust_flags(samples):
    """Flag the QC values far from those of comparable genomes.

    A genome is compared with its own lineage when that lineage has
    GROUP_MIN genomes or more, and with the whole run otherwise; a run of
    fewer than four genomes is too small to call anything unusual.
    """
    by_lineage = defaultdict(list)
    for s in samples:
        by_lineage[s['lineage']].append(s)
    for s in samples:
        s['flags'] = []
    if len(samples) < 4:
        return
    for s in samples:
        group = by_lineage[s['lineage']] if len(by_lineage[s['lineage']]) >= GROUP_MIN else samples
        others = [o for o in group if o is not s]
        label = f"the other {s['lineage']} genomes" if group is not samples else 'the other genomes'
        for key, direction, floor, what in QC_METRICS:
            x = s.get(key)
            values = [o.get(key) for o in others if o.get(key) is not None]
            if x is None or len(values) < 3:
                continue
            med = statistics.median(values)
            mad = statistics.median([abs(v - med) for v in values])
            gap = x - med
            kind, amount = floor
            minimum = amount * abs(med) if kind == 'rel' else amount
            if abs(gap) < max(minimum, 1e-9):
                continue
            scale = 1.4826 * mad
            z = gap / scale if scale > 0 else float('inf') * (1 if gap > 0 else -1)
            if (direction == 'high' and z < OUTLIER_Z) or (direction == 'low' and z > -OUTLIER_Z) \
                    or (direction == 'both' and abs(z) < OUTLIER_Z):
                continue
            s['flags'].append({
                'key': key, 'what': what, 'value': x, 'median': med,
                'z': None if scale == 0 else round(z, 1), 'group': label})


def deletion_blocks(ref, status_by_sample, min_genes=2):
    """Runs of consecutive reference genes missing from a genome, grouped
    across genomes where they overlap: regions of difference and lineage
    deletions show up as such blocks.
    """
    placed = [i for i, s in enumerate(ref['start']) if s is not None]
    order = {gi: k for k, gi in enumerate(placed)}
    runs = []
    for sample, status in status_by_sample.items():
        missing = sorted(order[i] for i, (code, _) in status.items() if code == MISSING and i in order)
        start = prev = None
        for k in missing + [None]:
            if k is not None and prev is not None and k == prev + 1:
                prev = k
                continue
            if start is not None and prev - start + 1 >= min_genes:
                runs.append((start, prev, sample))
            start = prev = k
    runs.sort()
    blocks = []
    for a, b, sample in runs:
        if blocks and a <= blocks[-1]['last']:
            blk = blocks[-1]
            blk['last'] = max(blk['last'], b)
            blk['runs'].append([sample, a, b])
        else:
            blocks.append({'first': a, 'last': b, 'runs': [[sample, a, b]]})
    out = []
    for blk in blocks:
        first, last = placed[blk['first']], placed[blk['last']]
        out.append({'first': first, 'last': last,
                    'runs': [[s, placed[a], placed[b]] for s, a, b in blk['runs']]})
    out.sort(key=lambda b: (-len({r[0] for r in b['runs']}), -(b['last'] - b['first'])))
    return out


def read_reference(path):
    """The reference genes as columns, and any warnings about the proteome."""
    genes, warnings = read_proteome(path)
    ref = {
        'names': [g['name'] for g in genes],
        'aa': [g['aa'] for g in genes],
        'placed': [g['start'] is not None for g in genes],
        'len': [(g['end'] - g['start'] + 1) if g['start'] is not None else None for g in genes],
        # The encoder needs a position for every gene; one with none gets 0.
        'start': [g['start'] or 0 for g in genes],
        'end': [g['end'] or 0 for g in genes],
    }
    ref['index'] = {n: i for i, n in enumerate(ref['names'])}
    return ref, warnings


def strand_votes(amap_files, ref):
    """The reference strand of each gene, by vote over every genome: where a
    genome runs collinear with the reference its genes keep the reference's
    strand, where it runs reversed they flip. A first pass over the tables,
    so that no genome's genes stay in memory while the others are read."""
    votes = [0] * len(ref['names'])
    for path in amap_files.values():
        last = {}
        for g in read_amap_table(path):
            if len(g['rvs']) != 1 or g['rvs'][0] not in ref['index']:
                continue
            ri = ref['index'][g['rvs'][0]]
            reverse = g['contig'] in last and ri < last[g['contig']]
            votes[ri] += 1 if (g['strand'] == '+') != reverse else -1
            last[g['contig']] = ri
    return [v >= 0 for v in votes]


def build(args):
    ref, warnings = read_reference(args.ref_proteome)
    n_ref = len(ref['names'])

    amap_files = _files(args.amap, '.AMAP.tsv')
    summaries = _files(args.cases, '-summary.tsv')
    case_files = {c: _files(args.cases, f'-case{c}.tsv') for c in '1234'}
    rvmissing_files = _files(args.cases, '-rvmissing.tsv')
    genome_files = {}
    for suffix in ('.fasta', '.fna', '.fa', '.fas'):
        for k, v in _files(args.genomes, suffix).items():
            genome_files.setdefault(k, v)
    bakta_files = _files(args.bakta, '.tsv')
    bakta_stats = _files(args.bakta, '.txt')
    bakta_gff = _files(args.bakta, '.gff3')
    master_files = _files(args.master, '-AMAPmastertable.tsv')
    myco = read_mycobrowser(args.mycobrowser) if args.mycobrowser and os.path.isfile(args.mycobrowser) else {}

    ref_strand = strand_votes(amap_files, ref)
    # The calls AMAP did not keep, one by one, for runs up to --call-detail
    # genomes; past that they are only counted, and CDS_CASES/ has them.
    call_detail = len(amap_files) <= args.call_detail

    names = sorted(set(amap_files) | set(summaries))
    notes, products = Strings(), Strings()
    samples, status_by_sample, tables = [], {}, {}
    # What the genomes' own annotation calls each reference gene, for its card
    said = defaultdict(lambda: {'gene_name': Counter(), 'product': Counter(), 'cog': Counter()})
    noref_groups = defaultdict(Counter)
    cog, kegg = {}, {}

    for name in names:
        summ = read_summary(summaries[name]) if name in summaries else {'counts': {}}
        lineage = summ.get('lineage') or 'NA'
        s = {'name': name, 'id': summ.get('sample') or name, 'lineage': lineage,
             'cases': summ.get('counts', {})}
        samples.append(s)
        if name not in amap_files:
            s['failed'] = True
            continue
        genes = read_amap_table(amap_files[name])
        bakta = read_bakta_products(bakta_files[name]) if name in bakta_files else {}
        master = read_master(master_files[name]) if name in master_files else {}
        for g in genes:
            m = master.get(g['id'], {})
            b = bakta.get(g['id']) or ('', '')
            g['product'] = m.get('Product_x') or m.get('Product') or m.get('Description') or b[1]
            g['gene_name'] = m.get('Preferred_name') or m.get('Gene') or b[0]
            g['cog'] = m.get('COG_category', '')
            g['kegg'] = m.get('KEGG_Pathway', '')

        # Genome: from its FASTA if given, otherwise from what Bakta says of it
        contig_list = []
        if name in genome_files:
            contig_list, gc = read_fasta_stats(genome_files[name])
            lengths = [l for _, l in contig_list]
            s.update({'size': sum(lengths), 'contigs': len(lengths), 'n50': n50(lengths), 'gc': gc})
        else:
            if name in bakta_stats:
                s.update(read_bakta_stats(bakta_stats[name]))
            if name in bakta_gff:
                contig_list = read_gff_regions(bakta_gff[name])
        contigs = [c for c, _ in contig_list] or list(dict.fromkeys(g['contig'] for g in genes))
        contig_len = [l for _, l in contig_list]
        lengths = dict(contig_list)
        for g in genes:
            g['contig_length'] = lengths.get(g['contig'])

        # Cases
        cases = {c: read_case_table(case_files[c][name]) for c in '1234' if name in case_files[c]}
        rvmissing = read_rvmissing(rvmissing_files[name]) if name in rvmissing_files else set()
        case_of = case_of_calls(cases)
        origins = origins_of(name, genes, case_of, rvmissing)
        final_ids = {g['id'] for g in genes}

        # Every reference gene's state in this genome
        by_ref = defaultdict(list)
        for g in genes:
            for r in g['rvs']:
                if r in ref['index']:
                    by_ref[ref['index'][r]].append(g)
        status, counts = {}, Counter()
        for i in range(n_ref):
            code, detail, best = gene_status(by_ref.get(i, []), ref['len'][i])
            if best is not None:
                counts['found'] += 1
                if len(best['rvs']) == 1:
                    for key in ('gene_name', 'product', 'cog'):
                        if best.get(key):
                            said[i][key][best[key]] += 1
            if code == PSEUDO:
                detail = str(notes(best['note']))
            if code:
                status[i] = (code, detail)
                counts[code] += 1
        status_by_sample[name] = status

        noref = [g for g in genes if not g['rvs']]
        for g in noref:
            noref_groups[g.get('product') or 'no product'][name] += 1
        s.update({
            'genes': len(genes),
            'bakta': s['cases'].get('total_cds'),
            'found': counts['found'],
            'found_pct': round(100.0 * counts['found'] / n_ref, 2) if n_ref else None,
            'missing': counts[MISSING],
            'pseudo': sum(1 for g in genes if g['pseudo']),
            'partial': sum(1 for g in genes if g['partial']),
            'copies': counts[COPIES],
            'lendiff': counts[LENGTH],
            'shared': counts[SHARED],
            'noref': len(noref),
        })

        # Flow: Bakta's calls through the cases to the final genes
        if case_of:
            flow = {}
            for c in '1234':
                tags = {tag for tag, case in case_of.items() if case == c}
                flow['c' + c] = len(tags)
                flow['k' + c] = len(tags & final_ids)
            flow['R'] = origins.count('R')
            flow['A'] = origins.count('A')
            flow['other'] = origins.count('?')
            flow['twice'] = sum(len(v) for v in cases.values()) - len(case_of)
            s['flow'] = flow

        # Functional profile: genes per COG category and KEGG pathway
        if name in master_files:
            c, k = Counter(), Counter()
            for g in genes:
                for x in [x for x in g['cog'] if x.isalpha()] or ['-']:
                    c[x] += 1
                for pw in g['kegg'].split(','):
                    if pw.strip().startswith('map'):
                        k[pw.strip()] += 1
            cog[name], kegg[name] = dict(c), dict(k)

        enc = encode_genes(name, genes, contigs, ref, ref_strand, origins, notes, products)
        enc['cn'] = contigs
        enc['cl'] = contig_len
        if call_detail and case_of:
            enc['gone'] = encode_gone(name, genes, cases, case_of, final_ids)
        tables[name] = enc

    samples.sort(key=lambda s: (lineage_key(s['lineage']), s['id'], s['name']))
    robust_flags([s for s in samples if not s.get('failed')])

    # The reference genes: names, products and categories, from Mycobrowser
    # first, then what eggNOG and Bakta said of the genes naming them.
    ann = {'name': [], 'product': [], 'cat': [], 'refpseudo': [], 'drugs': [], 'cog': [], 'src': []}
    for i, gname in enumerate(ref['names']):
        m = myco.get(gname, {})
        common = {k: (v.most_common(1) or [('', 0)])[0][0] for k, v in said[i].items()} if i in said else {}
        name_ = m.get('name') or common.get('gene_name', '')
        product = m.get('product') or common.get('product', '')
        ann['name'].append(name_ if name_ != gname else '')
        ann['product'].append(product)
        ann['cat'].append(m.get('category', ''))
        ann['refpseudo'].append({'Yes': 1, 'No': 0}.get(m.get('pseudo', ''), -1))
        ann['drugs'].append(m.get('drugs', ''))
        ann['cog'].append(common.get('cog', ''))
        ann['src'].append('m' if m.get('product') else ('e' if product else ''))

    # Status lists, as compact strings: "<index delta><state><detail>" per gene
    status_enc = {}
    for name, status in status_by_sample.items():
        out, prev = [], 0
        for i in sorted(status):
            code, detail = status[i]
            out.append(f'{i - prev}{code}{detail}')
            prev = i
        status_enc[name] = ','.join(out)

    placed_start = [s if p else None for s, p in zip(ref['start'], ref['placed'])]
    blocks = deletion_blocks({'start': placed_start}, status_by_sample)

    functional = {}
    if master_files:
        functional['cog'] = cog
        functional['kegg'] = kegg
    if myco:
        functional['myco_cats'] = sorted({v['category'] for v in myco.values() if v['category']})

    data = OrderedDict()
    data['meta'] = meta(args, warnings, samples, n_ref, ref, myco, bool(master_files))
    data['meta']['call_detail'] = call_detail
    data['ref'] = {'names': ref['names'], 'start': placed_start,
                   'end': [e if p else None for e, p in zip(ref['end'], ref['placed'])],
                   'aa': ref['aa'], 'strand': ''.join('+' if x else '-' for x in ref_strand),
                   **ann}
    data['samples'] = samples
    data['status'] = status_enc
    data['notes'] = notes.items
    data['tables'] = tables
    data['blocks'] = blocks
    data['noref'] = sorted(([p, dict(c)] for p, c in noref_groups.items()),
                           key=lambda x: (-sum(x[1].values()), x[0]))
    data['products'] = products.items
    data['functional'] = functional
    if args.pangene and os.path.isfile(args.pangene) and os.path.getsize(args.pangene) > 0:
        data['pangene'] = pangene_summary(read_rtab(args.pangene), samples)
    if args.communities and os.path.isfile(args.communities) and os.path.getsize(args.communities) > 1:
        data['communities'] = read_communities(args.communities)
    return data


def pangene_summary(rtab, samples):
    """What PANGENE found, per genome: how many of its genes each other
    genome carries. The merged matrix has one row per gene of every genome
    and one column per genome a genome's genes were aligned to, named
    <that genome>_<this genome>#0, so each genome's rows are read against
    its own columns only."""
    if not rtab:
        return None
    names = [s['name'] for s in samples if not s.get('failed')]
    out = {}
    for name in names:
        cols = [i for i, c in enumerate(rtab['columns']) if c.split('#')[0].endswith('_' + name)]
        if not cols:
            continue
        hist = Counter()
        n_genes = 0
        for gene, values in rtab['rows']:
            if not gene.startswith(name + '_'):
                continue
            n_genes += 1
            hist[sum(values[i] for i in cols if i < len(values))] += 1
        if n_genes:
            out[name] = {'others': len(cols), 'hist': {str(k): v for k, v in sorted(hist.items())}}
    return out or None


def meta(args, warnings, samples, n_ref, ref, myco, has_master):
    params, run = (read_params(args.params) if args.params and os.path.isfile(args.params)
                   else (OrderedDict(), OrderedDict()))
    versions = read_versions(args.versions) if args.versions and os.path.isfile(args.versions) else []
    placed = [e for e, p in zip(ref['end'], ref['placed']) if p]
    expected = []
    if args.expected and os.path.isfile(args.expected):
        with open(args.expected) as fh:
            expected = [l.strip() for l in fh if l.strip()]
    seen = {s['id'] for s in samples} | {s['name'] for s in samples}
    return {
        'title': args.title,
        'version': args.amap_version or params.get('version') or '',
        'generated': datetime.now().strftime('%Y-%m-%d %H:%M'),
        # Which reference, as the run's parameters say; without them it is
        # not named rather than guessed
        'reference': {
            'genes': n_ref,
            'length': max(placed) if placed else 0,
            'custom': (params.get('ref_annotation') not in (None, '', 'null')) if params else None,
            'annotation': params.get('ref_annotation', ''),
            'genome': params.get('ref_genome', ''),
        },
        'modes': {
            'functional': has_master,
            'mycobrowser': bool(myco),
            'cases': any('flow' in s for s in samples),
            'genomes': any(s.get('size') for s in samples),
            'typing': any(s['lineage'] not in NO_LINEAGE for s in samples),
            'skip_functional': params.get('skip_functional_annotation') == 'true',
            'skip_typing': params.get('skip_typing') == 'true',
            'skip_mycobrowser': params.get('skip_mycobrowser') == 'true',
        },
        'params': params,
        'run': run,
        'versions': versions,
        'expected_missing': [e for e in expected if e not in seen],
        'warnings': warnings,
    }


# --------------------------------------------------------------------------
# Writing
# --------------------------------------------------------------------------

def render(data, plain=False):
    """The page: the template with the styles, the script and the data in."""
    with open(TEMPLATE) as fh:
        page = fh.read()
    with open(STYLES) as fh:
        styles = fh.read()
    with open(SCRIPT) as fh:
        script = fh.read()
    text = json.dumps(data, separators=(',', ':'), ensure_ascii=False)
    if plain:
        payload, encoding = text.replace('</', '<\\/').replace('\u2028', '\\u2028').replace('\u2029', '\\u2029'), 'json'
    else:
        payload = base64.b64encode(zlib.compress(text.encode('utf-8'), 9)).decode('ascii')
        encoding = 'deflate-base64'
    title = data['meta'].get('title') or 'AMAP report'
    page = page.replace('{{TITLE}}', _html_escape(title))
    with open(LOGO, 'rb') as fh:
        page = page.replace('{{LOGO}}', 'data:image/png;base64,' + base64.b64encode(fh.read()).decode('ascii'))
    page = page.replace('/*{{STYLES}}*/', styles)
    page = page.replace('{{ENCODING}}', encoding)
    # Last: the data may contain any text, including the markers above.
    page = page.replace('/*{{SCRIPT}}*/', script.replace('</script', '<\\/script'))
    page = page.replace('{{DATA}}', payload)
    return page


def _html_escape(text):
    return (text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
            .replace('"', '&quot;'))


def load_embedded(page):
    """The data a written page carries, for the tests."""
    enc = re.search(r'<script id="amap-data" type="application/json" data-encoding="([^"]+)">(.*?)</script>',
                    page, re.S)
    if not enc:
        raise ValueError('no embedded data')
    if enc.group(1) == 'json':
        return json.loads(enc.group(2))
    return json.loads(zlib.decompress(base64.b64decode(enc.group(2))).decode('utf-8'))


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--amap', required=True, help='Folder with every <name>.AMAP.tsv.')
    p.add_argument('--ref-proteome', required=True, help='Reference proteome (headers NAME_START-END).')
    p.add_argument('--cases', help='Folder with the TYPE_CASES tables.')
    p.add_argument('--genomes', help='Folder with the genomes, <name>.fasta.')
    p.add_argument('--bakta', help="Folder with Bakta's outputs: <name>.tsv, .txt and .gff3, "
                                   "directly or in a folder per sample.")
    p.add_argument('--master', help='Folder with the master tables.')
    p.add_argument('--mycobrowser', help='Mycobrowser table of the reference (MTBC only).')
    p.add_argument('--pangene', help="Merged pangene matrix (REPORT_INPUTS' input_pangene.txt).")
    p.add_argument('--communities', help='homologous_communities.tsv from GET_ALLDUP.')
    p.add_argument('--versions', help='software_versions.tsv from SOFTWARE_VERSIONS.')
    p.add_argument('--params', help='params.yml from DUMP_PARAMS.')
    p.add_argument('--expected', help='The sample IDs of the sample sheet, one per line.')
    p.add_argument('--call-detail', type=int, default=100,
                   help='Embed every Bakta call AMAP did not keep for runs of up to this many '
                        'genomes (default 100); larger runs get the counts only.')
    p.add_argument('--amap-version', default='', help='AMAP version, for the header.')
    p.add_argument('--title', default='AMAP report', help='Title of the page.')
    p.add_argument('--plain', action='store_true', help='Embed the data as plain JSON (debugging).')
    p.add_argument('-o', '--out', default='amap_report.html', help='Output file.')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    data = build(args)
    page = render(data, plain=args.plain)
    with open(args.out, 'w', encoding='utf-8') as fh:
        fh.write(page)
    n = sum(1 for s in data['samples'] if not s.get('failed'))
    print(f"amap_report: {n} genome(s), {data['meta']['reference']['genes']} reference genes, "
          f"{len(page) / 1e6:.2f} MB -> {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
