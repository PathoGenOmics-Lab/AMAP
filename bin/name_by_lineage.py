#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Name a typed genome after its sample and its lineage: <lineage>_<sampleId>.

pathotypr types every record of a FASTA on its own and writes one summary row
per record that carries lineage markers, named after that record. A closed
genome is one record and one row. A draft is one row per contig with markers,
in no particular order, and taking the name from the first row named the
sample after one of its contigs (CDC1551d came out as L4_contig_14). So the
name is always built from the sample ID, and the lineage of a draft is the
one pathotypr would call over all its rows together.

    name_by_lineage.py SAMPLE_ID SUMMARY_TSV FASTA

Prints the name and writes the genome to <name>.fasta. A genome of one record
has that record renamed after the sample, as before; the contigs of a draft
keep their own names, so every gene can still be traced back to its contig.
Without a lineage the name is NA_<sampleId>, as before.
"""

import sys
from collections import Counter


def read_rows(summary):
    """(marker counts, pathotypr's call) for each row of a pathotypr summary.

    The counts column reads `L4:5 L4;L4.1:3`: every lineage path with the
    number of markers found for it. The file AMAP writes when pathotypr fails
    has a header and no rows.
    """
    rows = []
    try:
        with open(summary) as fh:
            lines = fh.read().splitlines()[1:]
    except FileNotFoundError:
        return rows
    for line in lines:
        fields = line.split('\t')
        if len(fields) < 3:
            continue
        counts = Counter()
        for item in fields[1].split():
            path, _, n = item.rpartition(':')
            if path and n.isdigit():
                counts[path] += int(n)
        rows.append((counts, fields[-1].strip()))
    return rows


def sample_lineage(rows):
    """The sample's lineage, or None when pathotypr found none.

    One row is a whole genome and keeps the call pathotypr made for it. Over
    several rows, the markers of all contigs are added up and the lineage is
    the call pathotypr makes for one record without --nested-classification:
    the path with the most markers, the first alphabetically on a tie, which
    is always a top-level lineage since its count includes all its sublineages.
    """
    if len(rows) == 1:
        call = rows[0][1]
    else:
        total = Counter()
        for counts, _ in rows:
            total.update(counts)
        call = min(total, key=lambda path: (-total[path], path)) if total else ''
    return call if call and call != 'Unclassified' else None


def write_fasta(src, dst, name):
    """Write the genome to dst, renaming its record after the sample if it has
    only one. Line by line, ending every line with a newline, the way the awk
    this replaces wrote it."""
    with open(src, 'rb') as fh:
        lines = fh.read().split(b'\n')
    if lines and lines[-1] == b'':
        lines.pop()
    headers = [i for i, line in enumerate(lines) if line.startswith(b'>')]
    if len(headers) == 1:
        lines[headers[0]] = b'>' + name.encode()
    with open(dst, 'wb') as out:
        out.writelines(line + b'\n' for line in lines)


def main(argv):
    if len(argv) != 4:
        print(f"Usage: {argv[0]} SAMPLE_ID SUMMARY_TSV FASTA", file=sys.stderr)
        return 2
    sample_id, summary, fasta = argv[1:]
    lineage = sample_lineage(read_rows(summary))
    if lineage is None:
        print(f"No lineage for {sample_id} in {summary}; naming it NA_{sample_id}.", file=sys.stderr)
        lineage = 'NA'
    name = f"{lineage}_{sample_id}"
    write_fasta(fasta, f"{name}.fasta", name)
    print(name)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
