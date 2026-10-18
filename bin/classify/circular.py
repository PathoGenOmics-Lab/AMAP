"""
Genes that run across the origin of a circular sequence.

Bakta joins the two ends of a sequence it takes for a complete replicon (one
whose FASTA description says "complete", as NCBI's do, or "circular=true", as
Unicycler's do) and calls a gene that runs across its origin from its start,
near the end, to its stop, near the beginning: stcE on plasmid pO157 of
E. coli O157:H7 Sakai is 92563..2502. AMAP's tables keep such a gene the way
Bakta's does, start after end. Measuring it takes the length of its sequence,
which classify_cases.py and the scripts after it read from the genome and keep
here.

A gene across the origin is handled as two copies of itself on an unrolled
circle, (start, end + L) and (start - L, end): each holds the whole gene, and
on 1..L each shows one of its two pieces, read in the gene's frame there.
Anything that sits on 1..L overlaps a copy exactly where it overlaps the gene.
"""

_lengths = {}


def set_lengths(lengths):
    """Keep the length of every sequence of the genome, by name."""
    _lengths.clear()
    _lengths.update({str(name): int(n) for name, n in lengths.items()})


def read_lengths(fasta_path):
    """The length of every sequence of a FASTA file, by its first word."""
    lengths, name = {}, None
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith('>'):
                name = line[1:].split()[0] if line[1:].strip() else ''
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths


def bakta_circular(header):
    """Whether Bakta takes a sequence for a complete, circular replicon, from
    its FASTA header without the '>': the way Bakta 1.12 does, by the words
    after the first space, "circular=true" (Unicycler) or "complete" without
    "complete=false" (NCBI's "complete genome", "complete sequence"). Those
    are the sequences it joins end to end, and calls genes across."""
    description = header.split(' ', 1)[1].lower() if ' ' in header else ''
    return 'circular=true' in description or \
        ('complete' in description and 'complete=false' not in description)


def wraps(start, end):
    """Whether start..end runs across the origin. A protein miniprot could
    not place has no place at all (1..0 once read), and does not."""
    return start > end > 0


def contig_length(contig):
    try:
        return _lengths[str(contig)]
    except KeyError:
        raise ValueError(f"a gene runs across the origin of {contig}, a sequence whose "
                         f"length is not known: give the genome it comes from") from None


def copies(start, end, contig=''):
    """The gene as stretches of one line: itself, or both ways of unrolling
    the circle for a gene across the origin."""
    if not wraps(start, end):
        return [(start, end)]
    n = contig_length(contig)
    return [(start, end + n), (start - n, end)]


def pieces(start, end, contig=''):
    """The stretches of the sequence the gene covers: two for a gene across
    the origin, its end and its beginning."""
    if not wraps(start, end):
        return [(start, end)]
    return [(start, contig_length(contig)), (1, end)]


def length(start, end, contig=''):
    """Bases of a gene, 1-based inclusive."""
    if not wraps(start, end):
        return end - start + 1
    return sum(e - s + 1 for s, e in pieces(start, end, contig))


def shared(a, b, contig=''):
    """Bases two features of one sequence have in common, either of them
    across its origin."""
    if not wraps(a[0], a[1]) and not wraps(b[0], b[1]):
        return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)
    return sum(max(0, min(e1, e2) - max(s1, s2) + 1)
               for s1, e1 in pieces(a[0], a[1], contig)
               for s2, e2 in pieces(b[0], b[1], contig))


def facing(coord, other, contig=''):
    """The copy of coord that lies over other (or nearest it), so that their
    reading frames can be compared from their 5' ends."""
    options = copies(coord[0], coord[1], contig)
    if len(options) == 1:
        return options[0]
    target = copies(other[0], other[1], contig)[0]
    return max(options, key=lambda c: min(c[1], target[1]) - max(c[0], target[0]))
