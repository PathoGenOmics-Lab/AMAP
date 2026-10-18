#!/usr/bin/env python3
"""Score an AMAP annotation of H37Rv against its own GenBank record.

Annotating the reference with itself has one right answer: every CDS in
AL123456.3, at exactly its coordinates and strand, named after the reference
protein it is. A gene counts as named right when its rvs column holds the
GenBank locus tag, or a reference protein that aligns to the GenBank
translation at 90% identity or more: the AMAP proteome carries a few curated
versions, corrected against the CP110619 assembly, under their own names
(Rv1450cnew for Rv1450c, Rv1003549new for Rv3343c).

Exits 1 if any CDS is missing, misplaced or misnamed, so it can gate a change
to the correction algorithm.

    score_selftest.py sequence.gb ref_proteome.faa L4_H37Rv.AMAP.tsv
"""
import re
import sys

import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner


def identity(a, b):
    """Identities in the best global alignment, over the shorter sequence.

    Over the shorter one, because PE/PPE genes differ between assemblies in
    how many repeats they carry: the curated PPE54 has 3099 residues and the
    AL123456.3 one 2523, and they are the same gene.
    """
    aligner = PairwiseAligner()
    aligner.match_score, aligner.mismatch_score = 1, 0
    aligner.open_gap_score = aligner.extend_gap_score = 0
    return aligner.score(a, b) / min(len(a), len(b)) if a and b else 0.0


def main(genbank, proteome, amap_tsv):
    record = next(SeqIO.parse(genbank, "genbank"))
    cds = {}
    for f in record.features:
        if f.type != "CDS":
            continue
        key = (int(f.location.start) + 1, int(f.location.end),
               "+" if f.location.strand == 1 else "-")
        cds[key] = (f.qualifiers.get("locus_tag", ["?"])[0],
                    f.qualifiers.get("translation", [""])[0])

    proteins = {re.sub(r"_\d+-\d+$", "", r.id): str(r.seq).rstrip("*")
                for r in SeqIO.parse(proteome, "fasta")}

    amap = pd.read_csv(amap_tsv, sep="\t", dtype={"rvs": str}).fillna({"rvs": ""})
    found = {}
    for _, r in amap.iterrows():
        found[(int(r.start), int(r.end), r.strand)] = [x.strip() for x in r.rvs.split(",") if x.strip()]

    missing, misnamed = [], []
    for key, (tag, translation) in sorted(cds.items()):
        if key not in found:
            missing.append((tag, key))
            continue
        names = found[key]
        if tag in names:
            continue
        if any(identity(proteins.get(n, ""), translation) >= 0.9 for n in names):
            continue
        misnamed.append((tag, key, names))

    extra = [k for k in found if k not in cds]
    print(f"GenBank CDS: {len(cds)}   AMAP genes: {len(found)}")
    print(f"  exact and named right: {len(cds) - len(missing) - len(misnamed)}")
    print(f"  not at exactly their coordinates: {len(missing)}")
    print(f"  named wrong: {len(misnamed)}")
    print(f"  AMAP genes that are no GenBank CDS: {len(extra)}")
    for tag, key in missing:
        print(f"    missing   {tag} {key[0]}-{key[1]} {key[2]}")
    for tag, key, names in misnamed:
        print(f"    misnamed  {tag} {key[0]}-{key[1]} {key[2]} as {', '.join(names) or '(nothing)'}")
    for key in sorted(extra):
        print(f"    extra     {key[0]}-{key[1]} {key[2]} [{', '.join(found[key])}]")
    return 1 if missing or misnamed or extra else 0


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    sys.exit(main(*sys.argv[1:]))
