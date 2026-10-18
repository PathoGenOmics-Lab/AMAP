#!/usr/bin/env python3
"""Truth, scoring and comparisons for tests/benchmark/benchmark.sh.

A truth table lists the genes whose right annotation is known, one per row:
name, start, end, strand, status. Only rows with status "ok" are scored.

  truth-genbank   genes of a reference from its own proteome and GenBank genome
  truth-liftover  those genes carried onto another assembly by a whole-genome
                  alignment, kept only where their sequence is identical
  truth-curated   add the curated H37Rv entries written in CP110619 coordinates
  draft           cut a complete genome into contigs, with a map back
  score           how many truth genes an AMAP table has exactly, with the name
  compare         the genes of one run that another run does not have

Everything here is independent of the algorithm being measured: the truth comes
from the reference annotation and from alignments of whole genomes.
"""

import argparse
import os
import random
import re
import subprocess
import sys
import tempfile
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# Curated H37Rv entries named after their accession, which score as the gene.
ALIAS = {"Rv1003549new": "Rv3343c", "Rv1003720new": "Rv3508"}


def _first_sequence(path, fmt="fasta"):
    return str(next(SeqIO.parse(path, fmt)).seq).upper()


def _translate(seq, start, end, strand):
    seg = Seq(seq[start - 1:end])
    if strand == "-":
        seg = seg.reverse_complement()
    return str(seg[:len(seg) - len(seg) % 3].translate(table=11))


def _write(rows, out, columns):
    pd.DataFrame(rows, columns=columns).to_csv(out, sep="\t", index=False)


TRUTH_COLUMNS = ["name", "start", "end", "strand", "status", "len_mod3", "prot_id"]


# --- truth from a reference annotation ----------------------------------------

def truth_genbank(args):
    """Every protein NAME_START-END of the proteome belongs at START..END. The
    strand is not in the header, so both are translated and the one that gives
    the protein back is kept."""
    genome = _first_sequence(args.genbank, "genbank")
    rows, stats = [], defaultdict(int)
    for rec in SeqIO.parse(args.proteome, "fasta"):
        m = re.match(r"^(.+)_(\d+)-(\d+)$", rec.id)
        if not m:
            continue
        name, s, e = m.group(1), int(m.group(2)), int(m.group(3))
        aa = str(rec.seq).rstrip("*")
        ok = [st for st in "+-" if _translate(genome, s, e, st).rstrip("*")[1:] == aa[1:]
              and len(_translate(genome, s, e, st).rstrip("*")) == len(aa)]
        guess = "-" if name.endswith("c") else "+"
        if len(ok) == 1:
            strand, status = ok[0], "ok"
        else:
            strand, status = guess, ("ambiguous" if ok else "untranslatable")
        stats[status] += 1
        rows.append((name, s, e, strand, status, (e - s + 1) % 3, rec.id))
    _write(rows, args.out, TRUTH_COLUMNS)
    print(f"truth from {os.path.basename(args.genbank)}: " + ", ".join(f"{k} {v}" for k, v in sorted(stats.items())))


def truth_curated(args):
    """The curated H37Rv entries (names ending in new, short or long) carry
    CP110619 coordinates, so they join the CP110619 truth as they are."""
    genome = _first_sequence(args.fasta)
    t = pd.read_csv(args.truth, sep="\t")
    t = t[t.status != "untranslatable"]
    add = []
    for rec in SeqIO.parse(args.proteome, "fasta"):
        m = re.match(r"^(.+)_(\d+)-(\d+)$", rec.id)
        if not m or not re.search(r"(new|short|long)$", m.group(1)):
            continue
        name, s, e = m.group(1), int(m.group(2)), int(m.group(3))
        aa = str(rec.seq).rstrip("*")
        ok = [st for st in "+-" if _translate(genome, s, e, st).rstrip("*")[1:] == aa[1:]]
        add.append(dict(name=name, start=s, end=e, strand=ok[0] if ok else "?",
                        status="ok" if len(ok) == 1 else "untranslatable",
                        len_mod3=(e - s + 1) % 3, prot_id=rec.id))
    pd.concat([t, pd.DataFrame(add)], ignore_index=True).to_csv(args.out, sep="\t", index=False)
    print(f"curated entries added: {len(add)}, {sum(a['status'] == 'ok' for a in add)} translate cleanly")


# --- truth carried onto another assembly ----------------------------------------

def _position_map(ref_seq, target_fasta, minimap2, threads):
    """For each 0-based reference position, the 0-based target position it
    aligns to (or -1), from the primary forward alignments of minimap2 asm5."""
    with tempfile.TemporaryDirectory() as tmp:
        ref_fa = os.path.join(tmp, "ref.fa")
        target_fa = os.path.join(tmp, "target.fa")
        with open(ref_fa, "w") as fh:
            fh.write(">ref\n" + ref_seq + "\n")
        with open(target_fa, "w") as fh:
            fh.write(">target\n" + _first_sequence(target_fasta) + "\n")
        paf = subprocess.run([minimap2, "-x", "asm5", "-c", "-t", str(threads), target_fa, ref_fa],
                             check=True, capture_output=True, text=True).stdout
    pos = [-1] * len(ref_seq)
    hits = 0
    for line in paf.splitlines():
        f = line.split("\t")
        tags = dict(t.split(":", 2)[0::2] for t in f[12:])
        if f[4] != "+" or tags.get("tp") != "P" or "cg" not in tags:
            continue
        hits += 1
        q, r = int(f[2]), int(f[7])
        for length, op in re.findall(r"(\d+)([MIDNSHP=X])", tags["cg"]):
            length = int(length)
            if op in "M=X":
                pos[q:q + length] = range(r, r + length)
                q += length
                r += length
            elif op == "I":     # in the reference only
                q += length
            elif op == "D":     # in the target only
                r += length
    return pos, hits


def truth_liftover(args):
    """A gene's coordinates in the target are where its reference coordinates
    land in the alignment, and it counts as truth only if its sequence there is
    identical, same length included."""
    fmt = "genbank" if args.reference.endswith((".gb", ".gbk", ".gbff")) else "fasta"
    ref = _first_sequence(args.reference, fmt)
    target = _first_sequence(args.fasta)
    pos, hits = _position_map(ref, args.fasta, args.minimap2, args.threads)
    t = pd.read_csv(args.truth, sep="\t")
    rows, counts = [], defaultdict(int)
    for g in t.to_dict("records"):
        if g["status"] != "ok":
            counts["excluded"] += 1
            rows.append(g)
            continue
        s, e = pos[g["start"] - 1], pos[g["end"] - 1]
        if s < 0 or e < 0:
            status = "unmapped"
        else:
            s, e = s + 1, e + 1
            same = (e - s) == (g["end"] - g["start"]) and target[s - 1:e] == ref[g["start"] - 1:g["end"]]
            status = "ok" if same else "changed"
            g = dict(g, start=s, end=e)
        counts[status] += 1
        rows.append(dict(g, status=status))
    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    print(f"lifted onto {os.path.basename(args.fasta)} ({hits} alignments): "
          + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))


# --- a draft assembly -----------------------------------------------------------

def draft(args):
    """Cut a complete genome at random points into contigs of very different
    sizes, about half of them reverse-complemented, with a fixed seed so every
    run cuts the same way."""
    seq = str(next(SeqIO.parse(args.fasta, "fasta")).seq)
    rng = random.Random(20261017)
    cuts = set()
    while len(cuts) < args.n - 1:
        c = rng.randrange(5_000, len(seq) - 5_000)
        if all(abs(c - x) >= 3_000 for x in cuts):
            cuts.add(c)
    bounds = [0] + sorted(cuts) + [len(seq)]
    with open(args.out_fasta, "w") as fa, open(args.out_map, "w") as mp:
        mp.write("contig\torig_start\torig_end\treverse\tlength\n")
        for i in range(args.n):
            s, e = bounds[i], bounds[i + 1]
            piece = seq[s:e]
            rev = rng.random() < 0.5
            if rev:
                piece = str(Seq(piece).reverse_complement())
            name = args.sample if i == 0 else f"contig_{i + 1}"
            fa.write(f">{name}\n" + "\n".join(piece[j:j + 80] for j in range(0, len(piece), 80)) + "\n")
            mp.write(f"{name}\t{s + 1}\t{e}\t{int(rev)}\t{e - s}\n")
    print(f"{args.n} contigs from {os.path.basename(args.fasta)}")


def _truth_on_contigs(truth, contig_map):
    """The truth of a complete genome in contig coordinates; genes that span a
    cut are left out."""
    cmap = pd.read_csv(contig_map, sep="\t")
    rows = []
    for t in truth.itertuples():
        hit = cmap[(cmap.orig_start <= t.start) & (cmap.orig_end >= t.end)]
        if hit.empty:
            continue
        c = hit.iloc[0]
        if c.reverse:
            s, e = c.orig_end - t.end + 1, c.orig_end - t.start + 1
            st = "-" if t.strand == "+" else "+"
        else:
            s, e, st = t.start - c.orig_start + 1, t.end - c.orig_start + 1, t.strand
        rows.append(dict(name=t.name, contig=c.contig, start=s, end=e, strand=st))
    return pd.DataFrame(rows, columns=["name", "contig", "start", "end", "strand"])


# --- scoring ----------------------------------------------------------------------

def _base(name):
    name = ALIAS.get(name, name)
    return re.sub(r"\d?new$", "", name)


def _load_amap(path):
    a = pd.read_csv(path, sep="\t", dtype={"rvs": str, "contig": str}).fillna({"rvs": ""})
    if "contig" not in a.columns:
        a["contig"] = ""
    return a


def _frame(start, end, strand):
    return (end if strand == "-" else start) % 3


def _doubles(df, contig_col):
    """Pairs of genes on one ORF: same contig and strand, same reading frame,
    overlapping by at least half of the shorter."""
    df = df.assign(_c=df[contig_col] if contig_col else "")
    n = 0
    for _, g in df.groupby(["_c", "strand"]):
        rows = list(g.sort_values("start")[["start", "end", "strand"]].itertuples(index=False))
        for i, a in enumerate(rows):
            for b in rows[i + 1:]:
                if b.start > a.end:
                    break
                ov = min(a.end, b.end) - max(a.start, b.start) + 1
                if (ov >= 0.5 * min(a.end - a.start + 1, b.end - b.start + 1)
                        and _frame(a.start, a.end, a.strand) == _frame(b.start, b.end, b.strand)):
                    n += 1
    return n


def score(args):
    """Truth genes found at exactly their coordinates and strand, labelled with
    their own name; plus the ORFs the table annotates twice."""
    truth = pd.read_csv(args.truth, sep="\t")
    truth = truth[truth.status == "ok"]
    amap = _load_amap(args.amap)
    on_contigs = args.contig_map is not None
    if on_contigs:
        truth = _truth_on_contigs(truth, args.contig_map)
        key = lambda r: (r.contig, int(r.start), int(r.end), r.strand)
    else:
        key = lambda r: (int(r.start), int(r.end), r.strand)
    names = defaultdict(set)
    for r in amap.itertuples():
        names[key(r)] |= {_base(x.strip()) for x in r.rvs.split(",") if x.strip()}
    res = defaultdict(list)
    for t in truth.itertuples():
        k = key(t)
        if k not in names:
            res["not exact"].append(t.name)
        elif _base(t.name) not in names[k]:
            res["wrong label"].append(t.name)
        else:
            res["exact"].append(t.name)
    n = len(truth)
    doubles = _doubles(amap, "contig")
    line = {"genome": args.label, "truth": n, "exact": len(res["exact"]),
            "wrong_label": len(res["wrong label"]), "not_exact": len(res["not exact"]),
            "genes": len(amap), "pseudogenes": int((amap.pseudogene == "yes").sum()),
            "orfs_twice": doubles}
    if on_contigs and args.fasta and args.faa:
        line.update(_own_contig(amap, args.fasta, args.faa))
    print(f"{args.label:<14} exact {len(res['exact'])}/{n} ({100 * len(res['exact']) / max(n, 1):.2f}%)"
          f"  wrong label {len(res['wrong label'])}  not exact {len(res['not exact'])}"
          f"  | genes {len(amap)}  ORFs twice {doubles}"
          + (f"  | proteins from their own contig {line['own_contig']}, from another {line['other_contig']}"
             if "own_contig" in line else ""))
    for k in ("wrong label", "not exact"):
        if res[k]:
            print(f"    {k}: {', '.join(res[k][:12])}{' ...' if len(res[k]) > 12 else ''}")
    if args.tsv:
        new = not os.path.exists(args.tsv)
        pd.DataFrame([line]).to_csv(args.tsv, sep="\t", index=False, mode="a", header=new)


def _own_contig(amap, fasta, faa):
    """Whether each protein in the .faa is the translation of its gene on the
    contig the table puts it on."""
    contigs = {r.id: str(r.seq).upper() for r in SeqIO.parse(fasta, "fasta")}
    proteins = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(faa, "fasta")}
    own = other = 0
    for r in amap.itertuples():
        rec = proteins.get(f"{r.id};{r.strand};{int(r.start) - 1}-{int(r.end)}")
        if rec is None or r.contig not in contigs:
            continue
        if _translate(contigs[r.contig], int(r.start), int(r.end), r.strand).rstrip("*") == rec:
            own += 1
        else:
            other += 1
    return {"own_contig": own, "other_contig": other}


def compare(args):
    """Genes of the first run missing from the second, split by whether a gene
    of the second sits on the same reading frame there, and the reverse."""
    a, b = _load_amap(args.before), _load_amap(args.after)
    for name, x, y in (("only before", a, b), ("only now", b, a)):
        ky = set(zip(y.contig, y.start, y.end, y.strand))
        lost, replaced = [], 0
        for r in x.itertuples():
            if (r.contig, r.start, r.end, r.strand) in ky:
                continue
            o = y[(y.contig == r.contig) & (y.strand == r.strand) & (y.start <= r.end) & (y.end >= r.start)]
            same = [s for s in o.itertuples() if _frame(s.start, s.end, s.strand) == _frame(r.start, r.end, r.strand)]
            if same:
                replaced += 1
            else:
                lost.append(f"{r.start}-{r.end}{r.strand}[{r.rvs or r.id}]")
        print(f"  {name}: {replaced + len(lost)} ({replaced} with a gene in the same frame there,"
              f" {len(lost)} with none{': ' + ', '.join(lost[:8]) if lost else ''})")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("truth-genbank")
    s.add_argument("genbank")
    s.add_argument("proteome")
    s.add_argument("out")
    s.set_defaults(func=truth_genbank)

    s = sub.add_parser("truth-liftover")
    s.add_argument("reference", help="reference genome, GenBank or FASTA")
    s.add_argument("fasta", help="genome to carry the truth onto (first sequence)")
    s.add_argument("truth", help="truth of the reference")
    s.add_argument("out")
    s.add_argument("--minimap2", default="minimap2")
    s.add_argument("--threads", type=int, default=4)
    s.set_defaults(func=truth_liftover)

    s = sub.add_parser("truth-curated")
    s.add_argument("truth")
    s.add_argument("fasta")
    s.add_argument("proteome")
    s.add_argument("out")
    s.set_defaults(func=truth_curated)

    s = sub.add_parser("draft")
    s.add_argument("fasta")
    s.add_argument("out_fasta")
    s.add_argument("out_map")
    s.add_argument("sample")
    s.add_argument("n", type=int)
    s.set_defaults(func=draft)

    s = sub.add_parser("score")
    s.add_argument("label")
    s.add_argument("truth")
    s.add_argument("amap")
    s.add_argument("--contig-map", help="score a draft: truth of the complete genome and its contig map")
    s.add_argument("--fasta", help="draft genome, to check each protein against its own contig")
    s.add_argument("--faa", help="the run's .AMAP.faa")
    s.add_argument("--tsv", help="append the scores to this table")
    s.set_defaults(func=score)

    s = sub.add_parser("compare")
    s.add_argument("before")
    s.add_argument("after")
    s.set_defaults(func=compare)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
