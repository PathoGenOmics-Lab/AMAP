#!/usr/bin/env python3
"""Write docs/assets/figure-data.js, the real genes the website draws.

Every drawing on the site is a window of an AMAP run on a genome whose right
answer is known: the seven genomes of tests/benchmark. This reads a benchmark
work directory and a run of it, and writes, for each window, every Bakta call
in it with its case, every placement of a reference protein the
classification saw, Bakta's known false positives there, what the
classification wrote and what AMAP kept, plus what the checks of the harder
steps need (the stretch of genome check_cds.py reads, the anchors of the
paralog step, the input of deduplicate_rvs.py).

tools/check_site.py runs every window back through bin/ and fails when the
code no longer gives what is drawn. When a change to the algorithm does that,
rerun the benchmark and this, then read the words of the drawings again:

    tests/benchmark/benchmark.sh benchmark_work final      # in the container
    python3 tools/figure_data.py benchmark_work final --reference /reference

--reference is where aa_H37Rv.translated.fas and sequence.gb are: /reference
inside the AMAP container. Needs pandas, intervaltree and biopython.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "bin"))
from classify.exact import place_exact_matches  # noqa: E402
from classify.io import read_bakta, read_miniprot, read_miniprot_unfiltered  # noqa: E402
from classify.synteny import drop_paralog_placements  # noqa: E402

OUT = ROOT / "docs" / "assets" / "figure-data.js"
HEADER = """/* The genes the drawings show. Every one is real: a window of AMAP run on a
   genome whose right answer is known (tests/benchmark), with every Bakta call
   in it ("bakta", with the case the run gave it), every placement of a
   reference protein the classification saw ("ref"), Bakta's known false
   positives aligned there ("bad"), how the case 1 calls' proteins match the
   reference genome ("refMatch"), what the classification wrote ("corrected")
   and the genes AMAP kept ("amap"). Where check_cds.py changed a gene, the
   window carries the stretch of genome it read ("translate"); a window with a
   "kind" or "exact" carries what that step's check needs. tools/check_site.py
   runs every window through bin/ and fails if anything comes out different.
   Written by tools/figure_data.py.

   Everything after the "=" is JSON. */
window.AMAP_FIGURE_DATA = """


def gene_names(genbank: Path) -> dict[str, str]:
    names = {}
    for rec in SeqIO.parse(genbank, "genbank"):
        for f in rec.features:
            if f.type == "CDS" and "locus_tag" in f.qualifiers:
                tag = f.qualifiers["locus_tag"][0]
                names[tag] = f.qualifiers.get("gene", [tag])[0]
    return names


class Genome:
    """One genome of the benchmark: its inputs, and a run of the algorithm on it."""

    def __init__(self, label: str, inputs: Path, sample: str, run: Path, proteome: Path, names: dict[str, str]):
        self.label, self.sample, self.names, self.run, self.proteome = label, sample, names, run, proteome
        self.inputs = inputs
        full = pd.read_csv(inputs / sample / f"{sample}.tsv", comment="#", sep="\t", header=None)
        full.columns = ["Sequence Id", "Type", "Start", "End", "Strand", "Locus Tag", "Gene", "Product", "DbXrefs"]
        self.gene = {t: (g if isinstance(g, str) else "") for t, g in zip(full["Locus Tag"], full["Gene"])}
        self.product = {t: (p if isinstance(p, str) else "") for t, p in zip(full["Locus Tag"], full["Product"])}
        self.bakta = read_bakta(inputs / sample / f"{sample}.tsv")
        self.raw = read_miniprot(inputs / f"{sample}.paf", ref_coords=True)
        self.exact = place_exact_matches(self.raw, str(inputs / f"{sample}.fasta"), str(proteome))
        self.syn = drop_paralog_placements(self.exact)
        self.case_of = {}
        for c in ("case1", "case2", "case3", "case4"):
            for tag in pd.read_csv(run / f"{sample}-{c}.tsv", sep="\t")["Locus Tag"]:
                self.case_of[tag] = int(c[-1])
        self.bad = read_miniprot_unfiltered(inputs / f"{sample}-badbakta.paf")
        self.refpaf = read_miniprot_unfiltered(inputs / f"{sample}-ref.paf")
        self.amap = pd.read_csv(run / f"{sample}.AMAP.tsv", sep="\t", dtype=str).fillna("")
        self.amap["start"] = self.amap.start.astype(int)
        self.amap["end"] = self.amap.end.astype(int)
        self.corr = pd.read_csv(run / f"{sample}_corrected.tsv", sep="\t", header=None, dtype=str,
                                names=["id", "st", "s", "e", "rvs", "contig"]).fillna("")
        self.corr["s"] = self.corr.s.astype(int)
        self.corr["e"] = self.corr.e.astype(int)
        self.seqs = {r.id: str(r.seq).upper() for r in SeqIO.parse(inputs / f"{sample}.fasta", "fasta")}
        self.proteins, self.spans = {}, {}
        for r in SeqIO.parse(proteome, "fasta"):
            m = re.match(r"^(.+)_(\d+)-(\d+)$", r.id)
            name = m.group(1) if m else r.id
            self.proteins.setdefault(name, str(r.seq))
            self.spans.setdefault(name, int(m.group(3)) - int(m.group(2)) + 1 if m else None)

    def window(self, contig: str, lo: int, hi: int, raw: bool = False, placed: bool = False) -> dict:
        b, syn, am = self.bakta, self.syn, self.amap
        out = {"genome": self.label, "base": self.sample, "contig": contig, "from": lo, "to": hi,
               "bakta": [], "ref": [], "bad": [], "refMatch": [], "amap": []}
        for r in b[(b.Contig == contig) & (b.End >= lo) & (b.Start <= hi)].to_dict("records"):
            tag = r["Locus Tag"]
            out["bakta"].append({"tag": tag, "s": int(r["Start"]), "e": int(r["End"]), "st": r["Strand"],
                                 "name": self.gene.get(tag, ""), "product": self.product.get(tag, ""), "case": self.case_of.get(tag)})
            if self.case_of.get(tag) == 1:
                for x in self.refpaf[self.refpaf["Prot Id"] == tag].to_dict("records"):
                    out["refMatch"].append({"tag": tag, "match": int(x["Match Nt"]), "num": int(x["Num nt"])})
        for r in syn[(syn.Contig == contig) & (syn.End >= lo) & (syn.Start <= hi)].to_dict("records"):
            out["ref"].append({"prot": r["Prot Id"], "name": self.names.get(r["Prot Id"], r["Prot Id"]), "s": int(r["Start"]),
                               "e": int(r["End"]), "st": r["Strand"], "len": int(r["Prot len"]), "match": int(r["Match Nt"]),
                               "num": int(r["Num nt"]), "complete": bool(r["Complete"])})
        if raw:
            out["miniprot"] = [{"prot": r["Prot Id"], "s": int(r["Start"]), "e": int(r["End"]), "st": r["Strand"],
                                "match": int(r["Match Nt"]), "num": int(r["Num nt"])}
                               for r in self.raw[(self.raw.Contig == contig) & (self.raw.End >= lo) & (self.raw.Start <= hi)].to_dict("records")]
        if placed:
            # Before the paralog step: every placement exact.py and miniprot made.
            ex = self.exact
            out["placed"] = [{"prot": r["Prot Id"], "s": int(r["Start"]), "e": int(r["End"]), "st": r["Strand"],
                              "match": int(r["Match Nt"]), "num": int(r["Num nt"])}
                             for r in ex[(ex.Contig == contig) & (ex.End >= lo) & (ex.Start <= hi)].to_dict("records")]
        for r in self.bad[(self.bad.Contig == contig) & (self.bad.End >= lo) & (self.bad.Start <= hi)].to_dict("records"):
            out["bad"].append({"prot": r["Prot Id"], "s": int(r["Start"]), "e": int(r["End"]), "st": r["Strand"]})
        for r in am[(am.contig == contig) & (am.end >= lo) & (am.start <= hi)].to_dict("records"):
            out["amap"].append({"id": r["id"], "rvs": r["rvs"], "s": int(r["start"]), "e": int(r["end"]), "st": r["strand"],
                                "pseudo": r["pseudogene"] == "yes", "partial": r["partial"], "note": r["note"]})
        # What the classification wrote, before check_cds.py translated it.
        cr = self.corr[(self.corr.contig == contig) & (self.corr.e >= lo) & (self.corr.s <= hi)]
        out["corrected"] = [{"id": r["id"], "rvs": r["rvs"], "s": int(r["s"]), "e": int(r["e"]), "st": r["st"]} for r in cr.to_dict("records")]
        changed = (sorted((c["id"], c["s"], c["e"]) for c in out["corrected"]) != sorted((a["id"], a["s"], a["e"]) for a in out["amap"])
                   or any(a["pseudo"] or a["partial"] or a["note"] for a in out["amap"]))
        if changed:
            # The stretch of genome check_cds.py reads, with room either side
            # for the stops it looks for; a contig's own ends stay its ends.
            seq = self.seqs[contig]
            first = max(1, min([c["s"] for c in out["corrected"]] + [lo]) - 1200)
            last = min(len(seq), max([c["e"] for c in out["corrected"]] + [hi]) + 1200)
            prots = {}
            for c in out["corrected"]:
                for name in [x.strip() for x in c["rvs"].split(",") if x.strip()]:
                    if name in self.proteins:
                        prots[name] = [self.proteins[name], self.spans[name]]
            out["translate"] = {"offset": first - 1, "dna": seq[first - 1:last], "contigLength": len(seq), "proteins": prots}
        return out


def paralog_anchors(cdc: Genome, proteins: list[str]) -> dict:
    """Both loci of two paralogs, and the five proteins placed once on each side
    of each, which is all drop_paralog_placements reads to tell them apart."""
    ex = cdc.exact[(cdc.exact.Start < cdc.exact.End) & (cdc.exact.Contig != "*")]
    n_place = ex.groupby("Prot Id")["Prot Id"].transform("size")
    anchors = ex[(n_place == 1) & ex["Ref Start"].notna()].copy()
    anchors["mid"] = (anchors.Start + anchors.End) / 2
    focal = ex[ex["Prot Id"].isin(proteins)]
    rows = [focal]
    for _, f in focal.iterrows():
        m = (f.Start + f.End) / 2
        a = anchors.sort_values("mid")
        rows += [a[a.mid < m].tail(5), a[a.mid >= m].head(5)]
    mini = pd.concat(rows).drop_duplicates(subset=["Prot Id", "Start", "End", "Strand"]).sort_values("Start").drop(columns=["mid"])
    kept = drop_paralog_placements(mini.reset_index(drop=True))
    keys = set(zip(kept["Prot Id"], kept.Start, kept.End))
    return {
        "genome": cdc.label, "base": cdc.sample, "contig": cdc.sample, "kind": "synteny",
        "placements": [{"prot": r["Prot Id"], "name": cdc.names.get(r["Prot Id"], r["Prot Id"]), "s": int(r["Start"]), "e": int(r["End"]),
                        "st": r["Strand"], "refStart": int(r["Ref Start"]), "refEnd": int(r["Ref End"]), "match": int(r["Match Nt"]),
                        "num": int(r["Num nt"]), "len": int(r["Prot len"]), "complete": bool(r["Complete"]),
                        "kept": (r["Prot Id"], r["Start"], r["End"]) in keys}
                       for r in mini.to_dict("records")],
    }


def copies(sakai: Genome, protein: str, name: str, loci: list[list[int]]) -> dict:
    """The copies of a protein, as check_cds.py hands them to deduplicate_rvs.py
    ("before") and as it keeps them ("genes")."""
    with tempfile.TemporaryDirectory() as tmp:
        corrected = Path(tmp) / f"{sakai.sample}_corrected.tsv"
        corrected.write_text((sakai.run / f"{sakai.sample}_corrected.tsv").read_text())
        subprocess.run([sys.executable, str(ROOT / "bin" / "check_cds.py"), str(sakai.inputs / f"{sakai.sample}.fasta"),
                        str(corrected), "--ref-proteome", str(sakai.proteome)], cwd=tmp, check=True, capture_output=True)
        pre = pd.read_csv(Path(tmp) / f"{sakai.sample}.AMAP.tsv", sep="\t", dtype=str).fillna("")
    pre["start"] = pre.start.astype(int)
    pre["end"] = pre.end.astype(int)

    def rows(table):
        out = []
        for lo, hi in loci:
            for r in table[(table.end >= lo) & (table.start <= hi)].to_dict("records"):
                out.append({"id": r["id"], "rvs": r["rvs"], "s": int(r["start"]), "e": int(r["end"]), "st": r["strand"],
                            "pseudo": r["pseudogene"] == "yes"})
        return out

    return {
        "genome": sakai.label, "base": sakai.sample, "contig": sakai.sample, "kind": "copies", "prot": protein, "name": name,
        "before": rows(pre), "refResidues": len(sakai.proteins[protein].rstrip("*")), "loci": loci, "genes": rows(sakai.amap),
        "bakta": [{"tag": r["Locus Tag"], "s": int(r["Start"]), "e": int(r["End"]), "st": r["Strand"]} for lo, hi in loci
                  for r in sakai.bakta[(sakai.bakta.End >= lo) & (sakai.bakta.Start <= hi)].to_dict("records")],
    }


def exact_hit(k12: Genome, protein: str, s: int, e: int) -> dict:
    """What exact.py needs to find a protein again: it, and the bases around it."""
    lo, hi = s - 30, e + 30
    return {
        "prot": protein, "protein": k12.proteins[protein], "offset": lo - 1, "dna": k12.seqs[k12.sample][lo - 1:hi],
        "elsewhere": [{"s": int(r["Start"]), "e": int(r["End"]), "st": r["Strand"], "match": int(r["Match Nt"]), "num": int(r["Num nt"])}
                      for r in k12.raw[k12.raw["Prot Id"] == protein].to_dict("records")],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("work", type=Path, help="the benchmark work directory (inputs/, runs/)")
    parser.add_argument("run", help="the run under work/runs to draw, or a directory holding one folder per genome")
    parser.add_argument("--reference", type=Path, default=Path("/reference"),
                        help="folder with aa_H37Rv.translated.fas and sequence.gb (default: /reference)")
    args = parser.parse_args()
    inputs = args.work / "inputs"
    runs = Path(args.run) if Path(args.run).is_dir() else args.work / "runs" / args.run
    mtb = args.reference / "aa_H37Rv.translated.fas"
    eco = inputs / "ecoli" / "custom_ref_ref_proteins.faa"
    mtb_names = gene_names(args.reference / "sequence.gb")
    eco_names = gene_names(inputs / "ecoli" / "k12.gbff")

    cdc = Genome("M. tuberculosis CDC1551", inputs / "mtb_cdc1551", "L4_CDC1551", runs / "CDC1551", mtb, mtb_names)
    c = cdc.sample
    w = {
        "hero": cdc.window(c, 1, 12400),
        "case2": cdc.window(c, 1980, 3275),
        "case3": cdc.window(c, 23780, 25650),
        "case4-start": cdc.window(c, 1, 1700),
        "case4-split": cdc.window(c, 9850, 11000),
        "case4-longer": cdc.window(c, 2267650, 2272700),
        "case1-kept": cdc.window(c, 1979450, 1985565),
        "case1-dropped": cdc.window(c, 11400, 11545),
        "missing": cdc.window(c, 33100, 34420),
        "paralogs": paralog_anchors(cdc, ["Rv0336", "Rv0515"]),
        "paralogs-1": cdc.window(c, 400170, 402300, placed=True),
        "paralogs-2": cdc.window(c, 607300, 609900, placed=True),
    }
    # Small proteins: K-12 against its own annotation, where miniprot put IbsC
    # (b4665) only on IbsD's locus and the exact search finds its own.
    k12 = Genome("E. coli K-12", inputs / "ecoli", "K12self", runs / "K-12", eco, eco_names)
    w["small-protein"] = k12.window(k12.sample, 3056150, 3058450, raw=True)
    w["small-protein"]["exact"] = exact_hit(k12, "b4665", 3056890, 3056949)
    # A contig end: the CDC1551 draft's contig_10.
    draft = Genome("M. tuberculosis CDC1551 in 40 contigs", inputs / "mtb_draft", "L4_CDC1551d", runs / "CDC1551-draft", mtb, mtb_names)
    w["contig-end"] = draft.window("contig_10", 100300, 102146)
    # A frame that runs into the next gene: M. bovis Rv0590 and Rv0590A.
    bovis = Genome("M. bovis AF2122/97", inputs / "mtb_mbovis", "La1_Mbovis", runs / "M.bovis", mtb, mtb_names)
    w["frameshift"] = bovis.window(bovis.sample, 688900, 690500)
    # Real copies: E. coli O157:H7 Sakai's two rhsD (b0497).
    sakai = Genome("E. coli O157:H7 Sakai", inputs / "ecoli", "Sakai", runs / "Sakai", eco, eco_names)
    w["copies"] = copies(sakai, "b0497", "rhsD", [[616500, 621800], [2045500, 2050800]])
    w["copies-1"] = sakai.window(sakai.sample, 616700, 621600)
    w["copies-2"] = sakai.window(sakai.sample, 2045700, 2050600)

    OUT.write_text(HEADER + json.dumps({"windows": w}, indent=1, ensure_ascii=False) + ";\n", encoding="utf-8")
    print(f"{OUT.relative_to(ROOT)}: {len(w)} windows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
