#!/usr/bin/env python3
"""Build a demonstration dataset for the AMAP report.

The values are synthetic, but the distributions are anchored on the real
M. tuberculosis CDC1551 run measured earlier: 4020 CDS, 65.61% GC, 91.5%
coding density, 11.2% hypothetical, and a case split of 2816 / 24 / 1097 /
11 / 97. Twelve samples across four MTBC lineages, so the report has
something to compare between groups.

Deterministic: a fixed seed, so re-rendering the demo gives the same figures.
"""
import os
import random
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else "demo_inputs"
os.makedirs(OUT, exist_ok=True)
rng = random.Random(20260725)

LINEAGES = ["L1", "L2", "L3", "L4"]
SAMPLES = []
for li, lin in enumerate(LINEAGES):
    for k in range(3):
        SAMPLES.append((lin, f"G{10000 + li * 100 + k * 7}"))

COGS = list("CDEFGHIJKLMNOPQSTUV")
COG_WEIGHT = {  # roughly how MTB distributes, so the bars are not flat
    "I": 0.10, "Q": 0.09, "E": 0.08, "C": 0.08, "K": 0.07, "G": 0.06,
    "M": 0.06, "P": 0.05, "L": 0.05, "J": 0.05, "O": 0.04, "T": 0.04,
    "H": 0.04, "F": 0.03, "V": 0.03, "U": 0.02, "D": 0.02, "N": 0.02,
    "S": 0.07,
}
MYCO = [
    ("intermediary metabolism and respiration", 0.225),
    ("conserved hypotheticals", 0.245),
    ("cell wall and cell processes", 0.185),
    ("lipid metabolism", 0.065),
    ("information pathways", 0.060),
    ("regulatory proteins", 0.050),
    ("PE/PPE", 0.042),
    ("virulence, detoxification, adaptation", 0.055),
    ("insertion seqs and phages", 0.073),
]
KEGG = ["ko00010", "ko00020", "ko00190", "ko00230", "ko00240", "ko00250",
        "ko00260", "ko00330", "ko00520", "ko00550", "ko00620", "ko00640",
        "ko00860", "ko01100", "ko02010", "ko03010"]


def jitter(base, pct):
    return int(round(base * (1 + rng.uniform(-pct, pct))))


def w(name, lines):
    with open(os.path.join(OUT, name), "w") as fh:
        fh.write("\n".join(lines) + "\n")


cog_rows, myco_rows, summary_rows, noref_rows, rvmiss_rows, qc_rows = [], [], [], [], [], []
master_rows = []
per_sample_cds = {}

qc_rows.append("\t".join([
    "sample", "genome_size", "n_contigs", "n50", "gc_pct", "total_cds",
    "pseudogene_count", "hypothetical_pct", "mean_gene_length",
    "median_gene_length", "coding_density_pct", "eggnog_annotation_pct"]))

# as REPORT_INPUTS writes it: the sample each gene is from, then its master table
master_header = ["sample", "simple_id", "lineage", "start", "end", "strand", "rvs",
                 "pseudogene", "COG_category", "Description", "Preferred_name",
                 "KEGG_Pathway", "Functional_Category"]
master_rows.append("\t".join(master_header))

for lin, gid in SAMPLES:
    sample = f"{lin}_{gid}"
    total_cds = jitter(4020, 0.02)
    per_sample_cds[sample] = total_cds

    # --- case summary ---
    case2 = jitter(2816, 0.04)
    case3 = jitter(24, 0.5)
    case4 = jitter(1097, 0.08)
    rvmiss = jitter(97, 0.25)
    case1_init = jitter(141, 0.3)
    case1_removed = int(case1_init * rng.uniform(0.85, 0.96))
    case1_kept = case1_init - case1_removed
    for t, n in [("case2", case2), ("case3", case3), ("case4", case4),
                 ("case1_initial", case1_init), ("rv_missing", rvmiss),
                 ("total_cds", case2 + case3 + case4 + case1_init),
                 ("case1_removed", case1_removed), ("case1_kept", case1_kept)]:
        # as TYPE_CASES writes it: the name, the lineage and the sample ID
        summary_rows.append(f"{t}\t{n}\t{sample}\t{lin}\t{gid}")

    # --- QC ---
    gsize = jitter(4403837, 0.012)
    qc_rows.append("\t".join(map(str, [
        sample, gsize, 1, gsize, round(rng.uniform(65.4, 65.8), 2), total_cds,
        jitter(13, 0.6), round(rng.uniform(10.2, 12.4), 1), jitter(1002, 0.03),
        jitter(846, 0.03), round(rng.uniform(90.6, 92.1), 1),
        round(rng.uniform(88.5, 94.0), 1)])))

    # --- COG counts ---
    assigned = int(total_cds * rng.uniform(0.88, 0.93))
    for c in COGS:
        n = max(1, int(assigned * COG_WEIGHT[c] * rng.uniform(0.88, 1.12)))
        cog_rows.append(f"{lin}\t{sample}\t{c}\t{n}\t{n / assigned:.6f}")

    # --- Mycobrowser counts ---
    for cat, frac in MYCO:
        n = max(1, int(total_cds * frac * rng.uniform(0.92, 1.08)))
        myco_rows.append(f"{lin}\t{sample}\t{cat}\t{n}\t{n / total_cds:.6f}")

    # --- genes with no reference match ---
    for _ in range(case1_kept):
        noref_rows.append(f"{gid}\t{lin}")

    # --- missing reference genes ---
    for _ in range(rvmiss):
        rv = f"Rv{rng.randint(1, 3924):04d}{rng.choice(['', '', '', 'c', 'A'])}"
        s = rng.randint(1, 4_400_000)
        ln = rng.randrange(150, 3000, 3)
        rvmiss_rows.append("\t".join(map(str, [
            f"{rv}_{s}-{s + ln}", ln // 3, rng.choice("+-"), s, s + ln,
            int(ln * rng.uniform(0.3, 0.95)), ln, sample])))

    # --- master table, a slice per sample so PE/PPE and KEGG have something ---
    for i in range(240):
        cat = rng.choices([m[0] for m in MYCO], [m[1] for m in MYCO])[0]
        s = rng.randint(1, 4_400_000)
        master_rows.append("\t".join([
            sample, f"{sample}_{i:05d}", lin, str(s), str(s + rng.randrange(150, 3000, 3)),
            rng.choice("+-"), f"Rv{rng.randint(1, 3924):04d}",
            rng.choice(["no"] * 24 + ["yes"]), rng.choice(COGS),
            "hypothetical protein" if cat == "conserved hypotheticals" else "putative enzyme",
            "-", rng.choice(KEGG + ["-"] * 6), cat]))

w("input_cog.txt", cog_rows)
w("input_myco.txt", myco_rows)
w("input_summary.txt", summary_rows)
w("input_noref.txt", noref_rows)
w("input_rvmissing.txt", rvmiss_rows)
w("input_qc.txt", qc_rows)
w("input_mastertables.txt", master_rows)

# --- pangenome presence/absence ---
names = [f"{lin}_{gid}" for lin, gid in SAMPLES]
pan = ["Gene\t" + "\t".join(names)]
for i in range(3600):
    if i < 2900:                       # core
        row = [1] * len(names)
    elif i < 3200:                     # softcore / shell
        row = [1 if rng.random() < 0.75 else 0 for _ in names]
    else:                              # cloud
        row = [1 if rng.random() < 0.18 else 0 for _ in names]
    pan.append(f"group_{i:05d}\t" + "\t".join(map(str, row)))
w("input_pangene.txt", pan)

# --- homologous communities across genomes, as GET_ALLDUP writes them: a row
# per pair, a community of one genome and its homolog in another, each with its
# size, its proteins' COG categories, its network metrics and the proteins.
# A protein eggNOG gave no category keeps the ID FIX_CDS wrote and counts as '-'.
FAMILY_COGS = ["L", "L", "L", "S", "K", "M", "I", "Q"]  # transposases first


def community(genome, n, cog):
    """A community of n proteins of a family in category cog: its size, COG
    counts, metrics and proteins."""
    nodes = []
    for locus in rng.sample(range(1, 4100), n):
        if rng.random() < 0.8:
            nodes.append(f"{genome}_{locus:05d};{cog}")
        else:
            s = rng.randint(1, 4_400_000)
            nodes.append(f"{genome}_{locus:05d};{rng.choice('+-')};{s}-{s + rng.randrange(150, 3000, 3)}")
    cats = [p.split(";")[1] if p.count(";") == 1 else "-" for p in nodes]
    counts = {c: cats.count(c) for c in sorted(set(cats), key=lambda c: (-cats.count(c), c))}
    if n == 2 or rng.random() < 0.7:   # every protein matches every other
        metrics = [float(n - 1), 1.0, 1.0 if n > 2 else 0.0, 1.0]
    else:
        density = round(rng.uniform(0.5, 0.9), 3)
        metrics = [round(density * (n - 1), 3), density,
                   round(rng.uniform(0.4, 0.9), 3), round(rng.uniform(1.1, 1.9), 3)]
    return [n, counts] + metrics + [", ".join(sorted(nodes))]


dup = ["\t".join(f"{c}{k}" for k in "12" for c in [
    "Genome", "Community", "Size", "Elements", "AvgDegree", "Density",
    "AvgClustering", "AvgPathLength", "Nodes"])]
for i in range(48):
    g1, g2 = rng.sample(names, 2)
    cog = rng.choice(FAMILY_COGS)
    s1 = rng.choice([2, 2, 3, 3, 4, 5, 6, 8, 11, 14, 19])
    s2 = max(2, s1 + rng.randint(-1, 3))
    dup.append("\t".join(map(str, [g1, f"Community_{rng.randint(0, 40)}"] + community(g1, s1, cog)
                              + [g2, f"Community_{rng.randint(0, 40)}"] + community(g2, s2, cog))))
w("homologous_communities.tsv", dup)

print(f"  {len(SAMPLES)} samples across {len(LINEAGES)} lineages")
print(f"  {sum(per_sample_cds.values())} CDS in total")
for f in sorted(os.listdir(OUT)):
    p = os.path.join(OUT, f)
    print(f"    {f:32s} {os.path.getsize(p):>9,} bytes")
