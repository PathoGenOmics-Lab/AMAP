#!/usr/bin/env python3
"""Write docs/assets/dup-data.js, the duplication family and the homologous
communities the duplications page draws.

They are real: the corrected proteomes of two benchmark genomes (CDC1551 and
M. bovis AF2122/97) put through GET_DUPLICATIONS and then GET_ALLDUP, Diamond
included, by tools/dup_run.sh in the AMAP container. That script runs the
commands of modules/local/get.nf but for rename_fastas.py, which only appends
each protein's eggNOG COG letter to its name and needs eggNOG's output.

    docker run --rm --platform linux/amd64 -v "$PWD:/amap" -w /amap \\
        paururo/amap:latest bash tools/dup_run.sh dup_work \\
        benchmark_work/runs/final/CDC1551/L4_CDC1551.AMAP.faa \\
        benchmark_work/runs/final/M.bovis/La1_Mbovis.AMAP.faa
    python3 tools/dup_data.py dup_work --protein Rv2885c \\
        --genbank /reference/sequence.gb \\
        --bakta benchmark_work/inputs/mtb_cdc1551/L4_CDC1551/L4_CDC1551.tsv \\
                benchmark_work/inputs/mtb_mbovis/La1_Mbovis/La1_Mbovis.tsv

The family is the one --protein belongs to in the first genome. With it go
every hit of Diamond's self search that touches a member, which is what
get_duplications.py needs to find the family again. The communities are the
pair of homologous communities that holds it, one in each genome, with every
hit within each genome, which is what network_genomes.py needs to build the
networks and split them again. tools/check_site.py reruns both scripts on
these hits and fails when they no longer give what is drawn. Needs pandas and
networkx.
"""

from __future__ import annotations

import argparse
import ast
import json
import sys
from pathlib import Path

import networkx as nx
import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "bin"))
import network_genomes  # noqa: E402

OUT = ROOT / "docs" / "assets" / "dup-data.js"
COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gaps", "qstart", "qend", "sstart", "send",
           "evalue", "bitscore", "qlen", "slen", "nident"]
LABELS = {"L4_CDC1551": "M. tuberculosis CDC1551", "La1_Mbovis": "M. bovis AF2122/97"}
HEADER = """/* The duplication family and the homologous communities the duplications
   page draws. They are real: the corrected proteomes of CDC1551 and M. bovis
   AF2122/97 (tests/benchmark) put through GET_DUPLICATIONS and GET_ALLDUP,
   Diamond included. "family" is one family of CDC1551 with every hit of the
   self search that touches a member ([query, subject, identical residues,
   query length, subject length], indices into "names"); "communities" is the
   pair of homologous communities that holds it, with every hit within each
   genome, the hits between the two and the genome of each sequence.
   tools/check_site.py runs them through bin/get_duplications.py and
   bin/network_genomes.py and fails if anything comes out different. Written
   by tools/dup_data.py.

   Everything after the "=" is JSON. */
window.AMAP_DUP_DATA = """


def read_hits(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#", header=None, names=COLUMNS)


def short(node: str, genome: str) -> str:
    return node.split(";")[0].removeprefix(genome + "_")


def products(paths: list[Path]) -> dict[str, str]:
    """Bakta's product for each of its locus tags."""
    out = {}
    for path in paths:
        table = pd.read_csv(path, sep="\t", comment="#", header=None)
        for tag, product in zip(table[5], table[7]):
            out[str(tag)] = product if isinstance(product, str) else ""
    return out


def reference_products(genbank: Path) -> dict[str, tuple[str, str]]:
    out = {}
    for record in SeqIO.parse(genbank, "genbank"):
        for f in record.features:
            if f.type == "CDS" and "locus_tag" in f.qualifiers:
                tag = f.qualifiers["locus_tag"][0]
                out[tag] = (f.qualifiers.get("gene", [""])[0], f.qualifiers.get("product", [""])[0])
    return out


def describe(node: str, genome: str, ref: dict, bakta: dict, length: int) -> dict:
    name = short(node, genome)
    gene, product = ref.get(name, ("", ""))
    if not product:
        product = bakta.get(f"{genome}_{name}", "")
    return {"id": name, "gene": gene if gene != name else "", "product": product, "aa": length}


def lengths(hits: pd.DataFrame) -> dict[str, int]:
    """Each protein's length, as Diamond reports it."""
    out = dict(zip(hits["qseqid"], hits["qlen"]))
    out.update(zip(hits["sseqid"], hits["slen"]))
    return {k: int(v) for k, v in out.items()}


def crossings(pos: dict, edges: list[tuple]) -> int:
    """How many pairs of lines cross, for picking a layout."""
    def ccw(a, b, c):
        return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])
    count = 0
    for i, (a, b) in enumerate(edges):
        for c, d in edges[i + 1:]:
            if len({a, b, c, d}) == 4:
                p, q, r, s = pos[a], pos[b], pos[c], pos[d]
                count += ccw(p, r, s) != ccw(q, r, s) and ccw(p, q, r) != ccw(p, q, s)
    return count


def grazes(pos: dict, edges: list[tuple]) -> int:
    """How many nodes sit on, or nearly on, a line between two others."""
    count = 0
    for a, b in edges:
        (ax, ay), (bx, by) = pos[a], pos[b]
        length = (bx - ax) ** 2 + (by - ay) ** 2 or 1e-9
        for n, (x, y) in pos.items():
            if n in (a, b):
                continue
            f = max(0.0, min(1.0, ((x - ax) * (bx - ax) + (y - ay) * (by - ay)) / length))
            count += ((x - ax - f * (bx - ax)) ** 2 + (y - ay - f * (by - ay)) ** 2) ** 0.5 < 0.09
    return count


def layout(graph: nx.Graph, lines: list[tuple]) -> dict:
    """Where each node sits, between 0 and 1: of 400 spring layouts with fixed
    seeds, the one where the lines drawn cross least and pass by the fewest
    nodes, and the nodes keep apart."""
    best = None
    for seed in range(400):
        pos = nx.spring_layout(graph, seed=seed, weight="weight")
        xs = [p[0] for p in pos.values()]
        ys = [p[1] for p in pos.values()]
        sx = (max(xs) - min(xs)) or 1
        sy = (max(ys) - min(ys)) or 1
        pos = {n: ((p[0] - min(xs)) / sx, (p[1] - min(ys)) / sy) for n, p in pos.items()}
        apart = min(((pos[a][0] - pos[b][0]) ** 2 + (pos[a][1] - pos[b][1]) ** 2) ** 0.5
                    for a in pos for b in pos if a < b)
        score = (crossings(pos, lines) + 2 * grazes(pos, lines), -round(apart, 2))
        if best is None or score < best[0]:
            best = (score, pos)
    return {n: (round(x, 3), round(y, 3)) for n, (x, y) in best[1].items()}


def family(work: Path, genome: str, protein: str, ref: dict, bakta: dict) -> dict:
    hits = read_hits(work / f"duplications_{genome}")
    groups = pd.read_csv(work / f"duplications_{genome}_groups.tsv", sep="\t")
    members = None
    for row in groups["Members"]:
        names = ast.literal_eval(row)  # a Python list, as get_duplications.py writes it
        if any(short(n, genome) == protein for n in names):
            members = set(names)
    if members is None:
        raise SystemExit(f"{protein} is in no duplication family of {genome}")
    touching = hits[hits["qseqid"].isin(members) | hits["sseqid"].isin(members)]
    nodes = sorted(set(touching["qseqid"]) | set(touching["sseqid"]), key=lambda n: (n not in members, short(n, genome)))
    index = {n: i for i, n in enumerate(nodes)}
    rows = [[index[q], index[s], int(n), int(ql), int(sl)]
            for q, s, n, ql, sl in zip(touching["qseqid"], touching["sseqid"], touching["nident"], touching["qlen"], touching["slen"])]
    # One line per pair of proteins: the best cover over both directions.
    pairs: dict[tuple[int, int], dict] = {}
    for q, s, n, ql, sl in rows:
        key = (min(q, s), max(q, s))
        cover = max(n / ql, n / sl) * 100
        best = pairs.get(key)
        if best is None or cover > best["cover"]:
            pairs[key] = {"a": key[0], "b": key[1], "cover": round(cover, 1), "nident": n,
                          "len": [ql, sl] if q == key[0] else [sl, ql]}
    kept_rows = touching[(touching["nident"] / touching["qlen"] >= 0.8) | (touching["nident"] / touching["slen"] >= 0.8)]
    kept = {(min(index[q], index[s]), max(index[q], index[s])) for q, s in zip(kept_rows["qseqid"], kept_rows["sseqid"])}
    # The family is laid out on its own; each protein outside it goes to the
    # nearer edge of the frame, level with the member it matches best.
    inside = [key for key in pairs if nodes[key[0]] in members and nodes[key[1]] in members]
    graph = nx.Graph()
    graph.add_nodes_from(i for i, n in enumerate(nodes) if n in members)
    for key in inside:
        graph.add_edge(*key, weight=1.0 if key in kept else 0.3)
    pos = {i: (round(0.14 + 0.72 * x, 3), y) for i, (x, y) in layout(graph, inside).items()}
    taken: list[tuple[float, float]] = []
    for i, n in enumerate(nodes):
        if n in members:
            continue
        partner = max((p for key, p in pairs.items() if i in key), key=lambda p: p["cover"])
        other = partner["a"] if partner["b"] == i else partner["b"]
        x = 0.0 if pos[other][0] < 0.5 else 1.0
        y = min(0.92, max(0.08, pos[other][1]))
        while any(abs(x - tx) < 0.01 and abs(y - ty) < 0.2 for tx, ty in taken):
            y = y + 0.22 if y < 0.6 else y - 0.22
        taken.append((x, y))
        pos[i] = (x, round(y, 3))
    size = lengths(touching)
    out_nodes = []
    for n in nodes:
        entry = describe(n, genome, ref, bakta, size[n])
        entry.update({"member": n in members, "x": pos[index[n]][0], "y": pos[index[n]][1]})
        out_nodes.append(entry)
    return {
        "genome": genome, "label": LABELS.get(genome, genome), "names": nodes, "nodes": out_nodes, "rows": rows,
        "pairs": [dict(p, kept=key in kept) for key, p in sorted(pairs.items())],
    }


def communities(work: Path, genome: str, protein: str, ref: dict, bakta: dict, placed: dict) -> dict:
    path = work / "all" / "all_vs_all.diamond"
    raw = read_hits(path)
    report = pd.read_csv(work / "all" / "homologous_communities.tsv", sep="\t")
    pick = None
    for _, row in report.iterrows():
        nodes1 = [n.strip() for n in row["Nodes1"].split(", ")]
        if row["Genome1"] == genome and any(short(n, genome) == protein for n in nodes1):
            pick = row
    if pick is None:
        raise SystemExit(f"no homologous pair of communities holds {protein} of {genome}")
    sides = []
    for k in ("1", "2"):
        g = pick[f"Genome{k}"]
        sides.append({
            "genome": g, "label": LABELS.get(g, g), "community": pick[f"Community{k}"],
            "members": [n.strip() for n in pick[f"Nodes{k}"].split(", ")],
            "metrics": {key: float(pick[f"{col}{k}"]) for key, col in (
                ("size", "Size"), ("avg_degree", "AvgDegree"), ("density", "Density"), ("avg_clustering", "AvgClustering"),
                ("avg_path_length", "AvgPathLength"))},
        })
    genomes = {s["genome"] for s in sides}
    # Which genome each sequence comes from, as GET_ALLDUP writes it from the
    # file each sample's sequences came in.
    genome_of = network_genomes.read_genomes(work / "all" / "genomes.tsv")
    gid = lambda n: genome_of[n]
    # Every hit within a genome, in the file's order (the networks are built
    # in that order, and Louvain's result depends on it), and the hits between
    # the two communities drawn.
    drawn = set(sides[0]["members"]) | set(sides[1]["members"])
    keep = [(gid(q) == gid(s) and gid(q) in genomes) or (q in drawn and s in drawn) for q, s in zip(raw["qseqid"], raw["sseqid"])]
    hits = raw[keep]
    names = list(dict.fromkeys(list(hits["qseqid"]) + list(hits["sseqid"])))
    index = {n: i for i, n in enumerate(names)}
    rows = [[index[q], index[s], int(n), int(ql), int(sl)]
            for q, s, n, ql, sl in zip(hits["qseqid"], hits["sseqid"], hits["nident"], hits["qlen"], hits["slen"])]
    # What the drawing needs: each community's links, and where its nodes sit:
    # where the family drawing put them, spread over the whole frame again.
    filtered = network_genomes.load_and_filter_blast_results(path)
    first = sides[0]
    at = {n: placed[short(n, genome)] for n in first["members"]}
    xs = [p[0] for p in at.values()]
    ys = [p[1] for p in at.values()]
    sx = (max(xs) - min(xs)) or 1
    sy = (max(ys) - min(ys)) or 1
    pos = {n: (round((p[0] - min(xs)) / sx, 3), round((p[1] - min(ys)) / sy, 3)) for n, p in at.items()}
    size = lengths(raw)
    # The second genome's proteins sit where their closest match of the first does.
    best: dict[str, tuple[int, str]] = {}
    for q, s, n in zip(raw["qseqid"], raw["sseqid"], raw["nident"]):
        if q in sides[1]["members"] and s in first["members"] and n > best.get(q, (-1, ""))[0]:
            best[q] = (n, s)
    cross = 0
    for q, s in zip(filtered["qseqid"], filtered["sseqid"]):
        if (q in sides[0]["members"] and s in sides[1]["members"]) or (q in sides[1]["members"] and s in sides[0]["members"]):
            cross += 1
    for side in sides:
        edges = set()
        for q, s in zip(filtered["qseqid"], filtered["sseqid"]):
            if q in side["members"] and s in side["members"]:
                edges.add(tuple(sorted((side["members"].index(q), side["members"].index(s)))))
        side["edges"] = sorted(edges)
        side["nodes"] = []
        for n in side["members"]:
            entry = describe(n, side["genome"], ref, bakta, size[n])
            anchor = n if side is first else best[n][1]
            entry["x"], entry["y"] = pos[anchor]
            entry["twin"] = short(anchor, first["genome"])
            side["nodes"].append(entry)
    return {"sides": sides, "cross": cross, "names": names, "genomes": [genome_of[n] for n in names], "rows": rows}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("work", type=Path, help="what tools/dup_run.sh wrote")
    parser.add_argument("--genome", default="L4_CDC1551", help="the genome whose family is drawn")
    parser.add_argument("--protein", required=True, help="a protein of that family, as AMAP names it (Rv2885c, 02906)")
    parser.add_argument("--genbank", type=Path, required=True, help="the reference GenBank, for gene names and products")
    parser.add_argument("--bakta", type=Path, nargs="+", required=True, help="each genome's Bakta table, for the products of its own calls")
    args = parser.parse_args()
    ref = reference_products(args.genbank)
    bakta = products(args.bakta)
    fam = family(args.work, args.genome, args.protein, ref, bakta)
    placed = {n["id"]: (n["x"], n["y"]) for n in fam["nodes"]}
    data = {"family": fam, "communities": communities(args.work, args.genome, args.protein, ref, bakta, placed)}
    OUT.write_text(HEADER + json.dumps(data, separators=(",", ":")) + ";\n", encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}: {OUT.stat().st_size:,} bytes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
