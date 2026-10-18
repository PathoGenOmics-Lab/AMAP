#!/usr/bin/env python3
"""Check the website in docs/ against the pipeline it describes.

The site is plain HTML with no build step, so nothing renders it and nothing
fails when it drifts. This does, in seven parts:

  links       every internal link and #anchor resolves, and no page loads a
              stylesheet, script or image from another host (the site has to
              open from a clone, offline);
  docs        the documentation pages are what tools/build_docs.py makes of the
              markdown today, and the markdown has no MkDocs syntax left that
              would print as text;
  options     every --option the pages or the command builder use is a
              parameter of main.nf or nextflow.config, and every -profile is
              one nextflow.config defines;
  drawings    every window of real genes in assets/figure-data.js, put through
              bin/classify, gives the cases and the kept genes it draws;
  thresholds  the numbers the drawings and the correction algorithm page give
              for each rule are the ones in bin/;
  metro       the metro map's stations are main.nf's processes, in the order
              its workflow runs them;
  readme      README.md's links to the website land on a page and an anchor
              that exist;
  dashes      no em-dashes.

Run from the repository root (the drawings need pandas and intervaltree):

    python3 tools/check_site.py

Exits 1 and lists every problem it found.
"""

from __future__ import annotations

import contextlib
import io
import json
import os
import re
import sys
import tempfile
from html.parser import HTMLParser
from pathlib import Path

SITE = Path("docs")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_docs  # noqa: E402  (tools/build_docs.py, next to this file)

FIGURE_DATA = SITE / "assets" / "figure-data.js"
FIGURES_JS = SITE / "assets" / "figures.js"
SITE_JS = SITE / "assets" / "site.js"
HOW_JS = SITE / "assets" / "how.js"
MAIN_NF = Path("main.nf")
CONFIG = Path("nextflow.config")
ALGORITHM_MD = SITE / "how-it-works" / "correction-algorithm.md"
README = Path("README.md")
SITE_URL = "https://pathogenomics-lab.github.io/AMAP/"
# Where Pages serves the site; the 404 page, served at any depth, links from it.
SITE_PATH = "/AMAP/"
# Nextflow's own options, which take one dash; -process.NAME sets any process
# setting from the command line.
NEXTFLOW_OPTIONS = {"-profile", "-resume", "-work-dir", "-c", "-bg", "-r", "-latest", "-with-report", "-with-trace", "-with-timeline"}


class Page(HTMLParser):
    """The parts of a page the checks read."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.ids: set[str] = set()
        self.links: list[tuple[str, str]] = []  # (tag, url)
        self.remote_assets: list[str] = []
        self.code_elems: list[tuple[str, str]] = []  # (tag, text) of each <code> and <pre>
        self.text: list[str] = []
        self._code_open: list[tuple[str, list[str]]] = []
        self._skip = 0

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if "id" in a:
            self.ids.add(a["id"])
        for key in ("href", "src"):
            if key in a and a[key] is not None:
                self.links.append((tag, a[key]))
        remote = lambda url: url and re.match(r"^(https?:)?//", url)
        if tag == "link" and "stylesheet" in (a.get("rel") or "") and remote(a.get("href")):
            self.remote_assets.append(a["href"])
        if tag in ("script", "img", "iframe") and remote(a.get("src")):
            self.remote_assets.append(a["src"])
        if tag in ("code", "pre"):
            self._code_open.append((tag, []))
        if tag in ("script", "style"):
            self._skip += 1

    def handle_endtag(self, tag):
        if tag in ("code", "pre") and self._code_open and self._code_open[-1][0] == tag:
            kind, parts = self._code_open.pop()
            self.code_elems.append((kind, "".join(parts)))
        if tag in ("script", "style"):
            self._skip = max(0, self._skip - 1)

    def handle_data(self, data):
        for _, parts in self._code_open:
            parts.append(data)
        if not self._skip:
            self.text.append(data)


def parse(path: Path) -> Page:
    page = Page()
    page.feed(path.read_text(encoding="utf-8"))
    return page


def all_pages() -> dict[str, Page]:
    """Every page of the site, by its path under docs/. The example report is
    the report's own output, rendered from the Rmd, and is not checked here."""
    out = {}
    for path in sorted(SITE.rglob("*.html")):
        rel = path.relative_to(SITE).as_posix()
        if rel != "examples/demo_report.html":
            out[rel] = parse(path)
    return out


# ── links ──


def check_links(pages: dict[str, Page]) -> list[str]:
    problems = []
    for name, page in pages.items():
        for url in page.remote_assets:
            problems.append(f"{name}: loads {url} from another host; the site must work offline")
        for tag, url in page.links:
            if re.match(r"^(https?:|mailto:)", url):
                continue
            target, _, anchor = url.partition("#")
            if not target:
                if anchor and anchor not in page.ids:
                    problems.append(f"{name}: link to #{anchor}, which no element on the page has")
                continue
            if target.startswith(SITE_PATH):
                target = "/" + target[len(SITE_PATH):]
                resolved = target.lstrip("/") or "index.html"
            else:
                resolved = os.path.normpath(os.path.join(os.path.dirname(name), target)).replace(os.sep, "/")
            if resolved.startswith("..") or not (SITE / resolved).exists():
                problems.append(f"{name}: <{tag}> points at {target}, which does not exist")
                continue
            if anchor and resolved in pages and anchor not in pages[resolved].ids:
                problems.append(f"{name}: link to {target}#{anchor}, which no element there has")
    return problems


# ── docs ──

# MkDocs admonitions, tabs, attribute lists, snippets, definitions and
# captions: GitHub and the Docs tab print them as text.
MKDOCS_LEFTOVER = re.compile(r"\{ ?[.#:]|\{ ?target=|^!!! |^\?\?\? |^=== \"|--8<--|^: |^///|markdown>|:material-|:octicons-", re.M)


def check_docs_built() -> list[str]:
    stale = [
        path.as_posix()
        for path, text in build_docs.build().items()
        if not path.exists() or path.read_text(encoding="utf-8") != text
    ]
    problems = [f"{p}: does not match its markdown; run python3 tools/build_docs.py" for p in stale]
    for _section, path, _title, _line in build_docs.PAGES:
        if path.endswith(".md"):
            in_fence = False
            for number, line in enumerate((SITE / path).read_text(encoding="utf-8").splitlines(), 1):
                if line.lstrip().startswith("```"):
                    in_fence = not in_fence
                if not in_fence and MKDOCS_LEFTOVER.search(line):
                    problems.append(f"docs/{path}:{number}: MkDocs syntax that prints as text: {line.strip()[:60]}")
    return problems


# ── options ──


def pipeline_params() -> set[str]:
    """Every parameter main.nf or nextflow.config defines, and every option
    the pipeline's own --help lists."""
    main = MAIN_NF.read_text(encoding="utf-8")
    names = set(re.findall(r"params\.([a-z_]+)\s*=", main))
    help_text = main[main.index("def helpMessage()"):]
    names |= set(re.findall(r"^\s+--([a-z_]+)\s", help_text[: help_text.index('"""', help_text.index('"""') + 3)], re.M))
    config = CONFIG.read_text(encoding="utf-8")
    names |= set(re.findall(r"params\.([a-z_]+)\s*=", config))
    block = re.search(r"^params\s*\{(.*?)^\}", config, re.M | re.S)
    if block:
        names |= set(re.findall(r"^\s*([a-z_]+)\s*=", block.group(1), re.M))
    return names


def profiles() -> set[str]:
    config = CONFIG.read_text(encoding="utf-8")
    block = config[config.index("profiles {"):]
    return set(re.findall(r"^    ([a-z]+)\s*\{", block, re.M))


def check_options(pages: dict[str, Page]) -> list[str]:
    problems = []
    params = pipeline_params()
    known_profiles = profiles()
    seen: dict[str, str] = {}
    # An option named on its own may be one of a script's in bin/ (the
    # validation step's --strict); one in a `nextflow run` line may not.
    script_options = set()
    for script in Path("bin").rglob("*.py"):
        script_options |= set(re.findall(r"add_argument\([^)]*?[\"'](--[a-z][\w-]*)[\"']", script.read_text(encoding="utf-8")))
    alone: dict[str, str] = {}

    def note(command: str, where: str) -> None:
        # What sits in quotes is for another program: a scheduler's own flags
        # in -process.clusterOptions, a species name.
        command = re.sub(r"'[^']*'|\"[^\"]*\"", "", command)
        for flag in re.findall(r"(?<![\w-])--([a-z][a-z_]*)", command):
            seen.setdefault("--" + flag, where)
        for option in re.findall(r"(?<![\w-])(-[a-z][\w.-]*)", command):
            if option not in NEXTFLOW_OPTIONS and not re.fullmatch(r"-process\.\w+", option):
                problems.append(f"{where}: {option} is not a Nextflow option the site knows")
        for listed in re.findall(r"-profile\s+([\w,]+)", command):
            for name in listed.split(","):
                if name not in known_profiles:
                    problems.append(f"{where}: -profile {name}, which nextflow.config does not define")

    # Command lines that run AMAP, and option names on their own in <code>.
    for name, page in pages.items():
        for tag, text in page.code_elems:
            if tag == "pre":
                for command in text.replace("\\\n", " ").splitlines():
                    command = command.strip().removeprefix("$ ")
                    if command.startswith("nextflow run"):
                        note(command, name)
            elif re.fullmatch(r"--[a-z_]+", text.strip()):
                alone.setdefault(text.strip(), name)
    for _section, path, _title, _line in build_docs.PAGES:
        if path.endswith(".md"):
            text = (SITE / path).read_text(encoding="utf-8")
            for block in re.findall(r"```[^\n]*\n(.*?)```", text, re.S):
                for command in block.replace("\\\n", " ").splitlines():
                    command = command.strip().removeprefix("$ ")
                    if command.startswith("nextflow run"):
                        note(command, f"docs/{path}")
    # The command builder writes its lines as [["--flag", "value"], "what it does"].
    js = SITE_JS.read_text(encoding="utf-8")
    for flags in re.findall(r"out\.push\(\[\[([^\]]*)\]", js):
        for flag in re.findall(r'"(--?[a-z][\w.]*)"', flags):
            if flag.startswith("--"):
                seen.setdefault(flag, "site.js")
    for profile in re.findall(r'profile: "([\w,]+)"', js):
        for name in profile.split(","):
            if name not in known_profiles:
                problems.append(f"site.js: profile {name}, which nextflow.config does not define")
    for flag, where in sorted(seen.items()):
        if flag[2:] not in params:
            problems.append(f"{where}: {flag} is not a parameter of main.nf or nextflow.config")
    for flag, where in sorted(alone.items()):
        if flag[2:] not in params and flag not in script_options:
            problems.append(f"{where}: {flag} is neither a parameter of the pipeline nor an option of a script in bin/")
    return problems


# ── drawings ──


def figure_data() -> dict:
    text = FIGURE_DATA.read_text(encoding="utf-8")
    start = text.index("window.AMAP_FIGURE_DATA =") + len("window.AMAP_FIGURE_DATA =")
    return json.loads(text[start:].strip().rstrip(";"))


def run_window(win: dict, contig: str, base: str) -> tuple[dict[str, int], list[tuple]]:
    """The cases and the corrected genes bin/classify gives a window, run the
    way classify_cases.py runs a genome from step 2 on."""
    import pandas as pd

    sys.path.insert(0, str(ROOT / "bin"))
    from classify.cases import classify_case1_case4, classify_case2_case3, filter_case1, find_rv_missing
    from classify.overlaps import build_overlap_tree, create_overlap_groups, resolve_multiple_overlaps, resolve_single_overlaps
    from classify.reconciliation import deduplicate, recover_missing_genes, write_corrected_tsv

    bakta = pd.DataFrame(
        [{"Locus Tag": c["tag"], "Start": c["s"], "End": c["e"], "Strand": c["st"], "Length_nt": c["e"] - c["s"] + 1, "Contig": contig}
         for c in win["bakta"]],
        columns=["Locus Tag", "Start", "End", "Strand", "Length_nt", "Contig"],
    )
    mp = pd.DataFrame(
        [{"Prot Id": r["prot"], "Prot len": r["len"], "Strand": r["st"], "Start": r["s"], "End": r["e"], "Match Nt": r["match"],
          "Num nt": r["num"], "Contig": contig, "Complete": r["complete"]} for r in win["ref"]],
        columns=["Prot Id", "Prot len", "Strand", "Start", "End", "Match Nt", "Num nt", "Contig", "Complete"],
    )
    bad = pd.DataFrame(
        [{"Prot Id": f["prot"], "Prot len": 0, "Strand": f["st"], "Start": f["s"], "End": f["e"], "Match Nt": 0, "Num nt": 0, "Contig": contig}
         for f in win["bad"]],
        columns=["Prot Id", "Prot len", "Strand", "Start", "End", "Match Nt", "Num nt", "Contig"],
    )
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        # filter_case1 reads Bakta's proteins aligned onto the reference genome
        # from miniprot's PAF; only the identity of each row counts.
        paf = Path(tmp) / "ref.paf"
        paf.write_text("".join(
            f"{m['tag']}\t100\t0\t100\t+\tref\t1000000\t0\t{max(m['num'], 1)}\t{m['match']}\t{m['num']}\t0\n" for m in win["refMatch"]
        ))
        case2, case3 = classify_case2_case3(bakta, mp)
        # A whole genome always has calls in case 1 or 4; a window of one gene
        # may not, and classify_case1_case4 cannot take an empty case (#39).
        exact = {(r.Start, r.End, r.Strand) for r in mp.itertuples()}
        if all((c["s"], c["e"], c["st"]) in exact for c in win["bakta"]):
            case4 = case1 = bakta.iloc[0:0]
        else:
            case4, case1 = classify_case1_case4(bakta, mp)
        # Nor can find_rv_missing take a window where nothing aligns.
        missing = find_rv_missing(mp, bakta) if len(mp) else mp
        keep, _removed = filter_case1(case1, bad, str(paf))
        cases = {}
        for number, frame in ((1, case1), (2, case2), (3, case3), (4, case4)):
            for tag in frame["Locus Tag"]:
                cases[tag] = number

        combined: dict = {}
        for frame, col in ((case2, "Locus Tag"), (case3, "Locus Tag"), (missing, "Prot Id"), (keep, "Locus Tag")):
            for _, row in frame.iterrows():
                combined.setdefault((row["Start"], row["End"], row.get("Strand", "+"), row.get("Contig", "")), row[col])
        if len(case4):
            mp_case4 = build_overlap_tree(mp, case4, "Prot Id", "Locus Tag")
            case4_mp = build_overlap_tree(case4, mp, "Locus Tag", "Prot Id")
            single = create_overlap_groups({k: v for k, v in mp_case4.items() if v["count"] == 1}, case4, mp)
            multiple = create_overlap_groups({k: v for k, v in mp_case4.items() if v["count"] > 1}, case4, mp)
            groups = {v["Prot Id"]: v for v in create_overlap_groups(case4_mp, case4, mp).values()}
            resolve_single_overlaps(single, mp, combined, groups)
            resolve_multiple_overlaps(multiple, mp, combined, groups)
        recover_missing_genes(combined, mp)
        deduplicate(combined)
        os.chdir(tmp)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                write_corrected_tsv(combined, mp, bakta, base)
        finally:
            os.chdir(cwd)
        kept = []
        for line in (Path(tmp) / f"{base}_corrected.tsv").read_text().splitlines():
            gene_id, strand, start, end, rvs, _contig = line.split("\t")
            kept.append((gene_id, strand, int(start), int(end), rvs))
    return cases, kept


def check_drawings() -> list[str]:
    problems = []
    windows = figure_data()["windows"]
    for name, win in windows.items():
        kind = win.get("kind", "window")
        try:
            if kind == "synteny":
                problems += check_synteny(name, win)
                continue
            if kind == "copies":
                problems += check_copies(name, win)
                continue
            cases, kept = run_window(win, win["contig"], win["base"])
        except ImportError as error:
            return [f"assets/figure-data.js: cannot run bin/ ({error}); pip install pandas intervaltree biopython"]
        for call in win["bakta"]:
            if cases.get(call["tag"]) != call["case"]:
                problems.append(
                    f"assets/figure-data.js: {name}: {call['tag']} is drawn as case {call['case']}, bin/classify says {cases.get(call['tag'])}"
                )
        corrected = sorted((c["id"], c["st"], c["s"], c["e"], c["rvs"]) for c in win["corrected"])
        if sorted(kept) != corrected:
            problems.append(f"assets/figure-data.js: {name}: has the classification keeping {corrected}, bin/classify keeps {sorted(kept)}")
        # Then check_cds.py translates every gene. Where it changed none, what
        # AMAP keeps is what the classification wrote.
        if "translate" in win:
            problems += check_translate(name, win)
        else:
            drawn = sorted((a["id"], a["st"], a["s"], a["e"], a["rvs"]) for a in win["amap"])
            if drawn != corrected or any(a["pseudo"] or a["partial"] or a["note"] for a in win["amap"]):
                problems.append(f"assets/figure-data.js: {name}: AMAP keeps something check_cds.py would have to change; add the genome it reads")
        if "exact" in win:
            problems += check_exact(name, win)
    return problems


def check_synteny(name: str, win: dict) -> list[str]:
    """The placements drawn as dropped are the ones drop_paralog_placements
    drops, given the proteins placed once around each locus."""
    import pandas as pd

    sys.path.insert(0, str(ROOT / "bin"))
    from classify.synteny import drop_paralog_placements

    df = pd.DataFrame([
        {"Prot Id": p["prot"], "Prot len": p["len"], "Strand": p["st"], "Start": p["s"], "End": p["e"], "Match Nt": p["match"],
         "Num nt": p["num"], "Contig": win["contig"], "Complete": p["complete"], "Ref Start": float(p["refStart"]),
         "Ref End": float(p["refEnd"])}
        for p in win["placements"]
    ])
    out = drop_paralog_placements(df)
    kept = set(zip(out["Prot Id"], out["Start"], out["End"]))
    problems = []
    for p in win["placements"]:
        if ((p["prot"], p["s"], p["e"]) in kept) != p["kept"]:
            problems.append(f"assets/figure-data.js: {name}: draws {p['prot']} at {p['s']} as {'kept' if p['kept'] else 'dropped'}, bin/classify says otherwise")
    return problems


def check_copies(name: str, win: dict) -> list[str]:
    """deduplicate_rvs.py keeps every copy drawn, and names the later ones _2."""
    sys.path.insert(0, str(ROOT / "bin"))
    import deduplicate_rvs

    # What check_cds.py hands it ("before"), where every copy only miniprot
    # placed carries the protein's name; what it keeps is drawn ("genes").
    genes = [
        {"id": g["id"], "strand": g["st"], "start": str(g["s"]), "end": str(g["e"]), "rvs": g["rvs"],
         "pseudogene": "yes" if g["pseudo"] else "no", "contig": win["contig"]}
        for g in win["before"]
    ]
    bakta_cds = {(c["s"], c["e"], c["st"], win["contig"]) for c in win["bakta"]}
    kept, _removed = deduplicate_rvs.deduplicate(genes, bakta_cds, {}, {win["prot"]: win["refResidues"]})
    for i, new_id in deduplicate_rvs.unique_ids(kept).items():
        kept[i]["id"] = new_id
    got = sorted((g["id"], int(g["start"]), int(g["end"])) for g in kept)
    want = sorted((g["id"], g["s"], g["e"]) for g in win["genes"])
    if got != want:
        return [f"assets/figure-data.js: {name}: draws {want}, deduplicate_rvs.py keeps {got}"]
    return []


def check_exact(name: str, win: dict) -> list[str]:
    """exact.py finds the protein where the drawing puts it."""
    sys.path.insert(0, str(ROOT / "bin"))
    from classify.exact import _six_frames, find_exact

    ex = win["exact"]
    with tempfile.TemporaryDirectory() as tmp:
        fasta = Path(tmp) / "piece.fasta"
        fasta.write_text(f">{win['contig']}\n{ex['dna']}\n")
        found = [(s + ex["offset"], e + ex["offset"], st) for _contig, s, e, st in find_exact(ex["protein"], _six_frames(fasta))]
    gene = next((a for a in win["amap"] if a["rvs"] == ex["prot"]), None)
    if gene is None or (gene["s"], gene["e"], gene["st"]) not in found:
        return [f"assets/figure-data.js: {name}: exact.py finds {ex['prot']} at {found}, the drawing keeps {gene}"]
    return []


def check_translate(name: str, win: dict) -> list[str]:
    """check_cds.py, run on the classification's genes over the stretch of
    genome the window carries, gives the coordinates, the pseudogenes, the
    partial ends and the notes drawn."""
    import pandas as pd

    sys.path.insert(0, str(ROOT / "bin"))
    import logging

    logging.disable(logging.CRITICAL)
    import check_cds

    tr = win["translate"]
    off = tr["offset"]
    rows = [{"id": c["id"], "strand": c["st"], "start": c["s"] - off, "end": c["e"] - off, "rvs": c["rvs"], "contig": win["contig"]}
            for c in win["corrected"]]
    refs = {p: (seq.upper(), span) for p, (seq, span) in tr["proteins"].items()}
    with tempfile.TemporaryDirectory() as tmp:
        fasta = Path(tmp) / "piece.fasta"
        fasta.write_text(f">{win['contig']}\n{tr['dna']}\n")
        _aa, _nt, out = check_cds.extract_and_translate_cds(str(fasta), pd.DataFrame(rows), refs)
    got = sorted((r["id"], int(r["start"]) + off, int(r["end"]) + off, r["pseudogene"] == "yes", r["partial"], r["note"])
                 for r in out.to_dict("records"))
    # A further copy's _2 is deduplicate_rvs.py's doing, checked with the copies.
    ids = {g[0] for g in got}
    plain = lambda gene_id: re.sub(r"_\d$", "", gene_id) if gene_id not in ids and re.sub(r"_\d$", "", gene_id) in ids else gene_id
    drawn = sorted((plain(a["id"]), a["s"], a["e"], a["pseudo"], a["partial"], a["note"]) for a in win["amap"])
    if got != drawn:
        return [f"assets/figure-data.js: {name}: draws {drawn}, check_cds.py gives {got}"]
    return []


# ── thresholds ──


def code_value(path: str, pattern: str) -> float | None:
    """The number a pattern captures in a file of bin/, or None."""
    m = re.search(pattern, (ROOT / path).read_text(encoding="utf-8"), re.M)
    return float(m.group(1).replace("_", "")) if m else None


# What each rule's number is in the code. A pattern that stops matching is a
# problem in its own right: the code moved and this table has to follow it.
RULES = {
    "single_bp": ("bin/classify/overlaps.py", r"abs\(bakta_len - rv_len\) > (\d+)"),
    "multiple_bp": ("bin/classify/overlaps.py", r"\(_range_length\(best\[0\]\) - _range_length\(miniprot_range\)\) < (\d+)"),
    "insert_overlap_pct": ("bin/classify/overlaps.py", r"if pct >= (\d+):"),
    "insert_margin_bp": ("bin/classify/overlaps.py", r"def _insert_with_dedup\(combined_dict, new_coord, new_id, threshold=(\d+)"),
    "false_positive_pct": ("bin/classify/cases.py", r"case1\['Overlap_Percentage'\] >= (\d+)"),
    "identity_pct": ("bin/classify/cases.py", r"matches\['match_per'\] >= (\d+)"),
    "account_cover": ("bin/classify/reconciliation.py", r"def recover_missing_genes\(combined_dict, miniprot_df, min_cover=([\d.]+)\)"),
    "synteny_bp": ("bin/classify/synteny.py", r"^SYNTENY_TOLERANCE = ([\d_]+)"),
    "anchors": ("bin/classify/synteny.py", r"^N_ANCHORS = (\d+)"),
    "claim_cover": ("bin/classify/synteny.py", r"^MIN_COVER = ([\d.]+)"),
    "exact_residues": ("bin/classify/exact.py", r"^MIN_RESIDUES = (\d+)"),
    "copy_length": ("bin/deduplicate_rvs.py", r"^MIN_COPY_LENGTH = ([\d.]+)"),
    "stop_search_bp": ("bin/check_cds.py", r"reach = (\d+)"),
    "past_reference_residues": ("bin/check_cds.py", r"^COMPLETION_SLACK = (\d+)"),
}


def check_thresholds() -> list[str]:
    problems = []
    values = {}
    for rule, (path, pattern) in RULES.items():
        value = code_value(path, pattern)
        if value is None:
            problems.append(f"tools/check_site.py: the rule {rule} is no longer where it was in {path}; update RULES")
        values[rule] = value
    # What the pages say each number is.
    said = {
        "single_bp": [("index.html", r"longer by more than ([\d,]+) bp")],
        "false_positive_pct": [("index.html", r"genome at (\d+)% identity")],
    }
    for rule, places in said.items():
        for page, pattern in places:
            for m in re.finditer(pattern, (SITE / page).read_text(encoding="utf-8")):
                if values.get(rule) is not None and float(m.group(1).replace(",", "")) != values[rule]:
                    problems.append(f"{page}: says {m.group(0)!r}, the code says {values[rule]:g}")
    if HOW_JS.exists():
        block = re.search(r"var RULES = (\{.*?\});", HOW_JS.read_text(encoding="utf-8"), re.S)
        if not block:
            problems.append("assets/how.js: no RULES block to check")
        else:
            for rule, number in json.loads(block.group(1)).items():
                if rule not in values:
                    problems.append(f"assets/how.js: RULES has {rule}, which tools/check_site.py does not know")
                elif values[rule] is not None and float(number) != values[rule]:
                    problems.append(f"assets/how.js: {rule} is {number}, the code says {values[rule]:g}")
    # The threshold table of the correction algorithm page.
    table = {
        ("Bakta length advantage", "single overlaps"): ("single_bp", "{:g} bp"),
        ("Bakta length advantage", "multiple overlaps"): ("multiple_bp", "{:g} bp"),
        ("False-positive overlap", "case 1 filtering"): ("false_positive_pct", "{:g}%"),
        ("Reference match identity", "case 1 filtering"): ("identity_pct", "{:g}%"),
        ("Overlap for deduplication", "on insert"): ("insert_overlap_pct", "{:g}%"),
        ("Length margin for displacement", "on insert"): ("insert_margin_bp", "{:g} bp"),
        ("Distance from the predicted position", "paralog placements"): ("synteny_bp", "{:g} kb", 1000),
        ("Proteins placed once, each side, for the local offset", "paralog placements"): ("anchors", "{:g}"),
        ("Cover to claim another protein's placement", "paralog placements"): ("claim_cover", "{:g}%", 0.01),
        ("Shortest protein looked for exactly", "small proteins"): ("exact_residues", "{:g} residues"),
        ("Cover for a gene to account for a placement", "rescuing missing genes"): ("account_cover", "{:g}%", 0.01),
        ("Length of a further copy, next to the first", "final deduplication"): ("copy_length", "{:g}%", 0.01),
        ("Residues an extension may run past the end of the reference protein", "quality control"): ("past_reference_residues", "{:g}"),
    }
    for line in ALGORITHM_MD.read_text(encoding="utf-8").splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) != 3 or (cells[0], cells[1]) not in table:
            continue
        spec = table[(cells[0], cells[1])]
        rule, form = spec[0], spec[1]
        scale = spec[2] if len(spec) > 2 else 1
        if values.get(rule) is None:
            continue
        expected = form.format(values[rule] / scale)
        if cells[2] != expected:
            problems.append(f"{ALGORITHM_MD}: {cells[0]} ({cells[1]}) reads {cells[2]!r}, the code says {expected!r}")
    return problems


# ── metro ──


def js_block(source: str, name: str) -> str:
    m = re.search(rf"var {name} = (\[.*?\]);", source, re.S)
    return m.group(1) if m else ""


def check_metro() -> list[str]:
    problems = []
    js = FIGURES_JS.read_text(encoding="utf-8")
    order = re.findall(r'"([A-Z_]+)"', js_block(js, "METRO_ORDER"))
    main = MAIN_NF.read_text(encoding="utf-8")
    workflow = main[main.index("\nworkflow {"):]
    workflow = workflow[: workflow.index("\n}\n")]
    defined = set(re.findall(r"include \{\s*([A-Z_; ]+?)\s*\}", main))
    names = set()
    for group in defined:
        names |= {n.strip() for n in group.split(";")}
    names |= set(re.findall(r"^process ([A-Z_]+)", main, re.M))
    helpers = {"SOFTWARE_VERSIONS", "DUMP_PARAMS", "TEST_GENOME", "CREATE_EMPTY_FILE", "REPORT_INPUTS"}
    called = []
    for name in re.findall(r"\b([A-Z][A-Z_]+)\(", workflow):
        if name not in names or name in helpers:
            continue
        name = "REFERENCE" if name in ("PREPARE_REFERENCE", "BUNDLED_REFERENCE") else name
        if name not in called:
            called.append(name)
    if order != called:
        problems.append(f"assets/figures.js: METRO_ORDER is {order}, main.nf's workflow runs {called}")
    stations = re.findall(r'^    ([A-Z_]+): \{\n\s+name: "[^"]*", proc: "([^"]*)"', js, re.M)
    for key, proc in stations:
        if key in ("SHEET",):
            continue
        if key == "REFERENCE":
            if proc != "BUNDLED_REFERENCE":
                problems.append(f"assets/figures.js: the Reference station shows {proc}")
            continue
        if key not in names:
            problems.append(f"assets/figures.js: station {key} is not a process main.nf includes")
        if proc != key:
            problems.append(f"assets/figures.js: station {key} is labelled {proc}")
    return problems


# ── readme ──


def check_readme(pages: dict[str, Page]) -> list[str]:
    problems = []
    text = README.read_text(encoding="utf-8")
    for url in re.findall(r"\((" + re.escape(SITE_URL) + r"[^)\s]*)\)", text):
        rest = url[len(SITE_URL):]
        target, _, anchor = rest.partition("#")
        target = target or "index.html"
        if target.endswith("/"):
            target += "index.html"
        if target not in pages:
            problems.append(f"README.md: links to {url}, and docs/{target} does not exist")
        elif anchor and anchor not in pages[target].ids:
            problems.append(f"README.md: links to {url}, and docs/{target} has no #{anchor}")
    if "mkdocs" in text.lower():
        problems.append("README.md: still mentions MkDocs")
    return problems


# ── dashes ──


EM_DASH = "\u2014"


def check_dashes() -> list[str]:
    """No em-dashes in what the site says. What a program prints is quoted as
    it prints it, so code blocks are left out (validate_inputs.py writes one),
    and so is the search index, which is the pages' own text again."""
    problems = []
    files = [
        p for p in SITE.rglob("*")
        if p.suffix in (".html", ".md", ".js", ".css") and p.name not in ("demo_report.html", "search-index.js")
    ]
    files += [README, Path("CHANGELOG.md"), Path("tools/build_docs.py"), Path("tools/check_site.py"), Path("tools/figure_data.py")]
    for path in files:
        text = path.read_text(encoding="utf-8")
        if path.suffix == ".html":
            text = re.sub(r"<pre\b.*?</pre>", lambda m: "\n" * m.group(0).count("\n"), text, flags=re.S)
        fenced = False
        for number, line in enumerate(text.splitlines(), 1):
            if path.suffix == ".md" and line.lstrip().startswith("```"):
                fenced = not fenced
            if EM_DASH in line and not fenced:
                problems.append(f"{path}:{number}: an em-dash")
    return problems


def main() -> int:
    pages = all_pages()
    checks = [
        ("links", lambda: check_links(pages)),
        ("docs", check_docs_built),
        ("options", lambda: check_options(pages)),
        ("drawings", check_drawings),
        ("thresholds", check_thresholds),
        ("metro", check_metro),
        ("readme", lambda: check_readme(pages)),
        ("dashes", check_dashes),
    ]
    failed = 0
    for name, check in checks:
        problems = check()
        status = "ok" if not problems else f"{len(problems)} problem(s)"
        print(f"{name:11s} {status}")
        for problem in problems:
            print(f"  {problem}")
        failed += len(problems)
    if failed:
        print(f"\n{failed} problem(s) in the website", file=sys.stderr)
        return 1
    print(f"\nThe website ({len(pages)} pages) matches the pipeline.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
