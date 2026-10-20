#!/usr/bin/env python3
"""Check the website in docs/ against the pipeline it describes.

The site is plain HTML with no build step, so nothing renders it and nothing
fails when it drifts. This does, in nine parts:

  links       every internal link and #anchor resolves (the stories how.js
              opens and the sections of the example report included), no page
              loads a stylesheet, script or image from another host (the site
              has to open from a clone, offline), and every page carries the
              card a shared link shows;
  docs        the documentation pages are what tools/build_docs.py makes of the
              markdown and the code today, the markdown has no MkDocs syntax
              left that would print as text, and docs/.nojekyll is there to
              keep GitHub Pages from running Jekyll over the markdown;
  options     every --option the pages or the command builder use is a
              parameter of main.nf or nextflow.config, and every -profile is
              one nextflow.config defines;
  drawings    every window of real genes in assets/figure-data.js, put through
              bin/classify, gives the cases and the kept genes it draws, and the
              duplication family and communities of assets/dup-data.js are the
              ones bin/get_duplications.py and bin/network_genomes.py find;
  thresholds  the numbers the pages give for each rule, each validation limit,
              each step's resources, the caps on a task and each tool's
              version are the ones in bin/, nextflow.config and the Dockerfile;
  benchmark   the front page's scores are tests/benchmark's expected results;
  metro       the metro map's stations are main.nf's processes, in the order
              its workflow runs them;
  readme      README.md's links land on a page and an anchor that exist, its
              pictures are in the repository, and its scores are the
              benchmark's;
  dashes      no em-dashes.

Run from the repository root, with what requirements-test.txt installs (the
drawings run bin/):

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
DUP_DATA = SITE / "assets" / "dup-data.js"
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
        self.meta: dict[str, str] = {}  # og: properties, what a shared link shows
        self._code_open: list[tuple[str, list[str]]] = []
        self._skip = 0

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if "id" in a:
            self.ids.add(a["id"])
        if tag == "meta" and a.get("property"):
            self.meta[a["property"]] = a.get("content") or ""
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


def story_keys() -> set[str]:
    """The stories of how-it-works.html, which how.js opens from the #key in its address."""
    return set(re.findall(r'^\s+key: "([\w-]+)", group:', HOW_JS.read_text(encoding="utf-8"), re.M)) if HOW_JS.exists() else set()


def check_links(pages: dict[str, Page]) -> list[str]:
    problems = []
    stories = story_keys()
    report_ids = None
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
            if anchor and resolved == "how-it-works.html" and anchor in stories:
                continue
            if anchor and resolved == "examples/demo_report.html":
                # The example report is not checked as a page, but the
                # sections the site opens it at have to be there.
                if report_ids is None:
                    report_ids = parse(SITE / resolved).ids
                if anchor not in report_ids:
                    problems.append(f"{name}: link to {target}#{anchor}, which the example report has no section for")
                continue
            if anchor and resolved in pages and anchor not in pages[resolved].ids:
                problems.append(f"{name}: link to {target}#{anchor}, which no element there has")
        # What a link to the page shows where it is shared: the card has to be
        # there, and the address has to be the page's own.
        for prop in ("og:title", "og:description", "og:image", "og:url"):
            if not page.meta.get(prop):
                problems.append(f"{name}: no {prop}, so a shared link to it shows no card")
        image = page.meta.get("og:image", "")
        if image and not (image.startswith(SITE_URL) and (SITE / image[len(SITE_URL):]).exists()):
            problems.append(f"{name}: og:image {image} is not a file of the site")
        url = page.meta.get("og:url")
        if url and name != "404.html" and url != SITE_URL + ("" if name == "index.html" else name):
            problems.append(f"{name}: og:url is {url}, the page is at {SITE_URL}{name}")
    return problems


# ── docs ──

# MkDocs admonitions, tabs, attribute lists, snippets, definitions and
# captions: GitHub and the Docs tab print them as text.
MKDOCS_LEFTOVER = re.compile(r"\{ ?[.#:]|\{ ?target=|^!!! |^\?\?\? |^=== \"|--8<--|^: |^///|markdown>|:material-|:octicons-", re.M)


def check_docs_built() -> list[str]:
    try:
        built = build_docs.build()
    except SystemExit as error:
        return [str(error)]
    stale = [path.as_posix() for path, text in built.items() if not path.exists() or path.read_text(encoding="utf-8") != text]
    problems = [f"{p}: does not match its markdown; run python3 tools/build_docs.py" for p in stale]
    # GitHub Pages serves docs/ from main, and by default runs Jekyll over it
    # first, which makes a page of each markdown file with front matter at the
    # address of the page built from it.
    if not (SITE / ".nojekyll").exists():
        problems.append("docs/.nojekyll is missing: GitHub Pages would run Jekyll over docs/")
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
    # Piece: a placement of the other piece of a gene a contig end cut, which
    # recover_missing_genes must not take for a copy.
    mp = pd.DataFrame(
        [{"Prot Id": r["prot"], "Prot len": r["len"], "Strand": r["st"], "Start": r["s"], "End": r["e"], "Match Nt": r["match"],
          "Num nt": r["num"], "Contig": contig, "Complete": r["complete"], "Piece": bool(r.get("piece", False))} for r in win["ref"]],
        columns=["Prot Id", "Prot len", "Strand", "Start", "End", "Match Nt", "Num nt", "Contig", "Complete", "Piece"],
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
    return problems + check_duplications()


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


def dup_data() -> dict:
    text = DUP_DATA.read_text(encoding="utf-8")
    start = text.index("window.AMAP_DUP_DATA =") + len("window.AMAP_DUP_DATA =")
    return json.loads(text[start:].strip().rstrip(";"))


def write_hits(path: Path, names: list[str], rows: list[list[int]]) -> None:
    """Diamond's table, as GET_DUPLICATIONS writes it, from the hits a drawing
    carries: the columns the scripts read, and placeholders for the rest."""
    lines = ["# DIAMOND, rebuilt by tools/check_site.py"]
    for q, s, nident, qlen, slen in rows:
        lines.append("\t".join(map(str, [names[q], names[s], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qlen, slen, nident])))
    path.write_text("\n".join(lines) + "\n")


def check_duplications() -> list[str]:
    """The family drawn is the one get_duplications.py finds in its hits, with
    the links it keeps; the communities and their numbers are the ones
    network_genomes.py gives."""
    if not DUP_DATA.exists():
        return []
    sys.path.insert(0, str(ROOT / "bin"))
    try:
        import get_duplications
        import network_genomes
    except ImportError as error:
        return [f"assets/dup-data.js: cannot run bin/ ({error}); pip install -r requirements-test.txt"]
    data = dup_data()
    problems = []
    fam = data["family"]
    with tempfile.TemporaryDirectory() as tmp:
        hits = Path(tmp) / "duplications"
        write_hits(hits, fam["names"], fam["rows"])
        kept = get_duplications.read_and_create(str(hits))
        groups = get_duplications.group_related_pairs(kept) if len(kept) else None
    pairs = {tuple(sorted((fam["names"].index(q), fam["names"].index(s)))) for q, s in zip(kept["qseqid"], kept["sseqid"])}
    for p in fam["pairs"]:
        if ((p["a"], p["b"]) in pairs) != p["kept"]:
            problems.append(
                f"assets/dup-data.js: draws {fam['nodes'][p['a']]['id']} and {fam['nodes'][p['b']]['id']} as "
                f"{'a link' if p['kept'] else 'no link'}, get_duplications.py says otherwise"
            )
    members = {fam["names"][i] for i, n in enumerate(fam["nodes"]) if n["member"]}
    found = [set(m) for m in groups["Members"]] if groups is not None else []
    if members not in found:
        problems.append(f"assets/dup-data.js: the family drawn is not one get_duplications.py finds: it finds {found}")

    com = data["communities"]
    with tempfile.TemporaryDirectory() as tmp:
        hits = Path(tmp) / "all_vs_all.diamond"
        write_hits(hits, com["names"], com["rows"])
        with contextlib.redirect_stdout(io.StringIO()):
            df = network_genomes.load_and_filter_blast_results(hits)
    genome_of = dict(zip(com["names"], com["genomes"]))
    networks = network_genomes.build_intra_genome_networks(df, network_genomes.get_genome_ids(df, genome_of), genome_of)
    labelled = network_genomes.detect_and_label_communities(networks)
    homologous = network_genomes.find_homologous_communities(labelled, df)
    for side in com["sides"]:
        name = side["genome"]
        got = labelled.get(name, {}).get(side["community"])
        if got is None or set(got) != set(side["members"]):
            problems.append(f"assets/dup-data.js: network_genomes.py does not split {name} into the {side['community']} drawn")
            continue
        metrics = network_genomes.calculate_subgraph_metrics(networks[name].subgraph(got))
        metrics["size"] = len(got)
        for key, value in side["metrics"].items():
            if abs(float(metrics[key]) - value) > 1e-9:
                problems.append(f"assets/dup-data.js: {name} {side['community']}: {key} is drawn as {value}, network_genomes.py gives {metrics[key]}")
        edges = {tuple(sorted((side["members"].index(a), side["members"].index(b)))) for a, b in networks[name].subgraph(got).edges()}
        if edges != {tuple(e) for e in side["edges"]}:
            problems.append(f"assets/dup-data.js: {name} {side['community']}: the links drawn are not the network's")
    a, b = ((s["genome"], s["community"]) for s in com["sides"])
    if b not in homologous.get(a, []):
        problems.append(f"assets/dup-data.js: network_genomes.py does not pair {a} with {b}")
    one, two = (set(s["members"]) for s in com["sides"])
    cross = sum((q in one and s in two) or (q in two and s in one) for q, s in zip(df["qseqid"], df["sseqid"]))
    if cross != com["cross"]:
        problems.append(f"assets/dup-data.js: {com['cross']} hits are said to join the two communities, the filtered hits have {cross}")
    # The run behind the data searches the way the pipeline does.
    def searches(text: str) -> list[tuple]:
        out = []
        for call in re.findall(r"diamond blastp(.*?)>", text, re.S):
            call = call.replace("\\", " ")
            out.append((tuple(sorted(set(re.findall(r"(?<![\w-])(--?[a-z][\w-]*)", call)))), re.search(r"--outfmt([ \w]+)", call).group(1).split()))
        return out
    pipeline = searches((ROOT / "modules/local/get.nf").read_text(encoding="utf-8"))
    ours = searches((ROOT / "tools/dup_run.sh").read_text(encoding="utf-8"))
    if pipeline != ours:
        problems.append(f"tools/dup_run.sh: runs Diamond as {ours}, modules/local/get.nf as {pipeline}")
    return problems


# ── thresholds ──


def code_value(path: str, pattern: str) -> float | None:
    """The number a pattern captures in a file of bin/, or None."""
    m = re.search(pattern, (ROOT / path).read_text(encoding="utf-8"), re.M)
    return float(m.group(1).replace("_", "")) if m else None


# What each rule's number is in the code. A pattern that stops matching is a
# problem in its own right: the code moved and this table has to follow it.
RULES = {
    "single_bp": ("bin/classify/overlaps.py", r"abs\(bakta_len - rv_len\) > (\d+)"),
    "multiple_bp": ("bin/classify/overlaps.py", r"\(best_len - _range_length\(miniprot_range\)\) < (\d+)"),
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
    "near_max_residues": ("bin/classify/exact.py", r"^NEAR_MAX_RESIDUES = (\d+)"),
    "near_max_differences": ("bin/classify/exact.py", r"^NEAR_MAX_DIFFERENCES = ([\d.]+)"),
    "piece_edge_bp": ("bin/classify/pieces.py", r"^EDGE = (\d+)"),
    "piece_edge_bp_pass": ("bin/contig_end_pieces.py", r"^EDGE = (\d+)"),
    "piece_min_residues": ("bin/contig_end_pieces.py", r"^MIN_RESIDUES = (\d+)"),
    # A placement names a longer call when 2 * span >= 3 * residues: half the codons.
    "name_share": ("bin/classify/reconciliation.py", r"or (\d+) \* span >= 3 \* residues"),
    "copy_cover": ("bin/deduplicate_rvs.py", r"def _covered\(gene, others, share=([\d.]+)\)"),
    "pieces_length": ("bin/deduplicate_rvs.py", r"^MAX_PIECES_LENGTH = ([\d.]+)"),
    "nterm_identity": ("bin/check_cds.py", r"if identity\(ext_aa, ref\[:k\]\) < ([\d.]+):"),
    "dup_cover_pct": ("bin/get_duplications.py", r"\(df_filtered\['p_query'\] >= (\d+)\)"),
    "dup_cover_pct_all": ("bin/network_genomes.py", r"\(df_filtered\['query_coverage'\] >= (\d+)\)"),
}


def check_thresholds() -> list[str]:
    problems = []
    values = {}
    for rule, (path, pattern) in RULES.items():
        value = code_value(path, pattern)
        if value is None:
            problems.append(f"tools/check_site.py: the rule {rule} is no longer where it was in {path}; update RULES")
        values[rule] = value
    # The naming rule is written as a ratio: 2 * span >= 3 * residues.
    if values.get("name_share"):
        values["name_share"] = 1 / values["name_share"]
    if values.get("piece_edge_bp") != values.get("piece_edge_bp_pass"):
        problems.append(f"bin/: a piece is cut within {values.get('piece_edge_bp')} bp of a contig edge in classify/pieces.py, "
                        f"{values.get('piece_edge_bp_pass')} bp in contig_end_pieces.py")
    # What the pages say each number is.
    said = {
        "single_bp": [("index.html", r"longer by more than ([\d,]+) bp")],
        "false_positive_pct": [("index.html", r"genome at (\d+)% identity")],
        "dup_cover_pct": [("how-it-works/duplications.html", r"(\d+)% or more"), ("how-it-works/duplications.html", r"under (\d+)%")],
        "dup_cover_pct_all": [("how-it-works/duplications.html", r"(\d+)% or more")],
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
    # The tables of what fails a run and what only warns, on the inputs page.
    try:
        limits = build_docs.validation_ranges()
    except SystemExit as error:
        problems.append(str(error))
        limits = None
    if limits:
        mb = lambda v: f"{v / 1e6:g} Mb" if v >= 1e6 else f"{v / 1e3:g} kb"
        expected = {
            ("What fails the run", "Genome size"): f"below {mb(limits['size_fail_low'])}, or above {mb(limits['size_fail_high'])}",
            ("What fails the run", "N content"): f"above {limits['n_fail']:g}%",
            ("What fails the run", "GC content"): f"below {limits['gc_fail_low']:g}% or above {limits['gc_fail_high']:g}%",
            ("What only warns", "Genome size"): f"below {mb(limits['size_warn_low'])}, or above {mb(limits['size_warn_high'])}",
            ("What only warns", "N content"): f"above {limits['n_warn']:g}%",
            ("What only warns", "GC content"): f"below {limits['gc_warn_low']:g}% or above {limits['gc_warn_high']:g}%",
            ("What only warns", "Assembly fragmentation"): f"more than {limits['contigs_warn']:g} contigs, and again above {limits['contigs_warn_high']:g}",
        }
        inputs_md = SITE / "getting-started" / "inputs.md"
        heading, seen_rows = None, set()
        for line in inputs_md.read_text(encoding="utf-8").splitlines():
            if line.startswith("### "):
                heading = line[4:].strip()
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if len(cells) == 2 and (heading, cells[0]) in expected:
                seen_rows.add((heading, cells[0]))
                if cells[1] != expected[(heading, cells[0])]:
                    problems.append(f"{inputs_md}: {heading}, {cells[0]}: reads {cells[1]!r}, bin/validate_inputs.py says {expected[(heading, cells[0])]!r}")
        for missing in sorted(set(expected) - seen_rows):
            problems.append(f"{inputs_md}: no row for {missing[1]} under {missing[0]!r} to check")
    # What each step asks for, in the tutorial's table, and the caps the
    # profiles set, where the pages give them.
    config = CONFIG.read_text(encoding="utf-8")
    # The top-level process { } only: the slurm profile gives Bakta a block of
    # its own, with its QoS and nothing it asks for.
    scope = re.search(r"^process \{\n(.*?)^\}", config, re.S | re.M).group(1)
    blocks = dict(re.findall(r"withName: '(\w+)' \{(.*?)\n    \}", scope, re.S))

    def asks(process: str) -> tuple[int, int]:
        body = blocks.get(process, "")
        cpus = re.search(r"cpus\s*=\s*(\d+)", body)
        memory = re.search(r"memory\s*=\s*\{?\s*'?(\d+)\.GB", body)
        return (int(cpus.group(1)) if cpus else None, int(memory.group(1)) if memory else None)

    steps = {
        "eggNOG annotation": ["ANNOT_EGGNOGMAP_BATCH"], "eggNOG database build": ["EGGNOG_DB"], "Pangenome": ["PANGENE"],
        "Duplication analysis": ["GET_DUPLICATIONS", "GET_ALLDUP"], "Report": ["FINAL_REPORT"], "Bakta annotation": ["BAKTA_ANNOTATOR"],
    }
    tutorial = SITE / "getting-started" / "tutorial.md"
    for line in tutorial.read_text(encoding="utf-8").splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) == 3 and cells[0] in steps:
            got = [asks(name) for name in steps[cells[0]]]
            cpus = sorted({c for c, _ in got})
            memory = sorted({m for _, m in got})
            want_cpus = str(cpus[0]) if len(cpus) == 1 else f"{cpus[0]} to {cpus[-1]}"
            want_memory = " or ".join(f"{m} GB" for m in memory)
            if (cells[1], cells[2]) != (want_cpus, want_memory):
                problems.append(f"{tutorial}: {cells[0]} asks for {cells[1]} CPUs and {cells[2]}, nextflow.config says {want_cpus} and {want_memory}")
    # The caps nextflow.config sets after its profiles, where the pages give
    # them: each named group of a pattern is checked against that cap.
    caps = build_docs.caps()
    number = lambda value: int(re.match(r"\d+", value).group())
    running, quickstart, parameters = (SITE / "user-guide" / "running.md", SITE / "getting-started" / "quickstart.md",
                                       SITE / "user-guide" / "parameters.md")
    said_caps = [
        (running, r"(?P<cpus>\d+) CPUs, (?P<memory>\d+) GB and (?P<time>\d+) h with `highmem`", "highmem"),
        (running, r"(?P<cpus>\d+) CPUs, (?P<memory>\d+) GB and (?P<time>\d+) h with `test`", "test"),
        (running, r"(?P<cpus>\d+) CPUs, (?P<memory>\d+) GB and (?P<time>\d+) h with an executor profile", "scheduler"),
        (running, r"otherwise the machine's CPUs, (?P<memory>\d+) GB and (?P<time>\d+) h", "machine"),
        (quickstart, r"the machine's CPUs and (?P<memory>\d+) GB", "machine"),
        (quickstart, r"scheduler profiles at (?P<cpus>\d+) CPUs and (?P<memory>\d+) GB", "scheduler"),
        (parameters, r"\| `--max_memory` \|[^|]*\|[^|]*\| `(?P<memory>\d+)\.GB` \|", "machine"),
        (parameters, r"\| `--max_time` \|[^|]*\|[^|]*\| `(?P<time>\d+)\.h` \|", "machine"),
        (parameters, r"HPC profiles raise them to (?P<cpus>\d+) CPUs and (?P<memory>\d+) GB", "scheduler"),
        (parameters, r"`-profile highmem` to (?P<cpus>\d+) CPUs, (?P<memory>\d+) GB and (?P<time>\d+) h", "highmem"),
        (parameters, r"`-profile test` sets (?P<cpus>\d+) CPUs, (?P<memory>\d+) GB and (?P<time>\d+) h", "test"),
    ]
    for page, pattern, kind in said_caps:
        m = re.search(pattern, re.sub(r"\s+", " ", page.read_text(encoding="utf-8")))
        if not m:
            problems.append(f"{page}: no longer says {pattern!r}; update tools/check_site.py")
            continue
        for cap, said in m.groupdict().items():
            if int(said) != number(caps[kind][cap]):
                problems.append(f"{page}: says {m.group(0)!r}, nextflow.config gives {caps[kind][cap]} ({kind})")
    # The versions the containers page gives are the image's.
    try:
        image = build_docs.image_contents()
    except SystemExit as error:
        problems.append(str(error))
        image = None
    if image:
        page = SITE / "deployment" / "containers.md"
        text = page.read_text(encoding="utf-8")
        for name, version in image["tools"] + image["runtimes"]:
            if not version:
                continue
            said = re.search(rf"\b{re.escape(name)} ([\d.]+)\b", text, re.I)
            if not said:
                problems.append(f"{page}: does not give {name}'s version ({version} in docker/Dockerfile)")
            elif said.group(1).rstrip(".") != version:
                problems.append(f"{page}: says {name} {said.group(1)}, docker/Dockerfile installs {version}")
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
        ("Longest protein looked for with substitutions", "small proteins"): ("near_max_residues", "{:g} residues"),
        ("Residues after the first that may differ", "small proteins"): ("near_max_differences", "{:g}%", 0.01),
        ("Distance from a contig's edge for a placement, or a piece, to be cut there", "pieces of cut genes"): ("piece_edge_bp", "{:g} bp"),
        ("Shortest part of a protein looked for at the other contig ends", "pieces of cut genes"): ("piece_min_residues", "{:g} residues"),
        ("Share of its protein an alignment spans to name a longer call", "rescuing missing genes"): ("name_share", "{:g}%", 0.01),
        ("Cover by another reference gene that rules a copy out", "final deduplication"): ("copy_cover", "{:g}%", 0.01),
        ("Length of a gene's pieces together, next to the reference gene", "final deduplication"): ("pieces_length", "{:g}%", 0.01),
        ("Identity to restore a clipped N-terminus", "quality control"): ("nterm_identity", "{:g}%", 0.01),
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


# ── benchmark ──


BENCHMARK = Path("tests/benchmark/README.md")


def check_benchmark(pages: dict[str, Page]) -> list[str]:
    """The front page's scores are the benchmark's expected results, row by row."""
    text = BENCHMARK.read_text(encoding="utf-8")
    table = text[text.index("## Expected results"):]
    expected = []
    for line in table.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) == 5 and re.match(r"\d+ / \d+", cells[1]):
            pair = lambda c: tuple(int(n) for n in c.split(" / ")) if c else None
            expected.append((cells[0], pair(cells[1]), pair(cells[2])))
    front = (SITE / "index.html").read_text(encoding="utf-8")
    rows = re.findall(r'<div class="pf-row"><span class="pf-name"><b>([^<]+)</b>.*?</div>', front, re.S)
    bars = re.findall(r'<div class="pf-row"><span class="pf-name">.*?</div>', front, re.S)
    problems = []
    if len(bars) != len(expected):
        return [f"index.html: shows {len(bars)} genomes of the benchmark, {BENCHMARK} has {len(expected)}"]
    for (genome, exact, changed), row, label in zip(expected, bars, rows):
        got = {cls: (int(f), int(n)) for cls, f, n in re.findall(r'class="pf-bar (pf-\w+)" data-found="(\d+)" data-total="(\d+)"', row)}
        if got.get("pf-exact") != exact or got.get("pf-changed") != changed:
            problems.append(f"index.html: {label} shows {got}, {BENCHMARK} expects exact {exact} and changed {changed} ({genome})")
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
    # Steps that only move files for another: no station of their own.
    helpers = {"SOFTWARE_VERSIONS", "DUMP_PARAMS", "TEST_GENOME", "CREATE_EMPTY_FILE", "REPORT_INPUTS", "PANGENE_GENOME"}
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
    # Its pictures are files of the repository.
    for src in re.findall(r'(?:src|srcset)="([^"]+)"', text) + re.findall(r"!\[[^\]]*\]\(([^)\s]+)", text):
        if not re.match(r"^https?:", src) and not Path(src).exists():
            problems.append(f"README.md: shows {src}, which is not in the repository")
    # GitHub takes a <picture> out of a link and then shows one picture in
    # both themes, so none may sit inside an <a>.
    for anchor in re.findall(r"<a\b[^>]*>(.*?)</a>", text, re.S):
        if "<picture" in anchor:
            problems.append("README.md: a <picture> inside a link loses its dark version on GitHub")
    # The benchmark's scores, as tests/benchmark/README.md expects them.
    bench = BENCHMARK.read_text(encoding="utf-8")
    pair = lambda c: tuple(int(n) for n in c.split(" / ")) if c.strip() else None
    expected = [
        (pair(cells[1]), pair(cells[2]))
        for cells in ([c.strip() for c in line.strip().strip("|").split("|")] for line in bench[bench.index("## Expected results"):].splitlines())
        if len(cells) == 5 and re.match(r"\d+ / \d+", cells[1])
    ]
    shown = [
        (pair(cells[1]), pair(cells[2]))
        for cells in ([c.strip() for c in line.strip().strip("|").split("|")] for line in text[text.index("## Does it work?"):].splitlines())
        if len(cells) == 3 and re.match(r"\d+ / \d+", cells[1])
    ] if "## Does it work?" in text else []
    if shown and shown != expected:
        problems.append(f"README.md: the benchmark table shows {shown}, {BENCHMARK} expects {expected}")
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
    files += [README, Path("CHANGELOG.md")] + sorted(Path("tools").glob("*.*"))
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
        ("benchmark", lambda: check_benchmark(pages)),
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
