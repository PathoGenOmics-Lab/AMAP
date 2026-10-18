#!/usr/bin/env python3
"""Build the website's documentation pages from the markdown in docs/.

The markdown stays the source: GitHub renders it, and this writes an .html page
next to each file with the site's header, theme, sidebar, page contents and
search, plus docs/documentation.html (the Docs tab) and the search index. The
pages are committed, so the site opens from a clone with no server and no build
step; --check fails when they no longer match the markdown, which is what CI
runs.

    python3 tools/build_docs.py           # write the pages
    python3 tools/build_docs.py --check   # exit 1 if any page is stale

Only the markdown these pages use is supported: headings, paragraphs, lists
(nested, loose or tight, with code inside), tables, block quotes and GitHub's
alerts (> [!NOTE]), fenced code (`bash title="run.sh"` gives the block a
caption, and a run of **Label** paragraphs each followed by a code block
becomes tabs), rules, inline code, emphasis, links, a few inline tags and raw
HTML blocks. A comment <!-- figure: NAME --> on a line of its own puts that
drawing from assets/figures.js there; GitHub shows nothing in its place.
Anything else shows up as text, and tools/check_site.py catches what that
breaks.
"""

from __future__ import annotations

import argparse
import html
import json
import re
import sys
from pathlib import Path

DOCS = Path("docs")
REPO = "https://github.com/PathoGenOmics-Lab/AMAP"

# The Docs tab, in reading order: (section, path, title, one line). A path to
# a .md file is built into a page; an .html path is a page written by hand and
# only gets its place in the list.
PAGES = [
    ("Getting started", "getting-started/tutorial.md", "Tutorial", "From an empty folder to a corrected genome, with the numbers of a real run."),
    ("Getting started", "getting-started/installation.md", "Installation", "Nextflow, a container engine, and running without one."),
    ("Getting started", "getting-started/quickstart.md", "Quickstart", "From an installed Nextflow to an open report."),
    ("Getting started", "getting-started/inputs.md", "Preparing inputs", "The sample sheet, the genomes, and the checks run on them."),
    ("User guide", "user-guide/parameters.md", "Parameters", "Every option, its default, and when to change it."),
    ("User guide", "user-guide/running.md", "Running AMAP", "Profiles, schedulers and larger runs."),
    ("User guide", "user-guide/outputs.md", "Outputs", "Every folder and file, and what a healthy run looks like."),
    ("User guide", "user-guide/report.md", "The report", "What each section of final_report.html shows."),
    ("User guide", "user-guide/custom-reference.md", "Custom reference", "Any organism, against a trusted annotation of your own."),
    ("How it works", "how-it-works.html", "Step by step, illustrated", "One gene at a time, from Bakta's call to what AMAP keeps."),
    ("How it works", "how-it-works/index.md", "Overview", "What goes into a decision, and the five outcomes."),
    ("How it works", "how-it-works/workflow.md", "The workflow", "Each step of the pipeline, on the metro map."),
    ("How it works", "how-it-works/correction-algorithm.md", "The correction algorithm", "Every rule and threshold, drawn on real genes."),
    ("How it works", "how-it-works/duplications.md", "Duplication metrics", "What the network numbers mean."),
    ("Deployment", "deployment/containers.md", "Containers", "The AMAP image, what is in it, and building your own."),
    ("Deployment", "deployment/hpc.md", "HPC clusters", "Schedulers, scratch, and where things should live."),
    ("Deployment", "deployment/troubleshooting.md", "Troubleshooting", "Organised by the error you are looking at."),
    ("About", "about/changelog.md", "Changelog", "What changed in each version."),
    ("About", "about/citation.md", "Citation", "How to cite AMAP and the tools it runs."),
    ("About", "about/contributing.md", "Contributing", "Working on the code, the tests and these pages."),
]
INDEX = "documentation.html"
SEARCH_INDEX = DOCS / "assets" / "search-index.js"
# The changelog page is the repository's CHANGELOG.md, so there is one source.
CHANGELOG = ("about/changelog.md", Path("CHANGELOG.md"))

GH_PATH = (
    "M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94"
    "-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07"
    "-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0"
    " 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25"
    ".54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"
)
GH_SVG = f'<svg viewBox="0 0 16 16" fill="currentColor" aria-hidden="true"><path d="{GH_PATH}"/></svg>'
SUN_SVG = (
    '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" aria-hidden="true">'
    '<circle cx="12" cy="12" r="4.5"/><path d="M12 2v2M12 20v2M4.9 4.9l1.4 1.4M17.7 17.7l1.4 1.4M2 12h2M20 12h2M4.9 19.1'
    'l1.4-1.4M17.7 6.3l1.4-1.4"/></svg>'
)
THEME_BOOT = (
    '<script>try { document.documentElement.dataset.theme = localStorage.getItem("amap.theme") === "light" ? '
    '"light" : "dark"; } catch (e) { document.documentElement.dataset.theme = "dark"; }</script>'
)
INLINE_TAGS = re.compile(r"</?(kbd|br|sup|sub|strong|em|b|i)\s*/?>", re.I)
ALERTS = {"NOTE": "Note", "TIP": "Tip", "IMPORTANT": "Important", "WARNING": "Warning", "CAUTION": "Caution"}
FIGURE_MARK = re.compile(r"^<!--\s*figure:\s*([\w-]+)\s*-->$")


# ── Figures ──
#
# What a <!-- figure: NAME --> becomes. The drawings themselves are in
# assets/figures.js, which every page with one of these loads.


def case_figure(name: str, caption: str) -> str:
    return (
        f'<figure class="fig case" data-fig="{name}"><div class="case-fig" data-fig-stage></div>'
        f"<figcaption>{caption}</figcaption></figure>"
    )


FIGURES = {
    "small-protein": case_figure("small-protein", "E. coli K-12 against its own annotation. miniprot placed IbsC (b4665), 19 residues, only on its neighbour IbsD's locus; the exact search finds it here, between fau and serA, where Bakta called nothing."),
    "paralogs-1": case_figure("paralogs-1", "CDC1551 near 400 kb: Rv0336's own locus. Its paralog Rv0515 is placed here too, and that placement is dropped."),
    "paralogs-2": case_figure("paralogs-2", "Near 608 kb: Rv0515's own locus, where Rv0336's placement is dropped."),
    "copies-1": case_figure("copies-1", "E. coli O157:H7 Sakai against K-12: the first of two whole copies of RhsD (b0497)."),
    "copies-2": case_figure("copies-2", "The second copy, 1.4 Mb away, is a further copy and stays, as b0497_2."),
    "frameshift": case_figure("frameshift", "M. bovis: Rv0590's frame runs on 67 residues past its protein, through Rv0590A, to the end of Bakta's one long call. It keeps its extent as a pseudogene, and Rv0590A is a gene of its own."),
    "contig-end": case_figure("contig-end", "CDC1551 in 40 contigs: Rv1239c's frame reaches the end of contig_10 with no stop codon, and the gene is marked partial at its 3' end."),
    "case2": case_figure("case2", "Case 2 in CDC1551: dnaN, where Bakta and the aligned protein agree base for base."),
    "case3": case_figure("case3", "Case 3: fhaA, with the same ends as Rv0020c but six residues fewer."),
    "case4-start": case_figure("case4-start", "A single overlap: Bakta's dnaA starts 33 bp late, and the reference's start is kept."),
    "case4-split": case_figure("case4-split", "Multiple overlaps: two Bakta calls on Rv0007, which comes back as one gene."),
    "case4-longer": case_figure("case4-longer", "A much longer call: 3,303 bp past an alignment that stops short of its protein's end, so Bakta's call is kept."),
    "case1-kept": case_figure("case1-kept", "Case 1 kept: four calls in a stretch CDC1551 has and H37Rv lacks."),
    "case1-dropped": case_figure("case1-dropped", "Case 1 removed: a call that Bakta's known false positives cover end to end."),
    "missing": case_figure("missing", "A missing reference gene: Rv0031, put back between two of Bakta's calls."),
    "metro": (
        '<figure class="fig metro-wrap" data-metro><div class="metro-opts"><span>Run with</span>'
        '<label><input type="checkbox" data-metro-opt="skipTyping"><code>--skip_typing</code></label>'
        '<label><input type="checkbox" data-metro-opt="skipFunctional"><code>--skip_functional_annotation</code></label>'
        '<label><input type="checkbox" data-metro-opt="custom">your own reference</label></div>'
        '<p class="metro-hint">Scroll sideways to follow the lines.</p><div class="metro-scroll" data-metro-stage></div><div class="metro-info" data-metro-info aria-live="polite"></div>'
        '<p class="metro-writes" data-metro-writes></p></figure>'
    ),
}


# Figures drawn from numbers the page itself gives, so the two cannot differ.

CASE_PARTS = [
    ("case2", "case 2", "agree"),
    ("case3", "case 3", "agree3"),
    ("case4", "case 4", "conflict"),
    ("case1_removed", "case 1, dropped", "dropped"),
    ("case1_kept", "case 1, kept", "novel"),
]


def case_bars(genomes: list[tuple[str, dict[str, int]]], caption: str) -> str:
    """Bakta's calls of each genome by case, and the reference genes put back."""
    longest = max(sum(n.get(k, 0) for k, _, _ in CASE_PARTS) + n.get("rv_missing", 0) for _, n in genomes)
    rows = []
    for name, n in genomes:
        segs = []
        for key, label, cls in CASE_PARTS:
            if key in n:
                width = 100 * n[key] / longest
                number = f'<span class="cb-n">{n[key]:,}</span>' if width >= 8 else ""
                segs.append(f'<span class="cb-seg cb-{cls}" style="width:{width:.2f}%" title="{label}: {n[key]:,}">{number}</span>')
        if "rv_missing" in n:
            width = 100 * n["rv_missing"] / longest
            number = f'<span class="cb-n">{n["rv_missing"]:,}</span>' if width >= 8 else ""
            segs.append(f'<span class="cb-seg cb-missing" style="width:{width:.2f}%" title="reference genes put back: {n["rv_missing"]:,}">{number}</span>')
        total = f'<span class="cb-total">{n["total_cds"]:,} Bakta calls</span>' if "total_cds" in n else ""
        values = "".join(
            f'<span><i class="cb-key cb-{cls}"></i>{label} <b>{n[key]:,}</b></span>'
            for key, label, cls in CASE_PARTS + [("rv_missing", "missing, put back", "missing")]
            if key in n
        )
        rows.append(
            f'<div class="cb-row"><span class="cb-name">{html.escape(name)}</span>{total}<div class="cb-bar" aria-hidden="true">{"".join(segs)}</div>'
            f'<div class="cb-values">{values}</div></div>'
        )
    return (
        f'<figure class="fig case-bars"><div class="cb-rows">{"".join(rows)}</div>'
        f"<figcaption>{caption}</figcaption></figure>"
    )


def results_map() -> str:
    """Each folder of results/, the processes that write it, and what leaves it out."""
    rows = [
        ("main", "PATHOTYPR/", "PATHOTYPR", "with <code>--skip_typing</code>"),
        ("main", "BAKTA/", "BAKTA_ANNOTATOR", ""),
        ("main", "CDS_CASES/", "TYPE_CASES", ""),
        ("main", "AMAP/", "FIX_CDS, then GET_OUTPUTS and GET_COUNT",
         "with <code>--skip_functional_annotation</code> only FIX_CDS's <code>*.AMAP.tsv</code>, <code>.faa</code> and <code>.fna</code> are left"),
        ("dup", "DUPLICATIONS/", "GET_DUPLICATIONS, GET_ALLDUP", "with <code>--skip_functional_annotation</code>"),
        ("pan", "PANGENE/", "PANGENE",
         "with <code>--skip_functional_annotation</code>; with <code>--skip_typing</code> its files are empty"),
        ("main", "final_report.html", "FINAL_REPORT", "with <code>--skip_functional_annotation</code>"),
        ("none", "pipeline_info/", "VALIDATE_INPUTS, SOFTWARE_VERSIONS, DUMP_PARAMS and Nextflow itself", ""),
    ]
    out = []
    for line, folder, by, when in rows:
        gone = f'<span class="rm-when">Left out {when}.</span>' if when else '<span class="rm-when">Always written.</span>'
        out.append(
            f'<div class="rm-row"><i class="rm-line rm-{line}"></i><code class="rm-folder">{folder}</code>'
            f'<div class="rm-what"><span class="rm-by">{by}</span>{gone}</div></div>'
        )
    return (
        '<figure class="fig results-map">'
        + "".join(out)
        + '<figcaption>Each folder under <code>results/</code>, the processes of the <a href="../how-it-works/workflow.html">metro map</a> that write it, '
        "and the options that leave it out. The colour is the metro line of its process.</figcaption></figure>"
    )


EXCERPT = (
    '<figure class="fig excerpt"><div class="excerpt-head"><div class="tabs" role="tablist" aria-label="What the files look like">'
    '<button class="tab" role="tab" aria-selected="true" aria-controls="ex-gff" id="ex-gff-tab">GFF3</button>'
    '<button class="tab" role="tab" aria-selected="false" aria-controls="ex-gbk" id="ex-gbk-tab">GenBank</button></div>'
    "<p>Real lines: CDC1551 cut into 40 contigs.</p></div>"
    '<div class="tabpanel" role="tabpanel" id="ex-gff" aria-labelledby="ex-gff-tab"><pre>'
    "contig_4   AMAP  CDS  164980  165406  .  -  0  ID=L4_CDC1551d_Rv0354c;&hellip;;Dbxref=AMAP_ids:Rv0354c;"
    "<mark>pseudo=true;Note=frameshift or lost stop codon: its frame runs on 41 residues past the end of the reference protein</mark>\n"
    "contig_10  AMAP  CDS  101478  102146  .  +  0  ID=L4_CDC1551d_Rv1239c;&hellip;;Dbxref=AMAP_ids:Rv1239c;"
    '<mark class="part">partial=true;end_range=102146,.</mark></pre>'
    '<div class="key"><span><mark>pseudo=true</mark> a pseudogene, the reason in <code>Note</code></span>'
    '<span><mark class="part">partial=true</mark> cut by its contig\'s end, the open end in <code>end_range</code> or <code>start_range</code></span></div></div>'
    '<div class="tabpanel" role="tabpanel" id="ex-gbk" aria-labelledby="ex-gbk-tab" hidden><pre>'
    "     CDS             complement(164980..165406)\n"
    '                     /locus_tag="L4_CDC1551d_Rv0354c"\n'
    '                     /db_xref="AMAP_ids:Rv0354c"\n'
    '                     <mark>/pseudo=""</mark>\n'
    '                     <mark>/note="frameshift or lost stop codon: its frame runs on 41\n'
    '                     residues past the end of the reference protein"</mark>\n'
    '     CDS             <mark class="part">101478..&gt;102146</mark>\n'
    '                     /locus_tag="L4_CDC1551d_Rv1239c"</pre>'
    '<div class="key"><span><mark>/pseudo</mark> no <code>/translation</code>, the reason in <code>/note</code></span>'
    '<span><mark class="part">&gt;102146</mark> the 3\' end runs off the contig (<code>&lt;</code> marks a cut 5\' end)</span></div></div>'
    "<figcaption><code>&hellip;</code> stands for the attributes left out; the marks are written the way NCBI writes them.</figcaption></figure>"
)


def prepare_reference(numbers: dict[str, int]) -> str:
    """What PREPARE_REFERENCE does with a reference, with the page's numbers."""
    steps = [
        ("1", "Extract the proteome", "Every CDS of your annotation becomes a protein named <code>ID_START-END</code>.",
         f'{numbers["proteins"]:,} proteins'),
        ("2", "Annotate the reference with Bakta", "The same Bakta that will annotate your samples, on a genome whose answer is known.",
         f'{numbers["bakta"]:,} Bakta calls'),
        ("3", "Keep what Bakta invents", "Every Bakta call that overlaps no curated gene: the calls it is prone to make on your samples too.",
         f'{numbers["false"]:,} false positives'),
    ]
    out = []
    for n, title, text, count in steps:
        out.append(
            f'<div class="pr-step"><span class="pr-n">{n}</span><b>{title}</b><p>{text}</p><span class="pr-count">{count}</span></div>'
        )
    return (
        '<figure class="fig prepare-ref"><div class="pr-in"><code>--ref_annotation</code> + <code>--ref_genome</code></div>'
        f'<div class="pr-steps">{"".join(out)}</div>'
        '<div class="pr-out"><span>The proteome is aligned onto every sample; the false positives rule out case 1 calls.</span></div>'
        "<figcaption>What happens to E. coli K-12 MG1655 given as the reference, once, before any sample is touched.</figcaption></figure>"
    )


# ── Markdown ──


def href(path: str) -> str:
    """Where a PAGES entry lives on the site, relative to docs/."""
    return path[: -len(".md")] + ".html" if path.endswith(".md") else path


def slugify(text: str, seen: dict[str, int]) -> str:
    """GitHub's anchor for a heading, so links written for GitHub land here too."""
    slug = re.sub(r"[^\w\- ]", "", text.strip().lower()).replace(" ", "-")
    if slug in seen:
        seen[slug] += 1
        slug = f"{slug}-{seen[slug]}"
    else:
        seen[slug] = 0
    return slug


def plain(markdown_inline: str) -> str:
    """The text of an inline span, without its markup."""
    text = re.sub(r"`([^`]*)`", r"\1", markdown_inline)
    text = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", text)
    text = INLINE_TAGS.sub("", text)
    # The emphasis inline() renders, and only that: an underscore inside a word
    # (skip_functional_annotation) is part of the word, and of its anchor.
    text = re.sub(r"\*\*(.+?)\*\*", r"\1", text)
    text = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"\1", text)
    return re.sub(r"(?<!\w)_(?!\s)(.+?)(?<!\s)_(?!\w)", r"\1", text)


def strip_front_matter(lines: list[str]) -> list[str]:
    if lines and lines[0].strip() == "---":
        for i in range(1, len(lines)):
            if lines[i].strip() == "---":
                return lines[i + 1 :]
    return lines


class Renderer:
    def __init__(self, source: str):
        self.source = source  # markdown path relative to docs/, e.g. user-guide/outputs.md
        self.seen: dict[str, int] = {}
        self.headings: list[tuple[int, str, str]] = []  # (level, text, slug)
        self.sections: list[dict] = []  # for the search index
        self.figures: list[str] = []

    # Links

    def link(self, target: str) -> str:
        if re.match(r"^(https?:|mailto:|#)", target):
            return target
        path, _, anchor = target.partition("#")
        resolved = (DOCS / Path(self.source).parent / path).resolve()
        try:
            inside = resolved.relative_to(DOCS.resolve())
        except ValueError:
            inside = None
        if inside is not None and inside.suffix == ".md" and inside.as_posix() in {p[1] for p in PAGES}:
            new = path[: -len(".md")] + ".html"
        elif inside is not None and inside.suffix != ".md":
            new = path
        else:
            repo_path = resolved.relative_to(Path.cwd().resolve())
            new = f"{REPO}/blob/main/{repo_path.as_posix()}"
        return new + (f"#{anchor}" if anchor else "")

    # Inline

    def inline(self, text: str) -> str:
        keep: list[str] = []

        def stash(fragment: str) -> str:
            keep.append(fragment)
            return f"\x00{len(keep) - 1}\x00"

        text = re.sub(r"(`+)(.+?)\1", lambda m: stash(f"<code>{html.escape(m.group(2).strip())}</code>"), text)
        text = re.sub(r"\\([\\`*_{}\[\]()#+\-.!|<>])", lambda m: stash(html.escape(m.group(1))), text)
        text = re.sub(r"<(https?://[^>\s]+)>", lambda m: stash(f'<a href="{html.escape(m.group(1))}">{html.escape(m.group(1))}</a>'), text)
        text = INLINE_TAGS.sub(lambda m: stash(m.group(0)), text)

        def link(m: re.Match) -> str:
            label, target = m.group(1), m.group(2).strip().split(" ")[0]
            return stash(f'<a href="{html.escape(self.link(target))}">') + label + stash("</a>")

        text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", link, text)
        text = re.sub(r"&(?!#?\w+;)", "&amp;", text).replace("<", "&lt;").replace(">", "&gt;")
        text = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", text)
        text = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"<em>\1</em>", text)
        text = re.sub(r"(?<![\w])_(?!\s)(.+?)(?<!\s)_(?![\w])", r"<em>\1</em>", text)
        while "\x00" in text:
            text = re.sub(r"\x00(\d+)\x00", lambda m: keep[int(m.group(1))], text)
        return text

    # Blocks

    def blocks(self, lines: list[str]) -> str:
        out: list[str] = []
        i = 0
        while i < len(lines):
            line = lines[i]
            stripped = line.strip()
            if not stripped:
                i += 1
                continue
            figure = FIGURE_MARK.match(stripped)
            if figure:
                name = figure.group(1)
                if name in DERIVED:
                    out.append(DERIVED[name](self, lines))
                elif name in FIGURES:
                    self.figures.append(name)
                    out.append(FIGURES[name])
                else:
                    raise SystemExit(f"docs/{self.source}: no figure called {name!r}")
                i += 1
                continue
            if stripped.startswith("<!--"):
                # A comment for whoever edits the markdown, not for the page.
                while i < len(lines) and "-->" not in lines[i]:
                    i += 1
                i += 1
                continue
            tabs = self.tab_group(lines, i)
            if tabs:
                i, rendered = tabs
                out.append(rendered)
                continue
            fence = re.match(r"^(\s*)(```+|~~~+)\s*(.*)$", line)
            if fence:
                indent, marker, info = len(fence.group(1)), fence.group(2), fence.group(3).strip()
                body = []
                i += 1
                while i < len(lines) and not lines[i].strip().startswith(marker):
                    body.append(lines[i][indent:] if lines[i][:indent].strip() == "" else lines[i])
                    i += 1
                i += 1
                out.append(self.code(info, body))
                continue
            heading = re.match(r"^(#{1,6})\s+(.*?)\s*#*$", line)
            if heading:
                level, text = len(heading.group(1)), heading.group(2)
                slug = slugify(plain(text), self.seen)
                self.headings.append((level, plain(text), slug))
                self.sections.append({"heading": plain(text), "slug": slug, "level": level, "text": []})
                anchor = f'<a class="anchor" href="#{slug}" aria-label="Link to this section">#</a>' if level > 1 else ""
                out.append(f'<h{level} id="{slug}">{self.inline(text)}{anchor}</h{level}>')
                i += 1
                continue
            if re.match(r"^(-{3,}|\*{3,}|_{3,})$", stripped):
                out.append("<hr>")
                i += 1
                continue
            if stripped.startswith("|") and i + 1 < len(lines) and re.match(r"^\|?\s*:?-{2,}", lines[i + 1].strip()):
                rows = []
                while i < len(lines) and lines[i].strip().startswith("|"):
                    rows.append(lines[i].strip())
                    i += 1
                out.append(self.table(rows))
                continue
            if stripped.startswith(">"):
                quoted = []
                while i < len(lines) and lines[i].strip().startswith(">"):
                    quoted.append(re.sub(r"^\s*> ?", "", lines[i]))
                    i += 1
                out.append(self.quote(quoted))
                continue
            if re.match(r"^<(picture|div|details|table|p|figure)\b", stripped):
                chunk = []
                while i < len(lines) and lines[i].strip():
                    chunk.append(lines[i])
                    i += 1
                out.append(self.raw_html("\n".join(chunk)))
                continue
            item = re.match(r"^(\s*)([-*+]|\d+\.)\s+", line)
            if item:
                i, rendered = self.list_block(lines, i)
                out.append(rendered)
                continue
            para = []
            while i < len(lines) and lines[i].strip() and not self.starts_block(lines, i):
                para.append(lines[i].strip())
                i += 1
            text = " ".join(para)
            self.note(text)
            out.append(f"<p>{self.inline(text)}</p>")
        return "\n".join(out)

    def starts_block(self, lines: list[str], i: int) -> bool:
        line = lines[i]
        s = line.strip()
        return bool(
            re.match(r"^(\s*)(```|~~~)", line)
            or re.match(r"^#{1,6}\s", line)
            or re.match(r"^(\s*)([-*+]|\d+\.)\s+", line)
            or s.startswith(">")
            or FIGURE_MARK.match(s)
            or (s.startswith("|") and i + 1 < len(lines) and re.match(r"^\|?\s*:?-{2,}", lines[i + 1].strip()))
        )

    def note(self, text: str) -> None:
        if self.sections:
            self.sections[-1]["text"].append(plain(text))

    def quote(self, quoted: list[str]) -> str:
        """A block quote, or one of GitHub's alerts (> [!WARNING]) as a callout
        whose first bold line, if it has one, is its title."""
        first = next((q.strip() for q in quoted if q.strip()), "")
        alert = re.match(r"^\[!(\w+)\]$", first)
        if not alert or alert.group(1).upper() not in ALERTS:
            return f"<blockquote>{self.blocks(quoted)}</blockquote>"
        kind = alert.group(1).upper()
        body = quoted[quoted.index(next(q for q in quoted if q.strip())) + 1 :]
        while body and not body[0].strip():
            body = body[1:]
        title = ALERTS[kind]
        if body and re.match(r"^\*\*(.+)\*\*$", body[0].strip()):
            title = self.inline(body[0].strip()[2:-2])
            self.note(plain(body[0].strip()))
            body = body[1:]
        return (
            f'<div class="callout callout-{kind.lower()}"><p class="callout-title">{title}</p>'
            f"{self.blocks(body)}</div>"
        )

    def tab_group(self, lines: list[str], i: int) -> tuple[int, str] | None:
        """Two or more **Label** paragraphs, each followed by a code block, as
        tabs: GitHub shows the labelled blocks one after the other."""
        tabs: list[tuple[str, str, list[str]]] = []
        j = i
        while j < len(lines):
            label = re.match(r"^\*\*([^*]+)\*\*$", lines[j].strip())
            if not label:
                break
            k = j + 1
            while k < len(lines) and not lines[k].strip():
                k += 1
            fence = re.match(r"^(```+|~~~+)\s*(.*)$", lines[k].strip()) if k < len(lines) else None
            if not fence:
                break
            body = []
            k += 1
            while k < len(lines) and not lines[k].strip().startswith(fence.group(1)):
                body.append(lines[k])
                k += 1
            tabs.append((label.group(1), fence.group(2).strip(), body))
            j = k + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
        if len(tabs) < 2:
            return None
        group = f"tabs-{re.sub(r'[^a-z0-9]+', '-', self.source.lower())}-{i}"
        buttons, panels = [], []
        for n, (label, info, body) in enumerate(tabs):
            tid = f"{group}-{n}"
            selected = "true" if n == 0 else "false"
            buttons.append(
                f'<button class="tab" role="tab" aria-selected="{selected}" aria-controls="{tid}" id="{tid}-tab">{html.escape(label)}</button>'
            )
            hidden = "" if n == 0 else " hidden"
            panels.append(f'<div class="tabpanel" role="tabpanel" id="{tid}" aria-labelledby="{tid}-tab"{hidden}>{self.code(info, body)}</div>')
            self.note(label)
        return j, f'<div class="tabset"><div class="tabs" role="tablist">{"".join(buttons)}</div>{"".join(panels)}</div>'

    def list_block(self, lines: list[str], i: int) -> tuple[int, str]:
        first = re.match(r"^(\s*)([-*+]|\d+\.)\s+", lines[i])
        base = len(first.group(1))
        ordered = first.group(2)[0].isdigit()
        start = int(first.group(2)[:-1]) if ordered else 1
        items: list[list[str]] = []
        loose = False
        while i < len(lines):
            line = lines[i]
            m = re.match(r"^(\s*)([-*+]|\d+\.)(\s+)", line)
            if m and len(m.group(1)) == base and (m.group(2)[0].isdigit()) == ordered:
                content_indent = len(m.group(1)) + len(m.group(2)) + len(m.group(3))
                items.append([line[content_indent:]])
                i += 1
                while i < len(lines):
                    nxt = lines[i]
                    if not nxt.strip():
                        # A blank line ends the list unless the next line carries on in it.
                        j = i + 1
                        while j < len(lines) and not lines[j].strip():
                            j += 1
                        if j < len(lines) and (len(lines[j]) - len(lines[j].lstrip())) > base:
                            items[-1].append("")
                            loose = loose or not re.match(r"^\s*([-*+]|\d+\.)\s", lines[j])
                            i += 1
                            continue
                        if j < len(lines) and re.match(rf"^\s{{{base}}}([-*+]|\d+\.)\s", lines[j]):
                            loose = True
                            i = j
                        break
                    indent = len(nxt) - len(nxt.lstrip())
                    if indent > base:
                        items[-1].append(nxt[min(indent, content_indent):])
                        i += 1
                        continue
                    if re.match(r"^(\s*)([-*+]|\d+\.)\s+", nxt):
                        break
                    items[-1].append(nxt.strip())  # lazy continuation
                    i += 1
                continue
            break
        tag = "ol" if ordered else "ul"
        attrs = f' start="{start}"' if ordered and start != 1 else ""
        rendered = []
        for content in items:
            body = self.blocks(content)
            if not loose and body.startswith("<p>") and body.count("<p>") == 1:
                body = body[3:].replace("</p>", "", 1)
            rendered.append(f"<li>{body}</li>")
        return i, f"<{tag}{attrs}>" + "".join(rendered) + f"</{tag}>"

    def table(self, rows: list[str]) -> str:
        def cells(row: str) -> list[str]:
            row = row.strip().strip("|")
            return [c.strip().replace("\\|", "|") for c in re.split(r"(?<!\\)\|", row)]

        head = cells(rows[0])
        aligns = cells(rows[1])
        body = [cells(r) for r in rows[2:]]
        self.last_table = (head, body)
        for row in body:
            self.note(" ".join(row))

        def cls(k: int) -> str:
            return ' class="num"' if k < len(aligns) and aligns[k].endswith(":") and not aligns[k].startswith(":") else ""

        out = ['<div class="table-wrap"><table><thead><tr>']
        out += [f"<th{cls(k)}>{self.inline(c)}</th>" for k, c in enumerate(head)]
        out.append("</tr></thead><tbody>")
        for row in body:
            out.append("<tr>" + "".join(f"<td{cls(k)}>{self.inline(c)}</td>" for k, c in enumerate(row)) + "</tr>")
        out.append("</tbody></table></div>")
        return "".join(out)

    def code(self, info: str, body: list[str]) -> str:
        self.last_code = body
        title = re.search(r'title="([^"]*)"', info)
        lang = info.split()[0] if info and not info.startswith("title=") else ""
        caption = f'<div class="code-title">{html.escape(title.group(1))}</div>' if title else ""
        data_lang = f' data-lang="{html.escape(lang)}"' if lang else ""
        text = html.escape("\n".join(body))
        self.note(" ".join(body)[:300])
        return f'<div class="code"{data_lang}>{caption}<pre>{text}</pre></div>'

    def raw_html(self, chunk: str) -> str:
        """Raw HTML passes through. A <picture> that swaps on the system colour
        scheme follows the site's own theme instead."""
        picture = re.search(
            r'<source media="\(prefers-color-scheme: dark\)" srcset="([^"]+)">\s*<img src="([^"]+)" alt="([^"]*)">', chunk
        )
        if picture:
            dark, light, alt = picture.groups()
            return (
                f'<p class="figure"><img class="only-dark" src="{dark}" alt="{alt}">'
                f'<img class="only-light" src="{light}" alt="{alt}"></p>'
            )
        return chunk


def counts_from_code(renderer: "Renderer", _lines: list[str]) -> str:
    """The case counts of the code block just above, as bars."""
    numbers = {}
    for line in getattr(renderer, "last_code", []):
        m = re.match(r"^(case\w+|rv_missing|total_cds)\s+(\d+)", line.strip())
        if m:
            numbers[m.group(1)] = int(m.group(2))
    numbers.pop("case1_initial", None)
    return case_bars([("CDC1551 against H37Rv", numbers)],
                     "The same numbers, drawn: most of Bakta's calls agree with H37Rv, about a quarter overlap it with other "
                     "coordinates and are settled, and very few are new.")


def counts_from_table(renderer: "Renderer", _lines: list[str]) -> str:
    """The case counts of the table just above, one bar per genome column."""
    head, body = renderer.last_table
    keys = {"case2": "case2", "case3": "case3", "case4": "case4", "case1_kept": "case1_kept", "rv_missing": "rv_missing", "total_cds": "total_cds"}
    genomes = []
    for col in range(1, len(head)):
        numbers = {}
        for row in body:
            key = re.match(r"^(\w+)", row[0]).group(1)
            if key in keys:
                numbers[keys[key]] = int(row[col].replace(",", ""))
        genomes.append((head[col], numbers))
    return case_bars(genomes, "The same table, drawn. The number that differs by two orders of magnitude is case 1: the "
                               "prophages and the LEE island that O157:H7 carries and K-12 does not.")


def reference_steps(renderer: "Renderer", lines: list[str]) -> str:
    """PREPARE_REFERENCE's three steps, with the counts the page reports."""
    text = "\n".join(lines)
    grab = lambda pattern: int(re.search(pattern + r"\s+(\d+)", text).group(1))
    return prepare_reference({
        "proteins": grab(r"reference proteins extracted from the K-12 GenBank"),
        "bakta": grab(r"CDS Bakta called on that same reference genome"),
        "false": grab(r"Bakta false positives derived from the difference"),
    })


DERIVED = {"case-counts": counts_from_code, "case-table": counts_from_table, "prepare-reference": reference_steps}
FIGURES["results-map"] = results_map()
FIGURES["excerpt"] = EXCERPT


# ── Pages ──


def header(root: str, current: str) -> str:
    def nav(target: str, label: str, key: str, optional: bool = False) -> str:
        cls = ' class="nav-optional"' if optional else ""
        here = ' aria-current="page"' if key == current else ""
        return f'      <a{cls} href="{target}"{here}>{label}</a>\n'

    return (
        '<header class="top">\n  <div class="wrap">\n'
        f'    <a class="brand" href="{root}index.html" aria-label="AMAP home">\n'
        f'      <img src="{root}assets/amap-logo.png" alt="" width="57" height="40"><span>AMAP</span>\n'
        "    </a>\n"
        '    <nav class="nav" aria-label="Site">\n'
        + nav(f"{root}index.html#start", "Start here", "start")
        + nav(f"{root}how-it-works.html", "How it works", "how", optional=True)
        + nav(f"{root}getting-started/tutorial.html", "Tutorial", "tutorial", optional=True)
        + nav(f"{root}{INDEX}", "Docs", "docs")
        + f'      <a class="gh-link" href="{REPO}">{GH_SVG}<span>GitHub</span></a>\n'
        + f'      <button class="icon-btn" type="button" data-theme-toggle aria-label="Switch between light and dark">\n        {SUN_SVG}\n      </button>\n'
        "    </nav>\n  </div>\n</header>\n"
    )


def footer(root: str) -> str:
    return (
        '<footer class="foot">\n  <div class="wrap">\n    <div>\n      <h3>Cite</h3>\n'
        "      <blockquote>Ruiz-Rodriguez, P., &amp; Coscolla, M. (2026). AMAP (Version 1.0.1) [Computer software].</blockquote>\n"
        '      <p class="doi">Software: <a href="https://doi.org/10.5281/zenodo.8024447">doi:10.5281/zenodo.8024447</a></p>\n    </div>\n'
        "    <div>\n      <h3>Learn more</h3>\n      <ul>\n"
        f'        <li><a href="{root}getting-started/tutorial.html">Tutorial</a></li>\n'
        f'        <li><a href="{root}{INDEX}">Documentation</a></li>\n'
        f'        <li><a href="{root}how-it-works.html">How it works</a></li>\n'
        f'        <li><a href="{root}user-guide/parameters.html">Parameters</a></li>\n'
        f'        <li><a href="{REPO}/blob/main/CHANGELOG.md">Changelog</a></li>\n'
        "      </ul>\n    </div>\n    <div>\n      <h3>Project</h3>\n      <ul>\n"
        f'        <li><a href="{REPO}">GitHub</a></li>\n'
        f'        <li><a href="{REPO}/issues">Report a bug</a></li>\n'
        '        <li><a href="https://hub.docker.com/r/paururo/amap">Docker image</a></li>\n'
        "        <li>GPL-3.0 &middot; PathoGenOmics Lab</li>\n      </ul>\n    </div>\n  </div>\n</footer>\n"
    )


def sidebar(root: str, current: str | None) -> str:
    out = [
        '<aside class="docs-nav" aria-label="Documentation">\n',
        '  <div class="docs-search">\n'
        '    <input type="search" placeholder="Search the docs" aria-label="Search the documentation" data-docs-search>\n'
        '    <div class="docs-results" data-docs-results hidden></div>\n  </div>\n',
        '  <details class="docs-menu" open>\n    <summary>Pages</summary>\n',
    ]
    section = None
    for sec, path, title, _ in PAGES:
        if sec != section:
            if section is not None:
                out.append("    </ul>\n")
            out.append(f'    <p class="docs-section">{sec}</p>\n    <ul>\n')
            section = sec
        here = ' aria-current="page"' if path == current else ""
        out.append(f'      <li><a href="{root}{href(path)}"{here}>{title}</a></li>\n')
    out.append("    </ul>\n  </details>\n</aside>\n")
    return "".join(out)


def page(root: str, title: str, description: str, body: str, figures: bool = False) -> str:
    return (
        "<!doctype html>\n"
        '<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"<title>{html.escape(title)}</title>\n"
        f'<meta name="description" content="{html.escape(description)}">\n'
        f'<link rel="icon" href="{root}assets/favicon.png" type="image/png">\n'
        f'<link rel="stylesheet" href="{root}assets/site.css">\n'
        f"{THEME_BOOT}\n</head>\n<body>\n"
        '<a class="skip" href="#content">Skip to the content</a>\n'
        "<!-- Generated by tools/build_docs.py from the markdown next to this file. Edit that, then run it. -->\n"
        + body
        + f'<script src="{root}assets/search-index.js"></script>\n'
        + f'<script src="{root}assets/site.js"></script>\n'
        + (f'<script src="{root}assets/figure-data.js"></script>\n<script src="{root}assets/figures.js"></script>\n' if figures else "")
        + "</body>\n</html>\n"
    )


def source_lines(path: str) -> list[str]:
    if path == CHANGELOG[0]:
        return CHANGELOG[1].read_text(encoding="utf-8").splitlines()
    return strip_front_matter((DOCS / path).read_text(encoding="utf-8").splitlines())


def doc_page(position: int) -> tuple[Path, str, list[dict]]:
    section, path, title, description = PAGES[position]
    renderer = Renderer(path)
    content = renderer.blocks(source_lines(path))
    root = "../" * path.count("/")
    toc = "".join(
        f'<li class="toc-{level}"><a href="#{slug}">{html.escape(text)}</a></li>'
        for level, text, slug in renderer.headings
        if level in (2, 3)
    )
    prev_page = PAGES[position - 1] if position > 0 else None
    next_page = PAGES[position + 1] if position + 1 < len(PAGES) else None
    pager = '<nav class="pager" aria-label="Previous and next page">'
    if prev_page:
        pager += f'<a class="prev" href="{root}{href(prev_page[1])}"><span>Previous</span>{prev_page[2]}</a>'
    if next_page:
        pager += f'<a class="next" href="{root}{href(next_page[1])}"><span>Next</span>{next_page[2]}</a>'
    pager += "</nav>"
    edit_path = CHANGELOG[1].as_posix() if path == CHANGELOG[0] else f"docs/{path}"
    current = "tutorial" if path == "getting-started/tutorial.md" else "docs"
    body = (
        header(root, current)
        + '<div class="wrap docs-layout">\n'
        + sidebar(root, path)
        + f'<main class="prose" id="content">\n<p class="crumb"><a href="{root}{INDEX}">Docs</a> / {section}</p>\n'
        + content
        + f'\n<p class="edit"><a href="{REPO}/blob/main/{edit_path}">{GH_SVG}Edit this page on GitHub</a></p>\n'
        + pager
        + "\n</main>\n"
        + (f'<aside class="docs-toc" aria-label="On this page"><p>On this page</p><ul>{toc}</ul></aside>\n' if toc else "")
        + "</div>\n"
        + footer(root)
    )
    entries = [
        {
            "page": title,
            "section": s["heading"] if s["level"] > 1 else "",
            "url": f"{path[:-3]}.html" + (f"#{s['slug']}" if s["level"] > 1 else ""),
            "text": " ".join(s["text"])[:600],
        }
        for s in renderer.sections
    ]
    html_page = page(root, f"{title} | AMAP docs", description, body, figures=bool(renderer.figures))
    return DOCS / f"{path[:-3]}.html", html_page, entries


def index_page() -> tuple[Path, str]:
    cards = []
    section = None
    for sec, path, title, line in PAGES:
        if sec != section:
            if section is not None:
                cards.append("</div>\n")
            cards.append(f'<h2 class="docs-group">{sec}</h2>\n<div class="deep">\n')
            section = sec
        cards.append(f'<a href="{href(path)}"><b>{title}</b><span>{line}</span></a>\n')
    cards.append("</div>\n")
    body = (
        header("", "docs")
        + '<div class="wrap docs-layout docs-home">\n'
        + sidebar("", None)
        + '<main class="prose" id="content">\n<h1>Documentation</h1>\n'
        + "<p>The long versions: every option, every file, and how each decision is made. "
        + 'New to AMAP? <a href="index.html#start">Start here</a>, then follow <a href="getting-started/tutorial.html">the tutorial</a>.</p>\n'
        + "".join(cards)
        + "</main>\n</div>\n"
        + footer("")
    )
    return DOCS / INDEX, page("", "Documentation | AMAP", "AMAP's documentation, page by page.", body)


def build() -> dict[Path, str]:
    files: dict[Path, str] = {}
    entries: list[dict] = []
    for position in range(len(PAGES)):
        if not PAGES[position][1].endswith(".md"):
            continue
        path, text, page_entries = doc_page(position)
        files[path] = text
        entries += page_entries
    index_path, index_text = index_page()
    files[index_path] = index_text
    files[SEARCH_INDEX] = (
        "/* Generated by tools/build_docs.py: the documentation, section by section, for the search box. */\n"
        "window.AMAP_DOCS_INDEX = " + json.dumps(entries, ensure_ascii=False, indent=0) + ";\n"
    )
    return files


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true", help="exit 1 if a generated page differs from what is committed")
    args = parser.parse_args()
    files = build()
    stale = [p for p, text in files.items() if not p.exists() or p.read_text(encoding="utf-8") != text]
    if args.check:
        if stale:
            print("These pages no longer match their markdown; run python3 tools/build_docs.py:", file=sys.stderr)
            for p in stale:
                print(f"  {p}", file=sys.stderr)
            return 1
        print(f"docs/: {len(files)} generated files match their markdown")
        return 0
    for p in stale:
        p.write_text(files[p], encoding="utf-8")
    print(f"docs/: wrote {len(stale)} of {len(files)} generated files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
