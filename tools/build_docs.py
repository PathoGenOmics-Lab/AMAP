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
    ("User guide", "user-guide/amap-report.md", "The interactive report", "amap_report.html, question by question, with or without eggNOG."),
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
    "contig-end": case_figure("contig-end", "CDC1551 in 40 contigs: Rv1239c's frame reaches the end of contig_10 with no stop codon, and this piece of the gene, Rv1239c_2, is marked partial at its 3' end. The other piece, at the end of contig_9, is marked partial at its 5' end."),
    "contig-end-2": case_figure("contig-end-2", "The other piece of Rv1239c in the CDC1551 draft: its last 144 residues, placed by the second pass at the very end of contig_9. Bakta's short call over it is case 4, and the piece wins."),
    "case2": case_figure("case2", "Case 2 in CDC1551: dnaN, where Bakta and the aligned protein agree base for base."),
    "case3": case_figure("case3", "Case 3: fhaA, with the same ends as Rv0020c but six residues fewer."),
    "case4-start": case_figure("case4-start", "A single overlap: Bakta's dnaA starts 33 bp late, and the reference's start is kept."),
    "case4-split": case_figure("case4-split", "Multiple overlaps: two Bakta calls on Rv0007, which comes back as one gene."),
    "case4-longer": case_figure("case4-longer", "A much longer call: 3,303 bp past an alignment that stops short of its protein's end, so Bakta's call is kept, and it takes the protein's name, Rv2024c."),
    "case1-kept": case_figure("case1-kept", "Case 1 kept: four calls in a stretch CDC1551 has and H37Rv lacks."),
    "case1-dropped": case_figure("case1-dropped", "Case 1 removed: a call that Bakta's known false positives cover end to end."),
    "missing": case_figure("missing", "A missing reference gene: Rv0031, put back between two of Bakta's calls."),
    "three-inputs": (
        '<figure class="fig case" data-inputs><div class="case-fig" data-fig-stage></div>'
        "<figcaption>The three inputs over the first 12.4 kb of CDC1551: Bakta's calls, the H37Rv proteins miniprot placed, and "
        "two of the calls Bakta makes on H37Rv itself that its curated annotation does not have, aligned here. Point at a gene "
        "for its coordinates.</figcaption></figure>"
    ),
    "outcomes": (
        '<figure class="fig outcomes" data-outcomes><div class="legend"><span><i class="chip bakta"></i>Bakta\'s call</span>'
        '<span><i class="chip ref"></i>reference protein, aligned</span><span><i class="chip amap"></i>what AMAP keeps</span></div>'
        '<div class="oc-grid">'
        + "".join(
            f'<a class="oc" href="{link}" data-mini="{key}"><span class="oc-head"><span class="case-tag {tag}">{case}</span><b>{title}</b></span>'
            f'<span class="oc-fig" data-fig-stage></span><span class="oc-text">{text}</span><span class="oc-more">{more}</span></a>'
            for key, tag, case, title, text, link, more in [
                ("case2", "c2", "Case 2", "Perfect match", "Same start, end, strand and length as an aligned protein. Kept, tied to it.", "../how-it-works.html#agree", "Step through it"),
                ("case3", "c3", "Case 3", "Imperfect match", "Same start, end and strand; the alignment's length differs. Kept, tied to the protein.", "correction-algorithm.html#1-classification", "The rule"),
                ("case4-start", "c4", "Case 4", "Conflicting overlap", "Overlaps an aligned protein with other coordinates. The reference's win, unless the call is much longer.", "../how-it-works.html#late-start", "Step through it"),
                ("case1-kept", "c1", "Case 1", "Potentially novel", "Nothing aligns. Kept, unless Bakta is known to call it wrongly or the reference genome has its protein.", "../how-it-works.html#new", "Step through it"),
                ("missing", "cm", "Missing", "Put back", "A reference protein aligns where Bakta called nothing. It goes in at the alignment.", "../how-it-works.html#missed", "Step through it"),
            ]
        )
        + '<div class="oc oc-harder"><span class="oc-head"><b>And the harder ones</b></span><ul>'
        + "".join(
            f'<li><a href="../how-it-works.html#{key}">{label}</a></li>'
            for key, label in [("paralog", "A protein on its paralog's locus"), ("small", "Proteins too small for miniprot"),
                               ("contig-end", "Genes a contig end cuts"), ("frameshift", "Frameshifts and lost stops"),
                               ("copies", "Two real copies of a gene")]
        )
        + "</ul></div>"
        + "</div><figcaption>The five outcomes, each on real genes of CDC1551: the same windows as the case pages, "
        "drawn small. Each card opens its story, step by step.</figcaption></figure>"
    ),
    "dup-family": (
        '<figure class="fig dup-fig" data-dup="family"><div class="case-fig" data-fig-stage></div>'
        '<div class="legend"><span><i class="chip node-in"></i>in the family</span>'
        '<span><i class="chip link-kept"></i>a link: identical residues cover 80% or more of one of the two</span>'
        '<span><i class="chip link-near"></i>a hit under 80%: no link</span>'
        '<span><i class="chip node-out"></i>a protein outside the family</span></div>'
        "<figcaption>CDC1551's corrected proteome against itself, run as GET_DUPLICATIONS runs it: six transposases make one "
        "family. Each line is a pair Diamond matched, with the larger of the two shares its identical residues cover; the "
        "size of a protein is its length. Point at a protein or a line for its numbers.</figcaption></figure>"
    ),
    "dup-communities": (
        '<figure class="fig dup-fig dup-pair" data-dup="communities"><div class="dup-sides">'
        '<div class="dup-side" data-net><p class="dup-head" data-net-head></p><div class="case-fig" data-fig-stage></div></div>'
        '<div class="dup-side" data-net><p class="dup-head" data-net-head></p><div class="case-fig" data-fig-stage></div></div>'
        '<p class="dup-bridge" data-net-bridge></p></div>'
        '<div class="table-wrap"><table class="dup-metrics" data-net-metrics></table></div>'
        '<div class="legend"><span><i class="chip link-kept"></i>a link in this genome</span>'
        '<span><i class="chip link-ghost"></i>a link only the other genome has</span></div>'
        "<figcaption>The same family in GET_ALLDUP, pooled with M. bovis AF2122/97's families: one community in each genome, "
        "joined by their hits, and the row of <code>homologous_communities.tsv</code> they make. Each protein of M. bovis sits "
        "where its closest match in CDC1551 does.</figcaption></figure>"
    ),
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
        ("pan", "PANGENE/", "PANGENE", "with <code>--skip_functional_annotation</code>"),
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
    "contig_9   AMAP  CDS  4160    4591    .  -  0  ID=L4_CDC1551d_Rv1239c;&hellip;;Dbxref=AMAP_ids:Rv1239c;"
    '<mark class="part">partial=true;end_range=4591,.</mark>\n'
    "contig_10  AMAP  CDS  101478  102146  .  +  0  ID=L4_CDC1551d_Rv1239c_2;&hellip;;Dbxref=AMAP_ids:Rv1239c;"
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
    '                     /locus_tag="L4_CDC1551d_Rv1239c_2"</pre>'
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


# ── The sample sheet and the checks on it ──
#
# Drawn from bin/validate_inputs.py itself: its messages come from running it
# on the sheets shown, and the ranges from its own comparisons, so a change
# there leaves these pages stale until they are rebuilt.

VALIDATOR = Path("bin/validate_inputs.py")
SHEET = [("sampleId", "fasta"), ("G26689", "/data/genomes/G26689.fasta"), ("G26690", "/data/genomes/G26690.fasta"),
         ("L4_G12345", "./assemblies/L4_G12345.fa")]


def run_validator(sheet: str) -> list[str]:
    """What validate_inputs.py prints, as a run of --validate_only shows it,
    for a sample sheet: its error lines and its verdict."""
    import contextlib
    import io
    import importlib.util
    import tempfile

    spec = importlib.util.spec_from_file_location("validate_inputs", VALIDATOR)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "samples.tsv"
        path.write_text(sheet)
        out = io.StringIO()
        argv = sys.argv
        sys.argv = ["validate_inputs.py", "--tsv", str(path)]
        try:
            with contextlib.redirect_stdout(out):
                module.main()
        except SystemExit:
            pass
        finally:
            sys.argv = argv
    lines = out.getvalue().replace(str(path), "samples.tsv").splitlines()
    return [line for line in lines if line.strip().startswith(("✗", "FAILED"))]


def validation_ranges() -> dict[str, float]:
    """The limits validate_inputs.py checks a genome against, read off its code."""
    code = VALIDATOR.read_text(encoding="utf-8")
    patterns = {
        "size_fail_low": r"(?<!el)if total_size < ([\d_]+):",
        "size_warn_low": r"elif total_size < ([\d_]+):",
        "size_fail_high": r"(?<!el)if total_size > ([\d_]+):",
        "size_warn_high": r"elif total_size > ([\d_]+):",
        "n_fail": r"(?<!el)if n_pct > ([\d.]+):",
        "n_warn": r"elif n_pct > ([\d.]+):",
        "gc_fail_low": r"(?<!el)if gc_pct < ([\d.]+) or gc_pct > [\d.]+:",
        "gc_fail_high": r"(?<!el)if gc_pct < [\d.]+ or gc_pct > ([\d.]+):",
        "gc_warn_low": r"elif gc_pct < ([\d.]+) or gc_pct > [\d.]+:",
        "gc_warn_high": r"elif gc_pct < [\d.]+ or gc_pct > ([\d.]+):",
        "contigs_warn_high": r"(?<!el)if n_contigs > ([\d_]+):",
        "contigs_warn": r"elif n_contigs > ([\d_]+):",
    }
    out = {}
    for key, pattern in patterns.items():
        m = re.search(pattern, code)
        if not m:
            raise SystemExit(f"tools/build_docs.py: {VALIDATOR} no longer checks {key} the way validation_ranges() reads it")
        out[key] = float(m.group(1).replace("_", ""))
    return out


# ── The container ──
#
# What the image holds, read off docker/Dockerfile, so the drawing cannot
# name a version the image does not have.

DOCKERFILE = Path("docker/Dockerfile")


def image_contents() -> dict:
    text = DOCKERFILE.read_text(encoding="utf-8")
    final = text[text.rindex("\nFROM "):]

    def version(pattern: str, where: str = text) -> str:
        m = re.search(pattern, where)
        if not m:
            raise SystemExit(f"tools/build_docs.py: {DOCKERFILE} no longer says {pattern!r}; update image_contents()")
        return m.group(1)

    tools = [
        ("Bakta", version(r"\bbakta=([\d.]+)", final)),
        ("miniprot", version(r"--branch v([\d.]+) [^\n]*lh3/miniprot")),
        ("Diamond", version(r"\bdiamond=([\d.]+)", final)),
        ("eggNOG-mapper", version(r"\beggnog-mapper=([\d.]+)", final)),
        ("pangene", version(r"pangene/releases/download/v([\d.]+)")),
        ("pathotypr", version(r"\bpathotypr=([\d.]+)")),
        ("gfatools", ""),
    ]
    runtimes = [("Python", version(r"\bpython=([\d.]+)", final)), ("R", version(r"\br-base=([\d.]+)", final))]
    reference = sorted(set(re.findall(r"-O /reference/([\w.]+)", text)) | set(re.findall(r"-P /reference \S+/([\w.]+)", text)))
    return {"tools": tools, "runtimes": runtimes, "reference": reference, "bakta_db": version(r"bakta_db download [^\n]*--type (\w+)")}


def image_figure(_renderer: "Renderer", _lines: list[str]) -> str:
    image = image_contents()
    nf = re.search(r"nextflowVersion\s*=\s*'>=([\d.]+?)(?:\.0)?'", Path("nextflow.config").read_text(encoding="utf-8")).group(1)
    pill = lambda name, v: f'<span class="img-pill"><b>{html.escape(name)}</b>{(" " + html.escape(v)) if v else ""}</span>'
    host = (
        '<div class="img-host"><p class="img-title">What you install</p><ul>'
        f"<li><b>Nextflow</b> {nf} or newer, on Java 17 or newer</li>"
        "<li><b>Docker</b>, or <b>Singularity</b> or <b>Apptainer</b> on a cluster</li></ul>"
        '<p class="img-title img-later">What the first run builds</p><ul>'
        "<li><b>The eggNOG database</b>, in the folder <code>--eggnog_db</code> names, kept for every run after</li></ul></div>"
    )
    box = (
        '<div class="img-box"><p class="img-title"><code>paururo/amap</code><span>linux/amd64</span></p>'
        '<div class="img-layer"><p>Tools</p><div>' + "".join(pill(n, v) for n, v in image["tools"]) + "</div></div>"
        '<div class="img-layer"><p>Runtimes</p><div>' + "".join(pill(n, v) for n, v in image["runtimes"]) + "</div></div>"
        '<div class="img-layer"><p>Data</p><div>'
        + pill("Bakta database", image["bakta_db"])
        + f'<span class="img-pill img-note"><b><code>/reference/</code></b> {len(image["reference"])} files: the H37Rv genome, proteome and GenBank, '
        "Bakta's false positives on it, Mycobrowser's table and the lineage markers</span></div></div></div>"
    )
    return (
        '<figure class="fig image-fig"><div class="img-grid">'
        + host
        + '<div class="img-arrow" aria-hidden="true"><span>pulled once</span><i class="img-left">&#8592;</i><i class="img-up">&#8593;</i></div>'
        + box
        + "</div><figcaption>Two things to install; everything AMAP runs comes in one image, its versions read off "
        "<code>docker/Dockerfile</code>. Only the eggNOG database stays outside: too large to ship, it is built on the first "
        "run that needs it.</figcaption></figure>"
    )


# ── A run on a cluster ──
#
# The scheduler's settings come from nextflow.config's profiles, and the
# launcher's from the job script the HPC page gives, so neither can drift.


def profile_settings(name: str) -> dict[str, str]:
    config = Path("nextflow.config").read_text(encoding="utf-8")
    block = re.search(rf"^    {name}\s*\{{(.*?)^    \}}", config[config.index("profiles {"):], re.M | re.S)
    if not block:
        raise SystemExit(f"tools/build_docs.py: nextflow.config has no profile {name}")
    return dict(re.findall(r"^\s+([\w.]+)\s*=\s*'?([^'\n]+?)'?\s*(?://.*)?$", block.group(1), re.M))


def local_processes() -> list[str]:
    config = Path("nextflow.config").read_text(encoding="utf-8")
    return re.findall(r"withName: '(\w+)' \{\s*executor = 'local'", config)


def hpc_figure(_renderer: "Renderer", lines: list[str]) -> str:
    slurm = profile_settings("slurm")
    others = {name: profile_settings(name) for name in ("sge", "pbs", "lsf")}
    script = "\n".join(lines)
    sbatch = dict(re.findall(r"^#SBATCH --([\w-]+)=(\S+)", script, re.M))
    rate = lambda s: s["executor.submitRateLimit"].replace("/1min", " a minute")
    other_rates = sorted({rate(s) for s in others.values()})
    memory = slurm["params.max_memory"].replace(".GB", " GB")
    hours = int(sbatch["time"].split(":")[0])
    launcher = (
        '<div class="hpc-box hpc-launch"><p class="hpc-title">The launcher</p>'
        f"<p><code>nextflow run</code>, in <code>tmux</code> or as a small job of its own: {sbatch['cpus-per-task']} CPUs, "
        f"{sbatch['mem'].replace('G', ' GB')}, {hours} h. It submits, waits and resubmits.</p>"
        '<p class="hpc-small">Runs here too: ' + ", ".join(f"<code>{p}</code>" for p in local_processes()) + "</p></div>"
    )
    queue = (
        '<div class="hpc-box hpc-queue"><p class="hpc-title">The scheduler</p>'
        f'<p>Up to <b>{slurm["executor.queueSize"]}</b> jobs at once, submitted at up to {rate(slurm)} with <code>slurm</code> '
        f'({" or ".join(other_rates)} with <code>sge</code>, <code>pbs</code> and <code>lsf</code>).</p></div>'
    )
    nodes = (
        '<div class="hpc-box hpc-nodes"><p class="hpc-title">The compute nodes</p>'
        f"<p>One job per task, at most {slurm['params.max_cpus']} CPUs and {memory} each, in the AMAP image, "
        + ("working on the node's own disk and copying back only the results.</p>" if slurm.get("params.scratch") == "true" else "</p>")
        + '<div class="hpc-tasks"><span>per sample: Bakta, miniprot, the correction, duplications, pangenome</span>'
        "<span>once for all: eggNOG, the genomes compared, the report</span></div></div>"
    )
    shared = (
        '<div class="hpc-box hpc-shared"><p class="hpc-title">Shared storage</p><div class="hpc-places">'
        '<span><code class="hpc-key">-work-dir</code>every task\'s files, which <code>-resume</code> needs</span>'
        '<span><code class="hpc-key">NXF_SINGULARITY_CACHEDIR</code>the image, converted once</span>'
        '<span><code class="hpc-key">--eggnog_db</code>the eggNOG database, built once</span>'
        '<span><code class="hpc-key">--outdir</code>the results</span></div></div>'
    )
    return (
        '<figure class="fig hpc-fig"><div class="hpc-flow">'
        + launcher
        + '<i class="hpc-arrow" aria-hidden="true"></i>'
        + queue
        + '<i class="hpc-arrow" aria-hidden="true"></i>'
        + nodes
        + "</div>"
        + shared
        + "<figcaption>A run on a cluster, with the settings of the <code>slurm</code> profile in <code>nextflow.config</code> "
        "and the launcher of the job script below.</figcaption></figure>"
    )


# ── The report, section by section ──
#
# Pictures of the example report, taken by tools/report_shots.mjs; each card
# opens the report at its section, and check_site.py checks those anchors.

REPORT_SECTIONS = [
    ("overview", "Overview", "Dataset composition, then sample quality and annotation quality: where an outlier sample shows itself."),
    ("amap-annotation-results", "Annotation results", "Genes per lineage, the genes the reference does not have, the correction summary and the reference genes Bakta missed."),
    ("functional-annotation", "Functional annotation", "COG and Mycobrowser categories, the comparison with the reference, PE/PPE families, KEGG pathways and pseudogenes."),
    ("gene-duplications", "Duplications", "Community sizes, the communities shared across genomes and their network metrics."),
    ("pangenome-analysis", "Pangenome", "Core and accessory genes, the gene frequency spectrum and the presence and absence of every gene."),
    ("amap-cds-stats", "CDS statistics", "How many genes fall in each case, lineage by lineage: box and violin plots with Dunn's tests, the share of each case, a heatmap and a table per case."),
    ("master-table-explorer", "Master table", "Every gene of every sample in one table to search, sort and export."),
]


def report_map(_renderer: "Renderer", _lines: list[str]) -> str:
    cards = "".join(
        f'<a class="rp-card" href="../examples/demo_report.html#{key}"><img src="../assets/report/{key}.webp" alt="The {title.lower()} section of the example report" '
        f'width="460" height="320" loading="lazy"><span class="rp-text"><b>{title}</b><span>{text}</span></span></a>'
        for key, title, text in REPORT_SECTIONS
    )
    return (
        '<figure class="fig report-map"><p class="rp-head"><span>final_report.html, section by section</span>'
        '<span class="tag-synthetic">synthetic data</span></p>'
        f'<div class="rp-grid">{cards}</div>'
        "<figcaption>The example report's seven sections, in their order. Each picture opens the report there; its twelve "
        "genomes are made up, so read the layout, not the numbers.</figcaption></figure>"
    )


# ── The command builder ──
#
# The front page's, taken from docs/index.html so the two are always the same.


def builder_figure(_renderer: "Renderer", _lines: list[str]) -> str:
    front = (DOCS / "index.html").read_text(encoding="utf-8")
    m = re.search(r'\n( *)<div class="builder" data-builder>(.*?)\n\1</div>\n', front, re.S)
    if not m:
        raise SystemExit("tools/build_docs.py: docs/index.html has no command builder to copy")
    inner = "\n".join(line[len(m.group(1)):] if line.startswith(m.group(1)) else line for line in m.group(2).splitlines())
    return (
        '<figure class="fig builder-fig"><div class="builder" data-builder>' + inner + "\n</div>"
        "<figcaption>The command for your machine and your genomes, with what each option does. "
        'The same as on the <a href="../index.html#start">front page</a>.</figcaption></figure>'
    )


# ── When a task fails ──


def failure_figure(_renderer: "Renderer", _lines: list[str]) -> str:
    """What nextflow.config's error strategy does with a failed task."""
    config = Path("nextflow.config").read_text(encoding="utf-8")
    rule = re.search(
        r"errorStrategy = \{ task\.exitStatus in \[([\d, ]+)\] \? 'retry' : \(task\.attempt <= (\d+) \? 'retry' : 'ignore'\) \}\s*maxRetries = (\d+)",
        config,
    )
    # Each withName block, up to the next one or the end of process { }.
    blocks = dict(re.findall(r"withName: '(\w+)' \{(.*?)\n    \}", config, re.S))
    bakta = re.search(r"memory\s*=\s*\{ (\d+)\.GB \* task\.attempt \}", blocks.get("BAKTA_ANNOTATOR", ""))
    terminate = [name for name, body in blocks.items() if "errorStrategy = 'terminate'" in body]
    eggnog = re.search(r"errorStrategy = \{ task\.exitStatus in \[[\d, ]+\] \? 'retry' : 'terminate' \}", blocks.get("ANNOT_EGGNOGMAP_BATCH", ""))
    if not (rule and bakta and eggnog):
        raise SystemExit("tools/build_docs.py: nextflow.config no longer handles failures the way failure_figure() reads it")
    codes = [c.strip() for c in rule.group(1).split(",")]
    retries = int(rule.group(3))
    base = int(bakta.group(1))
    ladder = " GB, then ".join(str(base * k) for k in range(1, retries + 2)) + " GB"
    code_list = ", ".join(codes[:-1]) + " or " + codes[-1]
    lane = lambda cls, title, steps: (
        f'<div class="fail-lane {cls}"><p class="fail-when">{title}</p><ol>'
        + "".join(f"<li>{s}</li>" for s in steps)
        + "</ol></div>"
    )
    return (
        '<figure class="fig fail-fig"><div class="fail-lanes">'
        + lane("fail-mem", f"Killed by the machine or the scheduler: exit {code_list}", [
            f"Tried again, up to {retries} more times, asking for more memory each time: Bakta {ladder}.",
            "Never past <code>--max_memory</code>: a request above it is lowered to it.",
            "Still failing after that: the run stops. Raise <code>--max_memory</code>, or the memory Docker Desktop has.",
        ])
        + lane("fail-other", "Any other exit, in a step run per sample", [
            f"Tried again, {int(rule.group(2))} more times.",
            "Still failing: that sample's task is set aside, and the run carries on with the others.",
            "Its work directory, named in the log, says why.",
        ])
        + "</div>"
        + '<p class="fail-except">The steps run once for all the samples are retried too, and then stop the run. Some stop it '
        + "at once: " + ", ".join(f"<code>{p}</code>" for p in terminate)
        + " on any failure, and the eggNOG annotation on anything but memory or time.</p>"
        + "<figcaption>What <code>nextflow.config</code> does with a task that fails.</figcaption></figure>"
    )


# The genomes of tests/benchmark, measured with validate_inputs.py's own
# parser: every one passes with no warning.
BENCHMARK_GENOMES = [
    # (label, size in bp, sequences, GC %, N %)
    ("M. tuberculosis H37Rv", 4411532, 1, 65.61, 0.0),
    ("M. tuberculosis CDC1551", 4403837, 1, 65.61, 0.001),
    ("CDC1551 cut into 40 contigs", 4403837, 40, 65.61, 0.001),
    ("M. bovis AF2122/97", 4349904, 1, 65.63, 0.0),
    ("E. coli K-12 MG1655", 4641652, 1, 50.79, 0.0),
    ("E. coli O157:H7 Sakai", 5594605, 3, 50.48, 0.0),
]


def sheet_figure() -> str:
    head, *rows = SHEET
    tab = '<td class="sheet-tab" aria-label="a tab">&#8677;</td>'
    body = "".join(f'<tr><td>{html.escape(a)}</td>{tab}<td class="sheet-path">{html.escape(b)}</td></tr>' for a, b in rows)
    table = (
        '<table class="sheet-table"><thead><tr><th><span class="sheet-n">1</span>sampleId</th>'
        f'{tab}<th><span class="sheet-n">3</span>fasta</th></tr></thead><tbody>{body}</tbody></table>'
    )
    notes = (
        '<ol class="sheet-notes">'
        "<li><b><code>sampleId</code></b>: unique, and carried into every folder and file name. Letters, digits and underscores.</li>"
        "<li><b>A real tab</b> between the two columns, on every line, header included. Spaces will not do.</li>"
        "<li><b><code>fasta</code></b>: the assembled genome, as an absolute path or one relative to where you launch Nextflow.</li></ol>"
    )
    good = "\n".join(f"{a}\t{b}" for a, b in SHEET) + "\n"
    spaced = "\n".join(f"{a:<12}{b}" for a, b in SHEET) + "\n"
    checks = [
        ("With spaces where the tabs go", run_validator(spaced)),
        ("With tabs, and paths that do not exist yet", run_validator(good)),
    ]
    blocks = "".join(
        f'<div class="sheet-check"><p class="sheet-when">{when}, <code>--validate_only</code> stops with:</p>'
        f'<pre>{html.escape(chr(10).join(line.strip() for line in lines))}</pre></div>'
        for when, lines in checks
    )
    return (
        '<figure class="fig sheet-fig"><div class="sheet"><div class="sheet-file"><p class="sheet-name">samples.tsv'
        '<span class="sheet-n sheet-n-tab">2</span></p>'
        f"{table}</div>{notes}</div>"
        f'<div class="sheet-checks">{blocks}</div>'
        "<figcaption>The sample sheet, with its tabs drawn. Below it, what the check prints for the two mistakes it "
        "catches most: both messages come from running <code>bin/validate_inputs.py</code> on this sheet.</figcaption></figure>"
    )


def ranges_figure() -> str:
    import math

    th = validation_ranges()

    def gauge(name: str, unit: str, lo: float, hi: float, log: bool, zones: list, ticks: list, marks: list, fmt) -> str:
        to = (lambda v: math.log10(v)) if log else (lambda v: v)
        pos = lambda v: 100 * (to(min(max(v, lo), hi)) - to(lo)) / (to(hi) - to(lo))
        out = [f'<div class="gauge"><p class="gauge-name"><b>{name}</b><span>{unit}</span></p><div class="gauge-track">']
        for a, b, kind in zones:
            out.append(f'<span class="gauge-zone gz-{kind}" style="left:{pos(a):.2f}%;width:{pos(b) - pos(a):.2f}%"></span>')
        for k, (value, label) in enumerate(marks):
            side = "up" if k % 2 == 0 else "down"
            # A label near either end of the scale reads inward from its mark.
            edge = " gm-start" if pos(value) < 12 else " gm-end" if pos(value) > 88 else ""
            out.append(
                f'<span class="gauge-mark gm-{side}{edge}" style="left:{pos(value):.2f}%"><i></i><em>{html.escape(label)}</em></span>'
            )
        out.append('</div><div class="gauge-ticks">')
        # A tick too close to the one before it is left out on a phone.
        shown = None
        for value in ticks:
            minor = shown is not None and pos(value) - shown < 12
            cls = ' class="minor"' if minor else ""
            out.append(f'<span{cls} style="left:{pos(value):.2f}%">{fmt(value)}</span>')
            if not minor:
                shown = pos(value)
        out.append("</div>")
        # On a phone the labels go under the scale, as a list.
        out.append('<ul class="gauge-list">' + "".join(f"<li>{html.escape(label)}</li>" for _value, label in marks) + "</ul></div>")
        return "".join(out)

    def size(v: float) -> str:
        return f"{v / 1e6:g} Mb" if v >= 1e6 else f"{v / 1e3:g} kb"

    mtb = [g for g in BENCHMARK_GENOMES if g[0].startswith(("M.", "CDC"))]
    coli = [g for g in BENCHMARK_GENOMES if g[0].startswith("E.")]
    def span(gs: list, i: int, f) -> str:
        low, high = f(min(g[i] for g in gs)), f(max(g[i] for g in gs))
        return low if low == high else f"{low} to {high}"

    mb = lambda v: f"{v / 1e6:.2f}"
    gauges = [
        gauge("Genome size", "log scale", 5e4, 3e7, True,
              [(5e4, th["size_fail_low"], "fail"), (th["size_fail_low"], th["size_warn_low"], "warn"),
               (th["size_warn_low"], th["size_warn_high"], "ok"), (th["size_warn_high"], th["size_fail_high"], "warn"),
               (th["size_fail_high"], 3e7, "fail")],
              [th["size_fail_low"], th["size_warn_low"], th["size_warn_high"], th["size_fail_high"]],
              [(min(g[1] for g in mtb), f"M. tuberculosis, M. bovis: {span(mtb, 1, mb)} Mb"),
               (min(g[1] for g in coli), f"E. coli: {span(coli, 1, mb)} Mb")], size),
        gauge("GC content", "% of A, C, G and T", 0, 100, False,
              [(0, th["gc_fail_low"], "fail"), (th["gc_fail_low"], th["gc_warn_low"], "warn"),
               (th["gc_warn_low"], th["gc_warn_high"], "ok"), (th["gc_warn_high"], th["gc_fail_high"], "warn"),
               (th["gc_fail_high"], 100, "fail")],
              [0, th["gc_fail_low"], th["gc_warn_low"], 50, th["gc_warn_high"], th["gc_fail_high"], 100],
              [(sum(g[3] for g in mtb) / len(mtb), f"tuberculosis complex: {span(mtb, 3, lambda v: f'{v:.1f}')}%"),
               (sum(g[3] for g in coli) / len(coli), f"E. coli: {span(coli, 3, lambda v: f'{v:.1f}')}%")],
              lambda v: f"{v:g}%"),
        gauge("N content", "% of the bases", 0, 15, False,
              [(0, th["n_warn"], "ok"), (th["n_warn"], th["n_fail"], "warn"), (th["n_fail"], 15, "fail")],
              [0, th["n_warn"], 5, th["n_fail"], 15],
              [(0, "all six: under 0.01%")], lambda v: f"{v:g}%"),
        gauge("Sequences", "contigs, log scale", 1, 2000, True,
              [(1, th["contigs_warn"], "ok"), (th["contigs_warn"], th["contigs_warn_high"], "warn"),
               (th["contigs_warn_high"], 2000, "warn2")],
              [1, 10, th["contigs_warn"], th["contigs_warn_high"], 2000],
              [(1, "closed genomes: 1"), (3, "Sakai: 3"), (40, "CDC1551 in 40 contigs")], lambda v: f"{v:,.0f}"),
    ]
    key = (
        '<div class="legend gauge-key"><span><i class="chip gz-ok"></i>passes</span><span><i class="chip gz-warn"></i>warns</span>'
        '<span><i class="chip gz-warn2"></i>warns: highly fragmented</span><span><i class="chip gz-fail"></i>fails the run</span></div>'
    )
    return (
        f'<figure class="fig gauges">{"".join(gauges)}{key}'
        "<figcaption>Where a genome has to sit, with the limits read off <code>bin/validate_inputs.py</code>, and where the six "
        "genomes of the accuracy benchmark sit: inside every range, with no warning.</figcaption></figure>"
    )


DERIVED = {
    "case-counts": counts_from_code,
    "case-table": counts_from_table,
    "prepare-reference": reference_steps,
    "sample-sheet": lambda _renderer, _lines: sheet_figure(),
    "validation-ranges": lambda _renderer, _lines: ranges_figure(),
    "image": image_figure,
    "hpc": hpc_figure,
    "failure": failure_figure,
    "report-map": report_map,
    "builder": builder_figure,
}
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


# Data a figure needs beyond figure-data.js, loaded only on the pages that draw it.
FIGURE_DATA = {"dup-family": "assets/dup-data.js", "dup-communities": "assets/dup-data.js"}


SITE_URL = "https://pathogenomics-lab.github.io/AMAP/"


def social(root: str, title: str, description: str, path: str) -> str:
    """What a link to the page shows where it is shared: its title, its line
    and the site's card, at the address Pages serves it from."""
    t, d = html.escape(title), html.escape(description)
    return (
        '<meta property="og:type" content="website">\n'
        '<meta property="og:site_name" content="AMAP">\n'
        f'<meta property="og:title" content="{t}">\n'
        f'<meta property="og:description" content="{d}">\n'
        f'<meta property="og:url" content="{SITE_URL}{path}">\n'
        f'<meta property="og:image" content="{SITE_URL}assets/social-card.png">\n'
        '<meta property="og:image:width" content="1200">\n<meta property="og:image:height" content="630">\n'
        '<meta property="og:image:alt" content="Bakta calls the genes. AMAP checks every one: real genes of M. tuberculosis CDC1551, corrected against H37Rv.">\n'
        '<meta name="twitter:card" content="summary_large_image">\n'
        '<meta name="theme-color" content="#0b0f0e">\n'
        f'<link rel="apple-touch-icon" href="{root}assets/touch-icon.png">\n'
    )


def page(root: str, title: str, description: str, body: str, figures: list[str] | None = None, path: str = "") -> str:
    data = list(dict.fromkeys(FIGURE_DATA[f] for f in figures or [] if f in FIGURE_DATA))
    return (
        "<!doctype html>\n"
        '<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"<title>{html.escape(title)}</title>\n"
        f'<meta name="description" content="{html.escape(description)}">\n'
        f'<link rel="icon" href="{root}assets/favicon.png" type="image/png">\n'
        + social(root, title, description, path)
        + f'<link rel="stylesheet" href="{root}assets/site.css">\n'
        f"{THEME_BOOT}\n</head>\n<body>\n"
        '<a class="skip" href="#content">Skip to the content</a>\n'
        "<!-- Generated by tools/build_docs.py from the markdown next to this file. Edit that, then run it. -->\n"
        + body
        + f'<script src="{root}assets/search-index.js"></script>\n'
        + f'<script src="{root}assets/site.js"></script>\n'
        + (f'<script src="{root}assets/figure-data.js"></script>\n' if figures else "")
        + "".join(f'<script src="{root}{d}"></script>\n' for d in data)
        + (f'<script src="{root}assets/figures.js"></script>\n' if figures else "")
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
    html_page = page(root, f"{title} | AMAP docs", description, body, figures=renderer.figures, path=href(path))
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
    return DOCS / INDEX, page("", "Documentation | AMAP", "AMAP's documentation, page by page.", body, path=INDEX)


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
