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
import functools
import html
import json
import re
import sys
import textwrap
from pathlib import Path

DOCS = Path("docs")
REPO = "https://github.com/PathoGenOmics-Lab/AMAP"

# The Docs tab, in reading order: (section, path, title, one line). A path to
# a .md file is built into a page; an .html path is a page written by hand and
# only gets its place in the list.
PAGES = [
    ("Getting started", "getting-started/installation.md", "Installation", "Nextflow, a container engine, and running without one."),
    ("Getting started", "getting-started/quickstart.md", "Quickstart", "The front page's four steps as a page: the test, a sample sheet, the check and the run."),
    ("Getting started", "getting-started/tutorial.md", "Tutorial", "One real genome, start to finish, with every number."),
    ("Getting started", "getting-started/inputs.md", "Preparing inputs", "The sample sheet, the genomes, and the checks run on them."),
    ("User guide", "user-guide/parameters.md", "Parameters", "Every option, its default, and when to change it."),
    ("User guide", "user-guide/running.md", "Running AMAP", "Profiles, schedulers and larger runs."),
    ("User guide", "user-guide/outputs.md", "Outputs", "Every folder and file, and what a healthy run looks like."),
    ("User guide", "user-guide/report.md", "The full report", "What each section of final_report.html shows."),
    ("User guide", "user-guide/amap-report.md", "The report by question", "amap_report.html, question by question, with or without eggNOG."),
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
    '<script>document.documentElement.classList.add("js"); try { document.documentElement.dataset.theme = '
    'localStorage.getItem("amap.theme") === "light" ? "light" : "dark"; } catch (e) { '
    'document.documentElement.dataset.theme = "dark"; }</script>'
)
INLINE_TAGS = re.compile(r"</?(kbd|br|sup|sub|strong|em|b|i)\s*/?>", re.I)
# What figures.js looks for to draw: a page with none of these needs neither it nor its data.
FIGURE_HOOK = re.compile(r"\bdata-(fig|inputs|outcomes|metro|dup)(?=[\s=>])")
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
                ("case3", "c3", "Case 3", "Imperfect match", "Same start, end and strand; the alignment's length differs. Kept, tied to the protein.", "../how-it-works.html#same-ends", "Step through it"),
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
        "<figcaption>The same family in GET_ALLDUP, pooled with <i>M.&nbsp;bovis</i> AF2122/97's families: one community in each genome, "
        "joined by their hits, and the row of <code>homologous_communities.tsv</code> they make. Each protein of <i>M.&nbsp;bovis</i> sits "
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
        ("3", "Keep what Bakta invents", "Every Bakta call that overlaps no curated gene of its own sequence: the calls it is prone to make on your samples too.",
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


# Anchors that links written before a heading was renamed still use: each
# page's old id, and the heading it now belongs to. The changelog keeps
# GitHub's anchor for every heading of CHANGELOG.md on its own (see
# Renderer.heading), since links to those are written on GitHub.
LEGACY_IDS = {
    "getting-started/quickstart.md": {
        "1-write-a-sample-sheet": "2-write-a-sample-sheet",
        "2-check-the-inputs": "3-check-the-inputs",
        "3-run-it": "4-run-it",
        "4-read-the-results": "5-read-the-results",
        "5-resume-an-interrupted-run": "6-resume-an-interrupted-run",
    },
    "how-it-works/workflow.md": {"10-report": "10-the-full-report", "11-interactive-report": "11-the-report-by-question"},
    "about/citation.md": {"citing-amap": "citation"},
    "deployment/hpc.md": {"running-on-a-cluster": "hpc-clusters"},
    "user-guide/custom-reference.md": {"using-your-own-reference": "custom-reference"},
    "user-guide/report.md": {"the-report": "the-full-report-final_reporthtml"},
    "user-guide/amap-report.md": {"the-interactive-report": "the-report-by-question-amap_reporthtml"},
    "how-it-works/index.md": {"how-it-works": "overview"},
}


def anchor_alias(old: str, new: str) -> str:
    """An old anchor that leads to a heading's current one; site.js puts the
    current one in the address."""
    return f'<span class="anchor-alias" id="{html.escape(old)}" data-to="{html.escape(new)}"></span>'


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


# ── Inline code ──
#
# A chip of code never splits inside a word. Up to SHORT_CODE characters it
# does not wrap at all; a longer one wraps only at a space or after a /, _ or
# ., so a path or a long option name comes apart where it reads in pieces.
# Each piece is a span that does not wrap, with a <wbr> between pieces: the
# browser's own rules would also break after a hyphen, '--|skip_typing'. A
# piece longer than LONG_TOKEN with nothing to break at, a commit hash, may
# break anywhere, or it would run off a phone's screen.

SHORT_CODE = 24


def code_pieces(text: str) -> str:
    out = []
    for token in re.split(r"( +)", text):
        if not token:
            continue
        if token.isspace():
            out.append(token)
            continue
        pieces = re.findall(r"[^/_.]*[/_.]+|[^/_.]+", token)
        out.append("<wbr>".join(f'<span class="{"nw" if len(p) <= LONG_TOKEN else "brk"}">{html.escape(p)}</span>' for p in pieces))
    return "".join(out)


def code_chip(text: str) -> str:
    if len(text) <= SHORT_CODE:
        return f'<code class="nw">{html.escape(text)}</code>'
    return f"<code>{code_pieces(text)}</code>"


def chips_in(fragment: str) -> str:
    """The same treatment for the <code> of a figure written as HTML."""
    return re.sub(r"<code>([^<]*)</code>", lambda m: code_chip(html.unescape(m.group(1))), fragment)


# ── Options ──
#
# Every --option of the parameters page has a row there with an id, and a
# chip of code naming one anywhere else in the docs links to it.

PARAMETERS = "user-guide/parameters.md"


@functools.lru_cache(maxsize=None)
def parameter_rows() -> list[dict]:
    """The --option rows of the parameters page: flag, section, cells."""
    rows, section, head = [], "", []
    for line in (DOCS / PARAMETERS).read_text(encoding="utf-8").splitlines():
        heading = re.match(r"^##\s+(.*)$", line)
        if heading:
            section = heading.group(1).strip()
            continue
        if not line.startswith("|"):
            head = []
            continue
        cells = [c.strip() for c in re.split(r"(?<!\\)\|", line.strip().strip("|"))]
        if not head:
            head = [plain(c) for c in cells]
            continue
        flag = re.fullmatch(r"`(--[a-z][a-z0-9_]*)`", cells[0])
        if flag:
            rows.append({"flag": flag.group(1), "section": section, "cells": dict(zip(head, (plain(c) for c in cells)))})
    return rows


def options() -> set[str]:
    return {row["flag"] for row in parameter_rows()}


def switches() -> set[str]:
    """Options that take no value: what follows one on a command line is not its value."""
    return {row["flag"] for row in parameter_rows() if row["cells"].get("Type") == "boolean"} | {"--rm", "--help", "--version"}


# Nextflow's own options that take a value.
NF_VALUED = {"-profile", "-work-dir", "-c", "-r"}


# ── Code blocks ──

SHELL = {"bash", "sh", "console", "shell"}
TOKEN = re.compile(r"""(?:[^\s'"]|'[^']*'|"(?:\\.|[^"\\])*")+|\S+""")
LONG_TOKEN = 36


def open_quote(line: str, quote: str | None) -> str | None:
    """The quote still open at the end of a line of shell, given the one open
    at its start."""
    i = 0
    while i < len(line):
        ch = line[i]
        if quote is None:
            if ch == "\\":
                i += 2
                continue
            if ch == "#" and (i == 0 or line[i - 1].isspace()):
                return None
            if ch in "'\"":
                quote = ch
        elif ch == quote:
            quote = None
        elif ch == "\\" and quote == '"':
            i += 2
            continue
        i += 1
    return quote


# A word of a command longer than this may come apart on a phone: a path
# after a /, anything else (a commit hash, a quoted sentence) where it must.
SHELL_LONG = 26


def shell_token(token: str, kind: str) -> str:
    """One word of a command, which never wraps inside itself unless it is
    too long for a phone's line: a path then breaks after a /, and anything
    else at its spaces or, failing those, anywhere."""
    if len(token) > SHELL_LONG and "/" in token and not token.startswith(("'", '"')):
        # After each /, and in a file name too long for a line after _ or .
        pieces = [q for p in re.findall(r"[^/]*/+|[^/]+", token)
                  for q in (re.findall(r"[^_.]*[_.]+|[^_.]+", p) if len(p) > SHELL_LONG else [p])]
        inner = "<wbr>".join(f'<span class="nw">{html.escape(p)}</span>' for p in pieces)
        return f'<span class="{kind} path">{inner}</span>'
    if len(token) > SHELL_LONG:
        return f'<span class="{kind} brk">{html.escape(token)}</span>'
    return f'<span class="{kind}">{html.escape(token)}</span>'


def shell_line(line: str) -> str:
    """A line of shell as words: flags (.f), their values (.v), other words
    (.t) and a comment (.c). The spaces stay between the words, so a long
    line wraps there and only there. A line's indent stays with its first
    word, and a closing backslash with its last, so neither ends up on a line
    of its own on a phone."""
    units, at, tokens = [], 0, list(TOKEN.finditer(line))
    k = 0
    while k < len(tokens):
        m = tokens[k]
        gap = html.escape(line[at:m.start()])
        word = m.group(0)
        if word.startswith("#"):
            units.append([gap, f'<span class="c">{html.escape(line[m.start():])}</span>'])
            at = len(line)
            break
        if re.match(r"^--?[A-Za-z]", word):
            nxt = tokens[k + 1] if k + 1 < len(tokens) else None
            valued = (
                (word.startswith("--") and "=" not in word and word not in switches()) or word in NF_VALUED or word.startswith("-process.")
            )
            if valued and nxt and not re.match(r"^(-|#|\\$|\||&&|;|>)", nxt.group(0)):
                value = nxt.group(0)
                pair = shell_token(word, "f") + html.escape(line[m.end():nxt.start()]) + shell_token(value, "v")
                units.append([gap, f'<span class="kv">{pair}</span>' if len(word) + len(value) < SHELL_LONG else pair])
                at = nxt.end()
                k += 2
                continue
            units.append([gap, shell_token(word, "f")])
        else:
            units.append([gap, shell_token(word, "t")])
        at = m.end()
        k += 1
    if len(units) > 1 and tokens and tokens[-1].group(0) == "\\" and at == tokens[-1].end():
        gap, slash = units.pop()
        # Only the last word goes with it: of a flag and its value, the value.
        body = units[-1][1]
        last = re.search(r'(<span class="v[^"]*">.*</span>)$', body) if body.startswith('<span class="f') else None
        head, tail = (body[: last.start()], last.group(1)) if last and not body.startswith('<span class="kv') else ("", body)
        if not re.search(r'class="[^"]*\b(path|brk)\b', tail):
            units[-1][1] = f'{head}<span class="nw">{tail}{gap}{slash}</span>'
        else:
            units.append([gap, slash])
    if units and units[0][0].strip() == "" and units[0][0]:
        # The indent goes with the first word: of a flag and its value, the flag.
        body = units[0][1]
        first = re.match(r'(<span class="f">[^<]*</span>)(.+)$', body, re.S)
        if first:
            units[0] = ["", f'<span class="nw">{units[0][0]}{first.group(1)}</span>{first.group(2)}']
        elif not re.match(r'<span class="[^"]*\b(path|brk)\b', body):
            units[0] = ["", f'<span class="nw">{units[0][0]}{body}</span>']
    return "".join(gap + body for gap, body in units) + html.escape(line[at:])


def shell_html(lines: list[str], prompts: bool) -> str:
    """Each line of a shell block as a .ln, which wraps with a hanging
    indent. A line that starts a command gets a $ prompt, which Copy leaves
    out; a line that carries one on, after a backslash or inside a quote or
    a loop, does not."""
    out = []
    carried = False
    quote = None
    depth = 0
    for n, line in enumerate(lines):
        stripped = line.strip()
        indent = len(line) - len(line.lstrip(" "))
        starts = prompts and stripped and not carried and quote is None and depth == 0 and not stripped.startswith("#")
        prompt = '<span class="p">$ </span>' if starts else ""
        # Inside a quote that spans lines is another program's text (an R
        # call, say): its words still wrap only between them, uncoloured.
        body = shell_line(line) if quote is None else re.sub(r"\S+", lambda m: shell_token(m.group(0), "t"), line)
        style = f' style="--i:{indent}"' if indent else ""
        end = "\n" if n + 1 < len(lines) else ""
        out.append(f'<span class="ln"{style}>{prompt}{body}{end}</span>')
        quote = open_quote(line, quote)
        carried = stripped.endswith("\\")
        if quote is None and not carried and stripped:
            first = stripped.split()[0]
            closes = re.search(r"\b(done|fi|esac)\b", stripped)
            if first in ("for", "while", "until", "if", "case") and not closes:
                depth += 1
            elif first in ("done", "fi", "esac"):
                depth = max(depth - 1, 0)
    return "".join(out)


def tsv_html(lines: list[str]) -> str:
    """A tab-separated file, with each tab drawn. Copy takes the real tabs."""
    return "\n".join(
        '<span class="tsv-tab">\t</span>'.join(html.escape(cell) for cell in line.split("\t")) for line in lines
    )


CHANGELOG_EMOJI = re.compile("[\U0001F000-\U0001FAFF☀-➿⬀-⯿️‍]")
MONTHS = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]


def version_heading(text: str) -> tuple[str, str, str]:
    """A version of the changelog as the page names it: '[1.0.1] - 2026-07-25'
    becomes ('1.0.1, 25 July 2026', 'v1-0-1', '1.0.1')."""
    if re.fullmatch(r"\[?Unreleased\]?", text.strip(), re.I):
        return "On main, not yet released", "unreleased", "Unreleased"
    m = re.match(r"^\[?([\d.]+)\]?(?:\s*-\s*(\d{4})(?:-(\d{2})-(\d{2}))?)?", text.strip())
    if not m:
        return text, slugify(text, {}), text
    version, year, month, day = m.groups()
    when = f"{int(day)} {MONTHS[int(month) - 1]} {year}" if day else (year or "")
    return (f"{version}, {when}" if when else version), "v" + version.replace(".", "-"), version


# The room a table has on a 360 px phone (326 px), with a margin for an
# estimate that runs about 10% wide: measured, tables estimated at 341 fit
# and those at 391 overflow. Then what a character takes there, code or
# text, and a cell's padding.
PHONE_TABLE = 360
PHONE_CHAR = 7.8
PHONE_CELL = 20


def phone_width(head: list[str], body: list[list[str]]) -> float:
    """About the narrowest a table can be: each column as wide as its widest
    word or unbreakable chip, the first column's chips whole."""
    columns = [0.0] * len(head)
    for row in [head] + body:
        for k, cell in enumerate(row[: len(head)]):
            units = []
            for m in re.finditer(r"`([^`]+)`", cell):
                code = m.group(1).strip()
                if k == 0 or len(code) <= SHORT_CODE:
                    units.append(len(code) + 1.5)
                else:
                    units.append(max(len(p) for t in code.split() for p in re.findall(r"[^/_.]*[/_.]+|[^/_.]+", t)) + 1.5)
            units += [len(w) for w in plain(re.sub(r"`[^`]+`", " ", cell)).split()]
            columns[k] = max([columns[k]] + units)
    return sum(c * PHONE_CHAR + PHONE_CELL for c in columns)


def first_sentence(text: str, limit: int = 160) -> str:
    m = re.match(r"(.+?[.:])(\s|$)", text)
    cut = m.group(1) if m else text
    return cut if len(cut) <= limit else cut[:limit].rsplit(" ", 1)[0] + "..."


class Renderer:
    def __init__(self, source: str):
        self.source = source  # markdown path relative to docs/, e.g. user-guide/outputs.md
        self.root = "../" * source.count("/")
        self.changelog = source == CHANGELOG[0]
        self.version: tuple[str, str] | None = None  # the changelog's version being read: (slug, label)
        self.seen: dict[str, int] = {}
        self.gh_seen: dict[str, int] = {}  # GitHub's anchors, for the changelog's aliases
        self.legacy = {new: old for old, new in LEGACY_IDS.get(source, {}).items()}
        self.headings: list[tuple[int, str, str]] = []  # (level, text, slug)
        self.sections: list[dict] = []  # for the search index
        self.extra: list[dict] = []  # search entries of their own: callouts
        self.figures: list[str] = []
        self.words = 0  # of running text, for the reading time

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

    def inline(self, text: str, links: bool = True, whole: bool = False) -> str:
        """An inline span as HTML. With links, a chip of code that names an
        option of the parameters page links to its row there; whole keeps
        a chip on one line, as the name of a table's row, unless it is longer
        than a phone's line holds (a long file name then breaks after a _ or .)."""
        keep: list[str] = []
        unlinked: dict[str, str] = {}

        def stash(fragment: str) -> str:
            keep.append(fragment)
            return f"\x00{len(keep) - 1}\x00"

        def code(m: re.Match) -> str:
            raw = m.group(2).strip()
            chip = f'<code class="nw">{html.escape(raw)}</code>' if whole and len(raw) <= LONG_TOKEN else code_chip(raw)
            flag = raw.split(" ")[0]
            if links and flag in options() and self.source != PARAMETERS:
                marker = stash(f'<a class="opt-link" href="{self.root}user-guide/parameters.html#opt-{flag[2:]}">{chip}</a>')
                unlinked[marker] = chip
                return marker
            return stash(chip)

        text = re.sub(r"(`+)(.+?)\1", code, text)
        text = re.sub(r"\\([\\`*_{}\[\]()#+\-.!|<>])", lambda m: stash(html.escape(m.group(1))), text)
        text = re.sub(r"<(https?://[^>\s]+)>", lambda m: stash(f'<a href="{html.escape(m.group(1))}">{html.escape(m.group(1))}</a>'), text)
        text = INLINE_TAGS.sub(lambda m: stash(m.group(0)), text)

        def link(m: re.Match) -> str:
            label, target = m.group(1), m.group(2).strip().split(" ")[0]
            # A link's own text keeps its chips plain: no link inside a link.
            for marker, chip in unlinked.items():
                label = label.replace(marker, stash(chip))
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
        # A text block right after a shell block is what that command prints.
        after_shell = False
        while i < len(lines):
            line = lines[i]
            stripped = line.strip()
            if not stripped:
                i += 1
                continue
            shell = False
            figure = FIGURE_MARK.match(stripped)
            if figure:
                name = figure.group(1)
                if name in DERIVED:
                    out.append(chips_in(DERIVED[name](self, lines)))
                elif name in FIGURES:
                    # Only a drawing that figures.js draws needs it loaded.
                    if FIGURE_HOOK.search(FIGURES[name]):
                        self.figures.append(name)
                    out.append(chips_in(FIGURES[name]))
                else:
                    raise SystemExit(f"docs/{self.source}: no figure called {name!r}")
                i += 1
                after_shell = False
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
                after_shell = False
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
                out.append(self.code(info, body, after_shell))
                after_shell = (info.split() or [""])[0] in SHELL and "title=" not in info
                continue
            after_shell = False
            heading = re.match(r"^(#{1,6})\s+(.*?)\s*#*$", line)
            if heading:
                out.append(self.heading(len(heading.group(1)), heading.group(2)))
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
            self.words += len(plain(text).split())
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

    def heading(self, level: int, text: str) -> str:
        """A heading with its id. Below the h1 it carries a # link, outside the
        heading itself so that the heading's name is only its text."""
        label = plain(text)
        github = slugify(label, self.gh_seen) if self.changelog else ""
        if self.changelog and level == 2:
            label, slug, short = version_heading(label)
            self.seen[slug] = 0
            self.version = (slug, short)
            rendered = html.escape(label)
        elif self.changelog and level > 2:
            label = CHANGELOG_EMOJI.sub("", label).strip()
            slug = slugify(f"{self.version[0]} {label}" if self.version else label, self.seen)
            rendered = html.escape(label)
        else:
            slug = slugify(label, self.seen)
            rendered = self.inline(text, links=False)
        self.headings.append((level, label, slug))
        section = label
        if self.changelog and level > 2 and self.version:
            section = f"{self.version[1]} / {label}"
        self.sections.append({"heading": section, "slug": slug, "level": level, "text": [], "code": []})
        aliases = ""
        if github and github != slug:
            aliases += anchor_alias(github, slug)
        if slug in self.legacy:
            aliases += anchor_alias(self.legacy.pop(slug), slug)
        if level == 1:
            return f'{aliases}<h1 id="{slug}">{rendered}</h1>'
        return (
            f'{aliases}<div class="hx hx-{level}"><h{level} id="{slug}">{rendered}</h{level}>'
            f'<a class="anchor" href="#{slug}" aria-label="Link to: {html.escape(label)}">#</a></div>'
        )

    def note(self, text: str) -> None:
        if self.sections:
            self.sections[-1]["text"].append(plain(text))

    def note_code(self, text: str) -> None:
        if self.sections:
            self.sections[-1]["code"].append(text)

    def quote(self, quoted: list[str]) -> str:
        """A block quote, or one of GitHub's alerts (> [!WARNING]) as a callout
        whose first bold line, if it has one, is its title. A callout says its
        kind in words, not only in its colour, and is found by search."""
        first = next((q.strip() for q in quoted if q.strip()), "")
        alert = re.match(r"^\[!(\w+)\]$", first)
        if not alert or alert.group(1).upper() not in ALERTS:
            return f"<blockquote>{self.blocks(quoted)}</blockquote>"
        kind = alert.group(1).upper()
        body = quoted[quoted.index(next(q for q in quoted if q.strip())) + 1 :]
        while body and not body[0].strip():
            body = body[1:]
        title = ALERTS[kind]
        label = ""
        if body and re.match(r"^\*\*(.+)\*\*$", body[0].strip()):
            title = self.inline(body[0].strip()[2:-2])
            self.note(plain(body[0].strip()))
            label = f'<span class="callout-kind">{ALERTS[kind]}</span>'
            body = body[1:]
        name = html.unescape(re.sub(r"<[^>]+>", "", title))
        cid = slugify("note " + name, self.seen)
        if label:
            self.extra.append({
                "section": name,
                "slug": cid,
                # The whole callout is in its section's entry; this one finds it by its title.
                "text": first_sentence(plain(" ".join(q.strip() for q in body if q.strip() and not q.strip().startswith("```")))),
            })
        return (
            f'<div class="callout callout-{kind.lower()}" id="{cid}"><p class="callout-title">{label}{title}</p>'
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
            # Each option of the parameters page has a search entry of its own.
            if not (self.source == PARAMETERS and re.fullmatch(r"`--[a-z][a-z0-9_]*`", row[0])):
                self.note(" ".join(row))
        labels = [plain(c) for c in head]

        def cls(k: int, extra: str = "") -> str:
            names = ["num"] if k < len(aligns) and aligns[k].endswith(":") and not aligns[k].startswith(":") else []
            names += [extra] if extra else []
            return f' class="{" ".join(names)}"' if names else ""

        def cell(k: int, c: str) -> str:
            label = f' data-label="{html.escape(labels[k])}"' if k < len(labels) else ""
            # A type is a word to read, not code to type.
            if k < len(labels) and labels[k].lower() == "type":
                return f"<td{cls(k, 'type')}{label}>{html.escape(plain(c))}</td>"
            return f"<td{cls(k)}{label}>{self.inline(c, whole=k == 0)}</td>"

        # A table of four columns or more, or one that cannot fit a phone
        # without a chip of code splitting, becomes a list of cards there.
        cards = len(head) >= 4 or phone_width(head, body) > PHONE_TABLE
        named = f' aria-labelledby="{self.sections[-1]["slug"]}"' if self.sections else ""
        wrap = "table-wrap cards-wrap" if cards else "table-wrap"
        table_class = ' class="cards"' if cards else ""
        # site.js makes it a named, focusable region only while it scrolls.
        out = [f'<div class="{wrap}"{named}><table{table_class}><thead><tr>']
        out += [f"<th{cls(k)}>{self.inline(c)}</th>" for k, c in enumerate(head)]
        out.append("</tr></thead><tbody>")
        for row in body:
            flag = re.fullmatch(r"`(--[a-z][a-z0-9_]*)`", row[0]) if row else None
            rid = f' id="opt-{flag.group(1)[2:]}"' if flag else ""
            out.append(f"<tr{rid}>" + "".join(cell(k, c) for k, c in enumerate(row)) + "</tr>")
        out.append("</tbody></table></div>")
        return "".join(out)

    def code(self, info: str, body: list[str], after_shell: bool = False) -> str:
        """A fenced block. A shell block is a terminal, its lines wrapping
        between words with a hanging indent; the text block right after one is
        what it prints; a file keeps its layout and says its name."""
        self.last_code = body
        title = re.search(r'title="([^"]*)"', info)
        lang = info.split()[0] if info and not info.startswith("title=") else ""
        data_lang = f' data-lang="{html.escape(lang)}"' if lang else ""
        self.note_code("\n".join(body))
        classes = ["code"]
        head = ""
        if title:
            head = f'<div class="code-title fname">{html.escape(title.group(1))}</div>'
        if lang in SHELL:
            classes.append("cmd")
            text = shell_html(body, prompts=not title)
            if not title:
                head = '<div class="code-title">Terminal</div>'
        elif lang == "text" and after_shell and not title:
            classes.append("out")
            head = '<div class="code-title">Output</div>'
            text = html.escape("\n".join(body))
        elif lang == "tsv":
            classes.append("tsv")
            text = tsv_html(body)
            if not title:
                head = '<div class="code-title">Tab separated</div>'
        else:
            text = html.escape("\n".join(body))
            # A tree of files with a comment after each: on a phone each
            # comment goes under its file, from where the name starts.
            if sum(1 for line in body if re.search(r"\S {2,}# ", line)) >= 2:
                classes.append("tree")
                rows = []
                for n, line in enumerate(body):
                    end = "\n" if n + 1 < len(body) else ""
                    m = re.match(r"^(.*?\S)( {2,})(# .*)$", line)
                    if not m:
                        rows.append(html.escape(line) + end)
                        continue
                    col = len(re.match(r"^[│├└─ ]*", m.group(1)).group(0))
                    # The line's end goes inside the comment: as a block on a
                    # phone, it would otherwise leave an empty line after it.
                    rows.append(
                        html.escape(m.group(1))
                        + f'<span class="tc" style="--c:{col}">{html.escape(m.group(2) + m.group(3))}{end}</span>'
                    )
                text = "".join(rows)
        if not head:
            classes.append("bare")
        return f'<div class="{" ".join(classes)}"{data_lang}>{head}<pre>{text}</pre></div>'

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


def sheet_error(title: str, values: dict[str, str]) -> list[str]:
    """The box main.nf stops with when it reads a sample sheet that has a
    mistake, before any task runs: the lines of the error whose box starts
    with title, with its ${...} filled in from values. The box's frame is
    left out, so the words wrap on a phone and nothing else does."""
    code = Path("main.nf").read_text(encoding="utf-8")
    for body in re.findall(r'error """\n(.*?)"""', code, re.S):
        lines = textwrap.dedent(body).strip().splitlines()
        if len(lines) > 1 and lines[1].strip("║ ").startswith(title):
            out = []
            for line in lines:
                if line.startswith(("╔", "╚")):
                    continue
                for expr, value in values.items():
                    line = line.replace("${" + expr + "}", value)
                if "${" in line:
                    raise SystemExit(f"tools/build_docs.py: main.nf's '{title}' box has {line.strip()}, which sheet_error() cannot fill")
                out.append(re.sub(r"^║ {0,2}", "", line).replace("\\\\t", "\\t").rstrip())
            return out
    raise SystemExit(f"tools/build_docs.py: main.nf no longer stops on a sample sheet with '{title}'")


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
# The scheduler's settings come from nextflow.config's profiles and the caps it
# sets after them, and the launcher's from the job script the HPC page gives, so
# neither can drift.


def profile_settings(name: str) -> dict[str, str]:
    config = Path("nextflow.config").read_text(encoding="utf-8")
    block = re.search(rf"^    {name}\s*\{{(.*?)^    \}}", config[config.index("profiles {"):], re.M | re.S)
    if not block:
        raise SystemExit(f"tools/build_docs.py: nextflow.config has no profile {name}")
    return dict(re.findall(r"^\s+([\w.]+)\s*=\s*'?([^'\n]+?)'?\s*(?://.*)?$", block.group(1), re.M))


def caps() -> dict[str, dict[str, str]]:
    """The caps nextflow.config sets after its profiles, by what decides them:
    {"highmem", "test", "scheduler", "machine"} -> {"cpus", "memory", "time"}.
    Each is a line `params.max_X = params.highmem ? a : params.test_data ? b :
    ... : otherwise`, and the first flag that is on picks the value."""
    config = Path("nextflow.config").read_text(encoding="utf-8")
    start = config.index("\nprofiles {")
    after = config[start + re.search(r"^\}", config[start + 1:], re.M).end() + 1:]
    flags = {"highmem": "highmem", "test": "test_data", "scheduler": "scheduler", "machine": None}
    found = {kind: {} for kind in flags}
    for cap, expr in re.findall(r"^params\.max_(cpus|memory|time)\s*=\s*(.+)$", after, re.M):
        branches = dict(re.findall(r"params\.(\w+)\s*\?\s*([^:]+?)\s*:", expr))
        otherwise = expr.rsplit(":", 1)[-1].strip()
        for kind, flag in flags.items():
            found[kind][cap] = branches.get(flag, otherwise).strip("'")
    if any(len(values) != 3 for values in found.values()):
        raise SystemExit("tools/build_docs.py: nextflow.config no longer sets the caps the way caps() reads them")
    return found


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
    cluster = caps()["scheduler"]
    memory = cluster["memory"].replace(".GB", " GB")
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
        f"<p>One job per task, at most {cluster['cpus']} CPUs and {memory} each, in the AMAP image, "
        + ("working on the node's own disk and copying back only the results.</p>" if slurm.get("params.scratch") == "true" else "</p>")
        + '<div class="hpc-tasks"><span>per sample: Bakta, miniprot, the correction, duplications</span>'
        "<span>once for all: eggNOG, the pangenome, the genomes compared, the report</span></div></div>"
    )
    # The singularity profile keeps the image under --outdir, and that wins
    # over NXF_SINGULARITY_CACHEDIR: one shared copy needs singularity.cacheDir
    # in a config passed with -c, as deployment/hpc.md says.
    cache = profile_settings("singularity").get("singularity.cacheDir", "")
    if "params.outdir" not in cache:
        raise SystemExit("tools/build_docs.py: the singularity profile no longer keeps its image under --outdir; update hpc_figure()")
    shared = (
        '<div class="hpc-box hpc-shared"><p class="hpc-title">Shared storage</p><div class="hpc-places">'
        '<span><code class="hpc-key">-work-dir</code>every task\'s files, which <code>-resume</code> needs</span>'
        '<span><code class="hpc-key">singularity.cacheDir</code>the image, converted once, from a <code>-c</code> config</span>'
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
    # Each withName block of the top-level process { }, not a profile's: the
    # slurm profile has one for Bakta of its own, with its QoS.
    scope = re.search(r"^process \{\n(.*?)^\}", config, re.S | re.M).group(1)
    blocks = dict(re.findall(r"withName: '(\w+)' \{(.*?)\n    \}", scope, re.S))
    bakta = re.search(r"memory\s*=\s*\{ (\d+)\.GB \* task\.attempt \}", blocks.get("BAKTA_ANNOTATOR", ""))
    terminate = [name for name, body in blocks.items() if "errorStrategy = 'terminate'" in body]
    eggnog = re.search(r"errorStrategy = \{ task\.exitStatus in \[[\d, ]+\] \? 'retry' : 'terminate' \}", blocks.get("ANNOT_EGGNOGMAP_BATCH", ""))
    # The pangenome runs once for the whole run but keeps the default rule, so
    # it is set aside like a step run per sample (main.nf goes on without it).
    pangene = blocks.get("PANGENE")
    if not (rule and bakta and eggnog and pangene is not None and "errorStrategy" not in pangene):
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
        + '<p class="fail-except">The steps run once for all the samples are retried too, and then stop the run, except the '
        + "pangenome: like a step run per sample, it is set aside after two more tries, and the reports go on without it. Some stop it "
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
    # A path comes apart only after a /, never inside a file name.
    path = lambda b: "<wbr>".join(f'<span class="nw">{html.escape(p)}</span>' for p in re.findall(r"[^/]*/+|[^/]+", b))
    body = "".join(f'<tr><td>{html.escape(a)}</td>{tab}<td class="sheet-path">{path(b)}</td></tr>' for a, b in rows)
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
    # The two mistakes main.nf catches as it reads the sheet, before any task
    # starts: spaces for tabs make the header one column, and the first path
    # that does not exist stops it.
    spaced_header = "".join(f"{a:<12}" for a in SHEET[0]).rstrip()
    first_id, first_path = SHEET[1]
    checks = [
        ("With spaces where the tabs go", sheet_error("Invalid TSV format!", {"row.keySet().join('  ')": spaced_header})),
        ("With tabs, and paths that do not exist yet",
         sheet_error("FASTA file not found!", {"row.sampleId": first_id, "row.fasta": first_path})),
    ]
    blocks = "".join(
        f'<div class="sheet-check"><p class="sheet-when">{when}, <code>nextflow run</code> stops before any task with:</p>'
        f'<pre>{html.escape(chr(10).join(lines))}</pre></div>'
        for when, lines in checks
    )
    return (
        '<figure class="fig sheet-fig"><div class="sheet"><div class="sheet-file"><p class="sheet-name">samples.tsv'
        '<span class="sheet-n sheet-n-tab">2</span></p>'
        f"{table}</div>{notes}</div>"
        f'<div class="sheet-checks">{blocks}</div>'
        "<figcaption>The sample sheet, with its tabs drawn. Below it, what a run prints for the two mistakes it catches "
        "most, <code>--validate_only</code> or not: both messages are <code>main.nf</code>'s own, filled in for this sheet. "
        "They stop the run before the checks start, so no <code>validation.log</code> is written.</figcaption></figure>"
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
        "<figcaption>Where a genome has to sit, with the limits read off <code>bin/validate_inputs.py</code>, and where six of "
        "the seven genomes of the accuracy benchmark sit: inside every range, with no warning. The seventh, CP110619, is "
        "another assembly of H37Rv and is not measured here.</figcaption></figure>"
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


def releases() -> list[tuple[str, bool]]:
    """The versions of CHANGELOG.md, newest first: (version, has notes). The
    first is 'Unreleased' while main is ahead of the last release."""
    out = []
    for block in re.split(r"^## ", CHANGELOG[1].read_text(encoding="utf-8"), flags=re.M)[1:]:
        name = re.match(r"\[?([\w.]+)\]?", block).group(1)
        out.append((name, bool(block.split("\n", 1)[1].strip()) if "\n" in block else False))
    return out


def released() -> str:
    """The last version released."""
    return next(name for name, _ in releases() if name[0].isdigit())


def header(root: str, current: str, this_page: str = "") -> str:
    """The site's header. The link to the section a page is in is current;
    it is the current page only when it leads to this page itself."""
    def nav(target: str, label: str, key: str, optional: bool = False) -> str:
        cls = ' class="nav-optional"' if optional else ""
        to_here = this_page and target == f"{root}{this_page}"
        here = ' aria-current="page"' if to_here else ' aria-current="true"' if key == current else ""
        return f'      <a{cls} href="{target}"{here}>{label}</a>\n'

    return (
        '<header class="top">\n  <div class="wrap">\n'
        f'    <a class="brand" href="{root}index.html" aria-label="AMAP home">\n'
        f'      <img class="logo-dark" src="{root}assets/amap-logo.png" alt="" width="57" height="40">'
        f'<img class="logo-light" src="{root}assets/amap-logo-light.png" alt="" width="57" height="40"><span>AMAP</span>\n'
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
    version = released()
    return (
        '<footer class="foot">\n  <div class="wrap">\n    <div>\n      <h3>Cite</h3>\n'
        f"      <blockquote>Ruiz-Rodriguez, P., &amp; Coscolla, M. (2026). AMAP (Version {version}) [Computer software].</blockquote>\n"
        '      <p class="doi">Software: <a href="https://doi.org/10.5281/zenodo.8024446">doi:10.5281/zenodo.8024446</a></p>\n    </div>\n'
        "    <div>\n      <h3>Learn more</h3>\n      <ul>\n"
        f'        <li><a href="{root}getting-started/tutorial.html">Tutorial</a></li>\n'
        f'        <li><a href="{root}{INDEX}">Documentation</a></li>\n'
        f'        <li><a href="{root}how-it-works.html">How it works</a></li>\n'
        f'        <li><a href="{root}user-guide/parameters.html">Parameters</a></li>\n'
        f'        <li><a href="{root}about/changelog.html">Changelog</a></li>\n'
        "      </ul>\n    </div>\n    <div>\n      <h3>Project</h3>\n      <ul>\n"
        f'        <li><a href="{REPO}">GitHub</a></li>\n'
        f'        <li><a href="{REPO}/issues">Report a bug</a></li>\n'
        '        <li><a href="https://hub.docker.com/r/paururo/amap">Docker image</a></li>\n'
        "      </ul>\n    </div>\n  </div>\n"
        f'  <div class="wrap foot-base"><span>PathoGenOmics Lab</span><span>GPL-3.0</span><span>AMAP {version}</span></div>\n'
        "</footer>\n"
    )


# Elements that close themselves: no </...> follows them.
VOID = {"br", "wbr", "img", "hr", "input", "source", "col"}


def lead_in(caption: str) -> str:
    """A caption with its first sentence, or what comes before its first
    colon, set bold, the way a paper sets a figure's: the rest explains it. A
    caption of one sentence stays as it is. A full stop after an initial (the
    M. of M. bovis) ends nothing, and neither does one inside an element."""
    depth = 0
    i = 0
    while i < len(caption):
        c = caption[i]
        if c == "<":
            j = caption.index(">", i)
            tag = re.match(r"</?([a-zA-Z0-9]+)", caption[i : j + 1])
            name = tag.group(1).lower() if tag else ""
            if caption[i + 1] == "/":
                depth -= 1
            elif name not in VOID and not caption[j - 1] == "/":
                depth += 1
            i = j + 1
            continue
        if depth == 0 and c in ".:" and caption[i + 1 : i + 2] == " ":
            word = re.search(r"([^\s>]*)$", caption[:i]).group(1)
            if c == "." and len(word) <= 1:
                i += 1
                continue
            if caption[i + 1 :].strip():
                return f'<b class="cap-lead">{caption[: i + 1]}</b>{caption[i + 1 :]}'
            return caption
        i += 1
    return caption


def lead_ins(body: str) -> str:
    """Every figure caption of a page, with its lead-in."""
    return re.sub(r"<figcaption>(.*?)</figcaption>", lambda m: f"<figcaption>{lead_in(m.group(1))}</figcaption>", body, flags=re.S)


def section_id(section: str) -> str:
    """The id of a section of the docs on the Docs tab, e.g. 'user-guide'."""
    return slugify(section, {})


def page_contents(headings: list[tuple[int, str, str]], levels: tuple[int, ...]) -> str:
    """The page's own contents for the sidebar, nested under its entry: each
    h2, and the h3s under it in a list of their own, always shown, so the
    list stays still as the reader scrolls (site.js marks the one being
    read)."""
    groups: list[tuple[str, str, list[tuple[str, str]]]] = []
    for level, text, slug in headings:
        if level not in levels:
            continue
        if level == 3 and groups:
            groups[-1][2].append((text, slug))
        else:
            groups.append((text, slug, []))
    if not groups:
        return ""
    out = ['\n        <ul class="page-toc" aria-label="On this page">\n']
    for text, slug, subs in groups:
        out.append(f'          <li class="toc-2"><a href="#{slug}">{html.escape(text)}</a>')
        if subs:
            out.append(
                "<ul>"
                + "".join(f'<li class="toc-3"><a href="#{s}">{html.escape(t)}</a></li>' for t, s in subs)
                + "</ul>"
            )
        out.append("</li>\n")
    out.append("        </ul>\n      ")
    return "".join(out)


def sidebar(root: str, current: str | None, contents: str = "") -> str:
    """The page list, its search, and under the current page its contents:
    the one place a page's sections are listed where there is a sidebar."""
    out = [
        '<div class="docs-nav">\n',
        '  <div class="docs-search">\n'
        '    <input type="search" placeholder="Search the docs" aria-label="Search the documentation" data-docs-search>\n'
        '    <div class="docs-results" data-docs-results hidden></div>\n  </div>\n',
        '  <nav aria-label="Documentation">\n  <details class="docs-menu" open>\n    <summary>Pages</summary>\n',
    ]
    section = None
    for sec, path, title, _ in PAGES:
        if sec != section:
            if section is not None:
                out.append("    </ul>\n")
            label = f"nav-{section_id(sec)}"
            out.append(f'    <p class="docs-section" id="{label}">{sec}</p>\n    <ul aria-labelledby="{label}">\n')
            section = sec
        if path == current:
            out.append(f'      <li><a href="{root}{href(path)}" aria-current="page">{title}</a>{contents}</li>\n')
        else:
            out.append(f'      <li><a href="{root}{href(path)}">{title}</a></li>\n')
    # On a phone the page list is folded before the first paint, so the page
    # itself comes first; with no script, or on a wider screen, it stays open.
    out.append(
        "    </ul>\n  </details>\n"
        '  <script>if (matchMedia("(max-width: 860px)").matches) document.currentScript.previousElementSibling.open = false;</script>\n'
        "  </nav>\n</div>\n"
    )
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
        '<meta name="theme-color" content="#0a0e0d">\n'
        f'<link rel="apple-touch-icon" href="{root}assets/touch-icon.png">\n'
    )


def page(root: str, title: str, description: str, body: str, figures: list[str] | None = None, path: str = "",
         body_class: str = "docs", note: str = "", tail: str = "") -> str:
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
        f'{THEME_BOOT}\n</head>\n<body class="{body_class}">\n'
        '<a class="skip" href="#content">Skip to the content</a>\n'
        + (note or "<!-- Generated by tools/build_docs.py from the markdown next to this file. Edit that, then run it. -->\n")
        + lead_ins(body)
        # The search index is not loaded here: site.js fetches it the first
        # time the reader goes to search.
        + f'<script src="{root}assets/site.js"></script>\n'
        + (f'<script src="{root}assets/figure-data.js"></script>\n' if figures else "")
        + "".join(f'<script src="{root}{d}"></script>\n' for d in data)
        + (f'<script src="{root}assets/figures.js"></script>\n' if figures else "")
        + tail
        + "</body>\n</html>\n"
    )


def source_lines(path: str) -> list[str]:
    if path == CHANGELOG[0]:
        return CHANGELOG[1].read_text(encoding="utf-8").splitlines()
    return strip_front_matter((DOCS / path).read_text(encoding="utf-8").splitlines())


def version_line(root: str) -> str:
    """Which AMAP the pages describe: main, while it is ahead of the last release."""
    ahead = any(name.lower() == "unreleased" and notes for name, notes in releases())
    version = released()
    if ahead:
        return (f'<p class="ver">Describes AMAP on main, ahead of {version}: '
                f'<a href="{root}about/changelog.html#unreleased">what changed</a></p>')
    return f'<p class="ver">Describes <a href="{root}about/changelog.html#v{version.replace(".", "-")}">AMAP {version}</a></p>'


def doc_page(position: int) -> tuple[Path, str, list[dict]]:
    section, path, title, description = PAGES[position]
    renderer = Renderer(path)
    content = renderer.blocks(source_lines(path))
    if renderer.legacy:
        raise SystemExit(f"tools/build_docs.py: {path} has no heading {sorted(renderer.legacy)} for LEGACY_IDS to lead to")
    if renderer.changelog:
        content = content.replace("documented in this file", "documented on this page")
    root = "../" * path.count("/")
    # The changelog's contents are its versions; its categories repeat in each.
    levels = (2,) if renderer.changelog else (2, 3)
    contents = page_contents(renderer.headings, levels)
    h2s = [(text, slug) for level, text, slug in renderer.headings if level == 2]
    # Under the h1: how long the page takes to read, and where there is no
    # sidebar to list them (a phone), the contents themselves, folded.
    minutes = max(1, round(renderer.words / 200))
    # Under the h1, what the page is for: its line from PAGES, the one the
    # Docs tab and a shared link show.
    under = f'<p class="standfirst">{html.escape(description)}</p>\n<p class="doc-meta">{minutes} min read</p>\n'
    if h2s:
        under += (
            '<details class="toc-inline"><summary>On this page</summary><ul>'
            + "".join(f'<li><a href="#{slug}">{html.escape(text)}</a></li>' for text, slug in h2s)
            + "</ul></details>\n"
        )
    end = content.find("</h1>")
    if end >= 0:
        content = content[: end + 5] + "\n" + under + content[end + 5 :]
    to_top = '<p class="to-top"><a href="#content">Back to the top</a></p>\n' if len(h2s) >= 6 else ""
    prev_page = PAGES[position - 1] if position > 0 else None
    next_page = PAGES[position + 1] if position + 1 < len(PAGES) else None
    pager = '<nav class="pager docs-pager" aria-label="Previous and next page">'
    if prev_page:
        pager += f'<a class="prev" href="{root}{href(prev_page[1])}"><span>Previous</span>{prev_page[2]}</a>'
    if next_page:
        pager += f'<a class="next" href="{root}{href(next_page[1])}"><span>Next</span>{next_page[2]}</a>'
    pager += "</nav>"
    edit_path = CHANGELOG[1].as_posix() if path == CHANGELOG[0] else f"docs/{path}"
    current = "tutorial" if path == "getting-started/tutorial.md" else "how" if path.startswith("how-it-works/") else "docs"
    crumb = (
        f'<div class="crumb-row"><p class="crumb"><a href="{root}{INDEX}">Docs</a> / '
        f'<a href="{root}{INDEX}#{section_id(section)}">{section}</a></p>{version_line(root)}</div>\n'
    )
    body = (
        header(root, current, href(path))
        + '<div class="wrap docs-layout">\n'
        + sidebar(root, path, contents)
        + f'<main class="prose" id="content">\n{crumb}'
        + content
        + f'\n{to_top}<p class="edit"><a href="{REPO}/edit/main/{edit_path}">{GH_SVG}Edit this page on GitHub</a></p>\n'
        + pager
        + "\n</main>\n"
        + "</div>\n"
        + footer(root)
    )
    url = href(path)
    entries = []
    for s in renderer.sections:
        kind = "page" if s["level"] == 1 else "changelog" if renderer.changelog else "section"
        entry = {
            "page": title,
            "section": s["heading"] if s["level"] > 1 else "",
            "url": url + (f"#{s['slug']}" if s["level"] > 1 else ""),
            "text": " ".join(s["text"]),
        }
        if s["code"]:
            # Each word of the section's code once: enough to find it by, and
            # never shown, as a snippet comes from the text.
            entry["code"] = " ".join(dict.fromkeys(w for w in " ".join(s["code"]).split() if re.search(r"\w", w)))
        entry["kind"] = kind
        entries.append(entry)
    for c in renderer.extra:
        entries.append({"page": title, "section": c["section"], "url": f"{url}#{c['slug']}", "text": c["text"], "kind": "section"})
    if path == PARAMETERS:
        entries += option_entries()
    html_page = page(root, f"{title} | AMAP docs", description, body, figures=renderer.figures, path=url)
    return DOCS / url, html_page, entries


def index_page() -> tuple[Path, str]:
    tasks = [
        ("deployment/hpc.html", "Run on a SLURM cluster"),
        ("user-guide/custom-reference.html", "Use another species"),
        ("deployment/troubleshooting.html", "My run failed"),
        ("user-guide/outputs.html", "Which file do I open"),
        ("user-guide/parameters.html", "Every option"),
    ]
    cards = [
        '<nav class="hub-tasks" aria-labelledby="common-tasks"><p class="hub-label" id="common-tasks">Common tasks</p><ul>'
        + "".join(f'<li><a href="{link}">{label}</a></li>' for link, label in tasks)
        + "</ul></nav>\n"
    ]
    section = None
    for sec, path, title, line in PAGES:
        if sec != section:
            if section is not None:
                cards.append("</div>\n")
            cards.append(f'<h2 class="docs-group" id="{section_id(sec)}">{sec}</h2>\n<div class="deep">\n')
            section = sec
        cards.append(f'<a href="{href(path)}"><b>{title}</b><span>{line}</span></a>\n')
    cards.append("</div>\n")
    body = (
        header("", "docs", INDEX)
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


# ── The page that is not there ──
#
# GitHub Pages serves 404.html for any address it has no file for, at any
# depth, so every link here starts from the site's own address. A script sends
# an address that names a page in another way (a trailing slash, no .html, a
# folder) to it, and offers the nearest page for anything else.

SITE_PATH = "/AMAP/"


def not_found_page() -> tuple[Path, str]:
    root = SITE_PATH
    known = [{"href": href(path), "title": title} for _sec, path, title, _ in PAGES]
    known += [{"href": "index.html", "title": "Front page"}, {"href": INDEX, "title": "Documentation"}]
    sections = sorted({section_id(sec) for sec, *_ in PAGES})
    script = (
        "<script>\n(function () {\n"
        f"  var ROOT = {json.dumps(root)};\n"
        f"  var PAGES = {json.dumps(known)};\n"
        f"  var SECTIONS = {json.dumps(sections)};\n"
        "  var path = decodeURIComponent(location.pathname);\n"
        "  var rel = path.indexOf(ROOT) === 0 ? path.slice(ROOT.length) : path.replace(/^\\/+/, \"\");\n"
        "  var stem = rel.replace(/\\/+$/, \"\").replace(/\\/index\\.html?$/, \"\").replace(/\\.html?$/, \"\").toLowerCase();\n"
        "  var base = function (h) { return h.replace(/\\.html$/, \"\"); };\n"
        "  var last = function (s) { return s.split(\"/\").pop(); };\n"
        "  // A page asked for without .html, or as a folder: go there.\n"
        "  if (!/\\.html?$/.test(rel) || /\\/$/.test(rel)) {\n"
        "    var exact = PAGES.filter(function (p) { return base(p.href) === stem; });\n"
        "    if (!exact.length && stem.indexOf(\"/\") < 0) exact = PAGES.filter(function (p) { return last(base(p.href)) === stem; });\n"
        "    if (exact.length === 1) { location.replace(ROOT + exact[0].href + location.search + location.hash); return; }\n"
        "    // A section of the docs: its place on the Docs tab.\n"
        "    if (SECTIONS.indexOf(stem) >= 0) { location.replace(ROOT + \"documentation.html#\" + stem); return; }\n"
        "  }\n"
        "  // Anything else: the page whose name is nearest.\n"
        "  var word = last(stem) || stem;\n"
        "  var distance = function (a, b) {\n"
        "    var row = [], i, j, prev, next;\n"
        "    for (j = 0; j <= b.length; j++) row.push(j);\n"
        "    for (i = 1; i <= a.length; i++) {\n"
        "      prev = row[0]; row[0] = i;\n"
        "      for (j = 1; j <= b.length; j++) {\n"
        "        next = row[j];\n"
        "        row[j] = Math.min(row[j] + 1, row[j - 1] + 1, prev + (a[i - 1] === b[j - 1] ? 0 : 1));\n"
        "        prev = next;\n"
        "      }\n"
        "    }\n"
        "    return row[b.length];\n"
        "  };\n"
        "  var best = null;\n"
        "  if (word) PAGES.forEach(function (p) {\n"
        "    var names = [last(base(p.href)), p.title.toLowerCase().replace(/[^a-z0-9]+/g, \"-\")];\n"
        "    names.forEach(function (n) {\n"
        "      var d = distance(word, n) / Math.max(word.length, n.length);\n"
        "      if (!best || d < best.d) best = { d: d, page: p };\n"
        "    });\n"
        "  });\n"
        "  var near = document.querySelector(\"[data-nearest]\");\n"
        "  if (near && best && best.d <= 0.5) {\n"
        "    var link = near.querySelector(\"a\");\n"
        "    link.href = ROOT + best.page.href;\n"
        "    link.textContent = best.page.title;\n"
        "    near.hidden = false;\n"
        "  }\n"
        "  // The search box starts with the words of the address, for the reader to\n"
        "  // correct; its results open when they do.\n"
        "  var box = document.querySelector(\".lost-search [data-docs-search]\");\n"
        "  if (box && word) box.value = word.replace(/[-_]+/g, \" \");\n"
        "})();\n</script>\n"
    )
    body = (
        header(root, "")
        + '<main id="content">\n  <section class="ref-hero lost-hero">\n    <div class="wrap">\n'
        "      <h1>That page is not here</h1>\n"
        "      <p>The link that brought you here points at something that has moved or never existed. Nothing is broken on your side.</p>\n"
        '      <p class="nearest" data-nearest hidden>The nearest page: <a href="' + root + 'documentation.html">Documentation</a></p>\n'
        '      <div class="docs-search lost-search">\n'
        '        <input type="search" placeholder="Search the docs" aria-label="Search the documentation" data-docs-search>\n'
        '        <div class="docs-results" data-docs-results hidden></div>\n      </div>\n'
        '      <figure class="lost" aria-label="A drawing in the site\'s style: the link you followed, drawn as a call, and no page of the site aligned under it.">\n'
        '        <div class="lost-scroll">\n'
        '        <svg class="trk" viewBox="0 0 560 96" role="img" aria-hidden="true">\n'
        '          <text class="row-label" x="0" y="28">Your link</text>\n'
        '          <line class="row-line" x1="96" x2="556" y1="24" y2="24"/>\n'
        '          <path class="g g-bakta" d="M190,14 H408 L420,24 L408,34 H190 Z"/>\n'
        '          <text class="gname" x="305" y="28" text-anchor="middle">this page</text>\n'
        '          <text class="row-label" x="0" y="74">The site</text>\n'
        '          <line class="row-line" x1="96" x2="556" y1="70" y2="70"/>\n'
        '          <text class="gsub" x="305" y="64" text-anchor="middle">no page aligns here</text>\n'
        "        </svg>\n        </div>\n"
        "        <figcaption>Case 1, for a page: nothing aligns where this link points.</figcaption>\n"
        "      </figure>\n"
        '      <div class="actions">\n'
        f'        <a class="btn btn-primary" href="{root}index.html">The front page</a>\n'
        f'        <a class="btn btn-ghost" href="{root}getting-started/tutorial.html">The tutorial</a>\n'
        f'        <a class="btn btn-ghost" href="{root}{INDEX}">All the docs</a>\n'
        "      </div>\n    </div>\n  </section>\n</main>\n"
        + footer(root)
    )
    note = (
        "<!-- Generated by tools/build_docs.py. GitHub Pages serves this for any missing page, at any\n"
        "     depth, so every link here starts from the site's own address. -->\n"
    )
    text = page(root, "Not here | AMAP", "The page you followed a link to has moved or never existed.", body,
                path="404.html", body_class="lost-page", note=note, tail=script)
    return DOCS / "404.html", text


# ── Search ──
#
# Besides the sections of every page: each option of the parameters page, the
# sections and questions of the front page, and the stories of how-it-works.html.


def option_entries() -> list[dict]:
    out = []
    for row in parameter_rows():
        cells = dict(row["cells"])
        cells.pop("Parameter", None)
        what = cells.pop("Description", "")
        rest = ". ".join(f"{k}: {v}" for k, v in cells.items() if v)
        out.append({
            "page": "Parameters",
            "section": row["flag"],
            "url": f"user-guide/parameters.html#opt-{row['flag'][2:]}",
            "text": f"{what} {rest}. In {row['section']}.".replace("..", "."),
            "kind": "option",
        })
    return out


def tag_text(fragment: str) -> str:
    """The words of an HTML fragment."""
    text = re.sub(r"<(script|style)\b.*?</\1>", " ", fragment, flags=re.S)
    return " ".join(html.unescape(re.sub(r"<[^>]+>", " ", text)).split())


FRONT_SECTIONS = ["start", "cases", "proof", "pipeline", "outputs", "report", "faq"]


def front_page_entries() -> list[dict]:
    front = (DOCS / "index.html").read_text(encoding="utf-8")
    description = re.search(r'<meta name="description" content="([^"]*)"', front)
    out = [{"page": "Front page", "section": "", "url": "index.html", "text": html.unescape(description.group(1)) if description else "", "kind": "page"}]
    for key in FRONT_SECTIONS:
        m = re.search(rf'<section\b[^>]*\bid="{key}"[^>]*>(.*?)</section>', front, re.S)
        if not m:
            raise SystemExit(f"tools/build_docs.py: docs/index.html has no section #{key} for the search index")
        block = m.group(1)
        title = re.search(r"<h2[^>]*>(.*?)</h2>", block, re.S)
        kicker = re.search(r'<p class="kicker">(.*?)</p>', block, re.S)
        name = tag_text(title.group(1)) if title else key
        lead = re.search(r'<div class="section-head">(.*?)</div>', block, re.S)
        text = tag_text(lead.group(1)) if lead else ""
        if kicker:
            text = text.replace(tag_text(kicker.group(1)), "", 1).strip()
        out.append({"page": "Front page", "section": name, "url": f"index.html#{key}", "text": text, "kind": "section"})
        if key == "faq":
            for q in re.finditer(r"<details\b[^>]*>\s*<summary>(.*?)</summary>(.*?)</details>", block, re.S):
                out.append({"page": "Front page", "section": tag_text(q.group(1)), "url": "index.html#faq",
                            "text": tag_text(q.group(2)), "kind": "faq"})
    return out


def story_entries() -> list[dict]:
    """Each story of how-it-works.html: its title and the words of its steps."""
    how = (DOCS / "assets" / "how.js").read_text(encoding="utf-8")
    nothing = dict(re.findall(r'^\s+(\w+): ("(?:\\.|[^"\\])*"),?$', how[how.index("var NOTHING"):how.index("var STORIES")], re.M))
    page_title = next(title for _sec, path, title, _ in PAGES if path == "how-it-works.html")
    out = [{"page": page_title, "section": "", "url": "how-it-works.html",
            "text": next(line for _sec, path, _t, line in PAGES if path == "how-it-works.html"), "kind": "page"}]
    stories = how[how.index("var STORIES"):]
    starts = list(re.finditer(r'^\s+key: "([\w-]+)", group: "([^"]*)", label: "([^"]*)", tag: "([^"]*)"', stories, re.M))
    if not starts:
        raise SystemExit("tools/build_docs.py: no stories found in docs/assets/how.js for the search index")
    for k, m in enumerate(starts):
        chunk = stories[m.end(): starts[k + 1].start() if k + 1 < len(starts) else len(stories)]
        steps = re.search(r"\btext: \[(.*?)\n\s+\]", chunk, re.S)
        words = []
        if steps:
            for item in re.finditer(r'"(?:\\.|[^"\\])*"|NOTHING\[(\w+)\]|NOTHING\.(\w+)', steps.group(1)):
                ref = item.group(1) or item.group(2)
                words.append(json.loads(nothing[ref]) if ref else json.loads(item.group(0)))
        key, group, label, tag = m.groups()
        out.append({"page": page_title, "section": label, "url": f"how-it-works.html#{key}",
                    "text": f"{tag}. {group}. " + " ".join(words), "kind": "story"})
    return out


# The index goes out packed: each page's title and address once, each kind
# once, and each entry as a row that points at them. The file unpacks itself
# into the entries the search box reads, {page, section, url, text, code?,
# kind}, as it loads; that keeps the whole of the docs under SEARCH_BUDGET.

SEARCH_BUDGET = 230_000

SEARCH_UNPACK = (
    "window.AMAP_DOCS_INDEX = (function (d) {\n"
    "  return d.rows.map(function (r) {\n"
    "    var p = d.pages[r[0]];\n"
    '    var e = { page: p[0], section: r[1], url: p[1] + (r[2] ? "#" + r[2] : ""), text: r[3] };\n'
    "    if (r[5]) e.code = r[5];\n"
    "    e.kind = d.kinds[r[4]];\n"
    "    return e;\n"
    "  });\n"
    "})("
)


def pack_index(entries: list[dict]) -> dict:
    pages: list[list[str]] = []
    kinds: list[str] = []
    rows = []
    for e in entries:
        base, _, fragment = e["url"].partition("#")
        key = [e["page"], base]
        if key not in pages:
            pages.append(key)
        if e["kind"] not in kinds:
            kinds.append(e["kind"])
        row = [pages.index(key), e["section"], fragment, e["text"], kinds.index(e["kind"])]
        if e.get("code"):
            row.append(e["code"])
        rows.append(row)
    return {"pages": pages, "kinds": kinds, "rows": rows}


def unpack_index(packed: dict) -> list[dict]:
    """What SEARCH_UNPACK makes of the packed index, for check_site.py."""
    out = []
    for r in packed["rows"]:
        page, base = packed["pages"][r[0]]
        e = {"page": page, "section": r[1], "url": base + ("#" + r[2] if r[2] else ""), "text": r[3]}
        if len(r) > 5 and r[5]:
            e["code"] = r[5]
        e["kind"] = packed["kinds"][r[4]]
        out.append(e)
    return out


def read_index(text: str) -> list[dict]:
    """The entries of a search-index.js as build_docs.py writes it."""
    return unpack_index(json.loads(text[text.index(SEARCH_UNPACK) + len(SEARCH_UNPACK):].rstrip().removesuffix(");")))


def search_index_js(entries: list[dict]) -> str:
    packed = pack_index(entries)
    if unpack_index(packed) != entries:
        raise SystemExit("tools/build_docs.py: the packed search index does not unpack to its entries")
    return (
        "/* Generated by tools/build_docs.py: the documentation, section by section, for the search box. */\n"
        + SEARCH_UNPACK
        + json.dumps(packed, ensure_ascii=False, separators=(",", ":"))
        + ");\n"
    )


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
    lost_path, lost_text = not_found_page()
    files[lost_path] = lost_text
    entries += front_page_entries() + story_entries()
    files[SEARCH_INDEX] = search_index_js(entries)
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
