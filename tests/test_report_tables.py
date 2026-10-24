#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The report reads the tables REPORT_INPUTS hands it with read.table, whose
default comment character is "#": a line stops at the first one, the header
line included. pangene names every genome <name>#0, so the report read the
pangenome matrix one genome column short, took the gene names for row names,
and never drew the heatmap. The tables a "#" can reach are read with no comment
character: the pangenome matrix, and the ones that list the reference
annotation's gene IDs, which extract_ref_proteins.py keeps as they are.
PANGENE now names each column after its sample, without the "#0".

The last test runs the report's own R code on the matrix PANGENE writes, a
column per sample, and counts the genes they all share. It needs Rscript and
the packages of the pangenome section (the AMAP container has both) and is
skipped without them, as in CI.
"""

import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
REPORT = ROOT / "bin" / "report" / "report_AMAP.Rmd"
MERGE_RTABS = ROOT / "bin" / "merge_rtabs.py"

# Rows of what PANGENE wrote for three typed samples, the first 2 Mb of
# M. bovis, CDC1551 and H37Rv against H37Rv: genes all three share, RD4, which
# M. bovis lacks, the phiRv1 prophage, which CDC1551 lacks, and two genes no
# reference gene matches, under the sample whose copy pangene kept.
RTABS = {
    "pangenome.Rtab": (
        "Gene\tA4_Mbovis\tL4_CDC1551\tL4_H37Rv\n"
        "Rv0001\t1\t1\t1\n"
        "Rv0002\t1\t1\t1\n"
        "Rv0003\t1\t1\t1\n"
        "Rv1506c\t0\t1\t1\n"
        "Rv1573\t1\t0\t1\n"
        "A4_Mbovis_01833\t1\t1\t0\n"
        "L4_CDC1551_01396\t0\t1\t0\n"),
}
GENOMES = ["A4_Mbovis", "L4_CDC1551", "L4_H37Rv"]

# params the report reads with read.table that a "#" can reach
HASH_TABLES = ["pangene_file", "rvmissing_file", "mastertable_file", "dup_file"]
READS = re.findall(r"read\.table\(params\$(\w+),([^\n]*)", REPORT.read_text())

PANGENE_CHUNKS = ["pangene_core", "pangene_spectrum", "pangene_heatmap"]
R_PACKAGES = ["dplyr", "ggplot2", "ggpubr", "pheatmap"]


@pytest.fixture
def merged(tmp_path):
    """input_pangene.txt, written from the Rtabs as REPORT_INPUTS does."""
    for name, text in RTABS.items():
        (tmp_path / name).write_text(text)
    subprocess.run([sys.executable, str(MERGE_RTABS), "--output", "input_pangene.txt",
                    *sorted(RTABS)], cwd=tmp_path, check=True, capture_output=True)
    return tmp_path / "input_pangene.txt"


def test_the_matrix_reaches_the_report_a_column_per_sample(merged):
    assert merged.read_text() == RTABS["pangenome.Rtab"]


@pytest.mark.parametrize("param", HASH_TABLES)
def test_read_with_no_comment_character(param):
    calls = [args for name, args in READS if name == param]
    assert calls, f"the report no longer reads params${param} with read.table"
    for args in calls:
        assert 'comment.char=""' in args.replace(" ", ""), (
            f"params${param} is read with the default comment character, so a '#' "
            f"in it cuts the line, header included")


def _rscript():
    exe = shutil.which("Rscript")
    if not exe:
        pytest.skip("needs Rscript (the AMAP container has it)")
    wanted = ", ".join(f'"{p}"' for p in R_PACKAGES)
    check = subprocess.run(
        [exe, "-e", f"quit(status = !all(sapply(c({wanted}), requireNamespace, quietly = TRUE)))"],
        capture_output=True)
    if check.returncode:
        pytest.skip(f"needs the R packages {', '.join(R_PACKAGES)}")
    return exe


def _pangene_code():
    """The report's R that loads the matrix, and the pangenome section's chunks."""
    text = REPORT.read_text()
    load = re.search(r"^# Pangene Rtab\n.*?^\}\n", text, re.S | re.M).group(0)
    chunks = dict(re.findall(r"^```\{r (\w+)[^}]*\}\n(.*?)^```", text, re.S | re.M))
    return load, [chunks[label] for label in PANGENE_CHUNKS]


def test_the_report_uses_every_genome_column(merged):
    exe = _rscript()
    load, chunks = _pangene_code()
    libraries = "; ".join(f"library({p})" for p in R_PACKAGES)
    (merged.parent / "pangene.R").write_text("\n".join([
        f"suppressPackageStartupMessages({{{libraries}}})",
        'params <- list(pangene_file = "input_pangene.txt")',
        load,
        'cat("columns:", colnames(pangene_raw), sep = "\\t"); cat("\\n")',
        *chunks,
        'cat("samples:", n_samples, "core:", sum(gene_class == "Core (100%)"), "\\n")',
    ]))
    run = subprocess.run([exe, "pangene.R"], cwd=merged.parent, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    assert "columns:\t" + "\t".join(["Gene", *GENOMES]) in run.stdout
    assert "samples: 3 core: 3" in run.stdout
    assert "No pangene data available" not in run.stdout
    assert "Could not generate heatmap" not in run.stdout
