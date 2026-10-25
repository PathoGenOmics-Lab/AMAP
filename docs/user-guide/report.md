# The full report: final_report.html

A run with the functional annotation ends in `results/final_report.html`, the
whole analysis in one file, built after eggNOG, the duplications and the
pangenome. Every run, with or without them, also writes
[the report by question](amap-report.md), `results/amap_report.html`. So a run
with `--skip_functional_annotation` has no `final_report.html`;
`amap_report.html` is still written, without its eggNOG sections.

It is self-contained: all the figures, tables and JavaScript are embedded, so
it opens with no network and survives being emailed to someone who has never
heard of Nextflow.

To see what the sections look like before running anything, open the
[example report](../examples/demo_report.html): a rendered report for twelve
genomes of the *Mycobacterium tuberculosis* complex (MTBC) across four lineages.

> [!NOTE]
> **The example uses generated data**
>
> The twelve samples in it are synthetic. The distributions are drawn around
> a real *M. tuberculosis* CDC1551 run, so the shapes are representative,
> but no conclusion should be drawn from the numbers themselves.
>
> `bin/report/demo/make_demo.py` builds that dataset if you want to render
> the report yourself without waiting for a pipeline run.

## What is in it

<!-- figure: report-map -->

### Overview

Dataset composition, then a per-sample quality panel: CDS counts, correction
rates and missing reference genes. The annotation quality tab carries coding
density, hypothetical percentage and genome size, which is where an outlier
sample declares itself before you have gone any further.

### Annotation results

CDS distribution per lineage, the genes with no counterpart in the reference,
and the correction summary: how many genes fell into each case and what the
adjudication did with them. Read this alongside
[the correction algorithm](../how-it-works/correction-algorithm.md).

Missing reference genes get both a heatmap, for spotting a gene absent from a
whole lineage rather than one sample, and a searchable table.

### Functional annotation

COG categories from eggNOG, Mycobrowser functional categories for MTBC, and a
comparison against the reference distribution. Each has its own per-category
breakdown in counts and percentages, with the statistical tests alongside.

PE/PPE gene families get their own section, since they are both the hardest
part of an MTBC genome to annotate and usually the most interesting. The KEGG
pathways follow: the most represented metabolic pathways, and how they differ
between lineages.

The pseudogenes close it: how many AMAP called in each sample, lineage by
lineage, a sample with none at 0, the same genes the QC table counts and the
GFF3 marks `pseudo=true`. The master table's `note` column says why each is a
pseudogene.

### Duplications

Community sizes, cross-genome homologous pairs and the network topology
metrics. [Duplication metrics](../how-it-works/duplications.md) explains what
each number means.

### Pangenome

Core, softcore, shell and cloud classification, the gene frequency spectrum,
and presence/absence heatmaps.

### CDS statistics

How many genes fall in each case, lineage by lineage: box and violin plots with
Dunn's tests between lineages, the share of each case in every lineage, a
heatmap of the means, and a table per case with the count, mean, median and
spread of each lineage.

### Master table explorer

Every gene of every sample in one searchable, sortable table, with the
functional annotation attached. Useful for answering a specific question
without leaving the report.

## What changes when steps are skipped

- **`--skip_typing`**: Everything still renders, but every sample is grouped under `Untyped`, so
  the per-lineage comparisons have a single group and the statistical tests
  drop out.
- **`--skip_mycobrowser`**: The Mycobrowser and PE/PPE sections are omitted. COG, KEGG, pangenome and
  duplications are unaffected. This is automatic for non-MTBC organisms.
- **`--skip_functional_annotation`**: The run stops after correction: no `final_report.html`;
  `amap_report.html` is still written, without its eggNOG sections.

## Rendering it yourself

The report is `bin/report/report_AMAP.Rmd`, rendered by the `FINAL_REPORT`
process. To iterate on it without running the pipeline:

```bash
python3 bin/report/demo/make_demo.py demo_inputs

docker run --rm -v "$PWD:/w" -w /w paururo/amap:latest \
  Rscript -e 'rmarkdown::render("bin/report/report_AMAP.Rmd",
    output_file="/w/demo.html",
    params=list(
      version="v1.0.1",
      cog_file="/w/demo_inputs/input_cog.txt",
      myco_file="/w/demo_inputs/input_myco.txt",
      summary_file="/w/demo_inputs/input_summary.txt",
      empty_file="/w/demo_inputs/input_noref.txt",
      rvmissing_file="/w/demo_inputs/input_rvmissing.txt",
      mastertable_file="/w/demo_inputs/input_mastertables.txt",
      pangene_file="/w/demo_inputs/input_pangene.txt",
      dup_file="/w/demo_inputs/homologous_communities.tsv",
      qc_file="/w/demo_inputs/input_qc.txt"))'
```

> [!WARNING]
> **Keep it self-contained**
>
> Nothing in the Rmd, `header.html` or `styles.css` may point at a remote
> URL. Pandoc fetches and inlines a remote resource when the machine
> rendering has a network, which hides the problem, and leaves it dangling
> when it does not. A report rendered on an isolated compute node then
> arrives with a broken image and a browser trying to reach the internet.
>
> Both logos live next to the Rmd for this reason, and the image zoom is
> plain DOM rather than jQuery from a CDN.
