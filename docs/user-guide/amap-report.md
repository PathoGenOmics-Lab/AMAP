# The report by question: amap_report.html

Every run writes `results/amap_report.html`, next to
[the full report](report.md), `final_report.html`, when the run has one: one
file that answers, section by section, the questions a run raises. Did each
genome annotate well? Where on the reference are genes missing or broken? What
happened to one gene in every genome? What did the correction change? Every
table filters, sorts and saves as CSV, every figure saves as SVG, and a click on
a genome or a gene opens it.

It is written with or without the functional annotation. It needs only what the
correction leaves, so a run with `--skip_functional_annotation` has no
`final_report.html`; `amap_report.html` is still written, without its eggNOG
sections, which say that the run skipped them.

It is self-contained: the data travels inside the file, compressed, and nothing
is fetched from anywhere, so it opens with no network and survives being
emailed. Fifteen genomes of the *Mycobacterium tuberculosis* complex (MTBC) make a file of about half a megabyte, against
4.5 MB for `final_report.html` with twelve. It is dark, like the website, with a
switch to light at the top right.

## The sections

### 1. Did each genome annotate well?

One row per genome: lineage, genome size, contigs and N50, Bakta's calls and
AMAP's genes, the share of the reference genes it has, and how many are missing,
pseudogenes, partial, found more than once or outside the reference. A value far
from those of comparable genomes is flagged (see [How a value is
flagged](#how-a-value-is-flagged)). A click on a genome opens it: its numbers,
and its genes in each state, each list searchable and saved as CSV.

### 2. Where on the reference are genes missing or broken?

Every reference gene at its place on the reference, one row per genome, the
genomes grouped by lineage. A missing gene fills its row; a pseudogene is a band
through the middle; a second copy a band on top. Deletions that genomes share
show as blocks: with MTBC genomes, the regions of difference. Drag across the
map to zoom into a region; its genes are then listed with how many genomes have
each in each state. Below the map, **Blocks of missing genes** lists every run
of two or more consecutive reference genes missing from a genome, grouped where
genomes overlap, and says which lineages lack each block in all their genomes.

### 3. What happened to one gene in every genome?

Type a locus, a gene name or part of a product. The card gives the reference's
own coordinates, length, product, functional category and, from Mycobrowser,
whether the reference calls it a pseudogene and which drugs it has resistance
mutations for. The table gives, for every genome, its state, locus tag, contig,
coordinates, strand, length against the reference, the note of a pseudogene and
how AMAP got the gene.

### 4. What did AMAP change in Bakta's annotation?

A waterfall from Bakta's calls to AMAP's genes: the case 1 calls the filter
removed, the case 4 calls the reference replaced, then the genes placed from the
reference over a call and those added where Bakta called nothing. For one genome
a click on a bar lists its genes; for the calls AMAP did not keep, Bakta's
coordinates next to those of the gene that took their place. [The correction
algorithm](../how-it-works/correction-algorithm.md) says what each case means.

### 5. Which genes are broken, and in which lineages?

Every gene AMAP flags as a pseudogene in at least one genome, against the
genomes, with the reason. First the genes broken in every genome of some
lineages and in no other, then those broken in every genome, then those broken
in all but one or two, then the rest. With Mycobrowser, each is compared with
the reference's own pseudogenes.

### 6. What do the genomes carry that the reference does not?

How many genomes carry each reference gene, with the genes some genomes lack;
the genes that name no reference gene, grouped by product across genomes; the
reference genes with more than one copy; and, when the run has them, the
pangenome and the homologous communities.

### 7. Do the lineages differ in function?

The genes each lineage loses (missing or pseudogene) by Mycobrowser category,
and with the functional annotation, genes by COG category and KEGG pathway. One
line per category, a dot per lineage, the categories that differ most between
lineages first.

### 8. How was this made?

AMAP's version, the command, the reference, every parameter
(`pipeline_info/params.yml`), every tool's version
(`pipeline_info/software_versions.tsv`), and what each state means.

The bar under the top one filters every section at once: a lineage's pill shows
or hides it (Alt or Cmd and a click keeps that lineage alone), and the box keeps
the genomes whose name contains what is typed.

## Gene states

Each reference gene gets one state in each genome, from the genes that name it
in the genome's final annotation, `AMAP/{lineage}_{sampleId}/{lineage}_{sampleId}.AMAP.tsv`.

| State | What it means |
| :--- | :--- |
| Intact | A gene names it, and nothing below applies |
| Missing | No gene names it |
| Pseudogene | The best gene naming it is flagged as a pseudogene: an internal stop, a lost stop or a frame that runs on past the reference protein |
| Partial | The best gene naming it runs into a contig end and is cut there |
| More copies | More than one gene names it alone |
| Length differs | One intact copy, of another length than the reference gene |
| Shared copy | Only a gene that names several reference genes at once names it, such as an IS6110 copy that matches all seventeen; counted as found |

When several apply, the first in the table wins: a pseudogene that is also cut
by a contig end is a pseudogene. "Reference genes found" counts every state but
missing.

## How a value is flagged

Each genome is compared with its own lineage when that lineage has five genomes
or more, and with every other genome otherwise. A value is flagged when it lies
more than 3.5 robust standard deviations (1.4826 times the median absolute
deviation) from their median, on the side that is bad: more missing genes,
pseudogenes, partial genes, contigs or genes outside the reference, fewer
reference genes found, a smaller N50, a genome size either way. Each also has to
be further from the median than a minimum (five genes, five contigs, half a
point of reference genes found, 1% of the genome size or 5% of the N50), so that
identical genomes do not flag a gene or two of difference. Fewer than four
genomes are not compared at all.

A flag says a genome is unusual, not that it is wrong. *M. bovis* among
*M. tuberculosis* genomes is flagged for its missing genes, which are its
regions of difference. A genome with twice the pseudogenes of the others is
worth a look at its assembly: frameshifts from sequencing errors look like
pseudogenes.

## What each section needs

| Input | From | Without it |
| :--- | :--- | :--- |
| `*.AMAP.tsv` | `FIX_CDS` | No report |
| The reference proteome | the reference | No report |
| Case summaries and tables | `TYPE_CASES` | No lineages, no correction section |
| Bakta's folder | `BAKTA_ANNOTATOR` | No genome sizes, contigs or products of Bakta's calls |
| Mycobrowser table | the reference, MTBC only | Names and products from the genomes' own annotation; no categories or drug resistance |
| Master tables | `GET_OUTPUTS`, functional annotation | No COG categories or KEGG pathways |
| Pangenome matrix | `REPORT_INPUTS`, functional annotation | No pangenome |
| Homologous communities | `GET_ALLDUP`, functional annotation | No communities |
| `params.yml`, `software_versions.tsv` | `DUMP_PARAMS`, `SOFTWARE_VERSIONS` | No parameters or versions |

A section without its inputs says which and why instead of drawing anything.
With a custom reference, under `--skip_typing` and `--skip_mycobrowser`, every
genome is `Untyped` and the names and products come from Bakta and eggNOG.

A genome that was classified but has no corrected genes, because the pipeline
set aside its correction after it kept failing, is named at the top of the
report. One whose Bakta, miniprot or classification task was set aside never
reaches the report: see [a sample is missing](../deployment/troubleshooting.md#a-sample-is-missing-from-the-results-ignored--1).

## Size

The data of each genome is a few kilobytes: its genes are stored as how far each
lies from where the previous gene and the reference predict it, which for a
genome collinear with its reference is mostly nothing. Up to 100 genomes, the
report also lists every Bakta call AMAP did not keep, with its coordinates; past
that it only counts them, to keep the file small, and `CDS_CASES/{lineage}_{sampleId}/` has
them for any run.

| Genomes | Report |
| :--- | :--- |
| 15 | 0.5 MB, measured |
| 100 | 1.5 to 1.7 MB |
| 300 | 1.9 to 2.1 MB |
| 500 | 2.9 to 3.3 MB |

The sizes past fifteen are worked out from the fifteen MTBC genomes of eight
lineages the report was measured on, the larger figure with each genome's data
compressed on its own and the smaller with the sharing measured between those
fifteen. In the browser, five hundred genomes open in about three seconds.

## Writing one by hand

`bin/report/amap_report.py` needs only Python 3.9 or later. Pointed at a
`results/` folder it reads the files where the pipeline publishes them:

```bash
python3 bin/report/amap_report.py \
    --amap results/AMAP \
    --cases results/CDS_CASES \
    --bakta results/BAKTA \
    --master results/AMAP \
    --ref-proteome ref_proteome.faa \
    --mycobrowser Mycobacterium_tuberculosis_H37Rv_txt_v4.txt \
    --versions results/pipeline_info/software_versions.tsv \
    --params results/pipeline_info/params.yml \
    --out amap_report.html
```

The reference proteome is the one the run used: for H37Rv, the AMAP container
has it at `/reference/aa_H37Rv.translated.fas`, and the Mycobrowser table at
`/reference/Mycobacterium_tuberculosis_H37Rv_txt_v4.txt`.

> [!NOTE]
> The page reads its data with the browser's own decompression, which every
> current browser has: Chrome and Edge 80, Firefox 113 and Safari 16.4, and
> later. An older one says so instead of drawing the report.
