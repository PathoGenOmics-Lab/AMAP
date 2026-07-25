---
title: Outputs
---

# Outputs

Everything lands under `--outdir`, `results/` by default, in one subdirectory
per pipeline stage.

## Where to look first

| What you want | File |
| :--- | :--- |
| The whole analysis, interactive | `results/final_report.html` |
| Corrected annotation for a sample | `results/AMAP/{sampleId}/*-AMAP.gff3` |
| Proteins for downstream tools | `results/AMAP/{sampleId}/*-AMAP.clean.faa` |
| Everything known about every gene | `results/AMAP/{sampleId}/*-AMAPmastertable.tsv` |
| Tool versions for your methods section | `results/pipeline_info/software_versions.yml` |

## The full tree

```text
results/
├── AMAP/
│   └── {sampleId}/
│       ├── *-AMAP.faa                  # Corrected proteins, all genes including pseudogenes
│       ├── *-AMAP.clean.faa            # No pseudogenes, no internal stops
│       ├── *-AMAP.fna                  # Corrected nucleotide sequences
│       ├── *-AMAP.gff3                 # Corrected annotation, GFF3
│       ├── *-AMAP.gbk                  # Corrected annotation, GenBank
│       ├── *-AMAP.qc.tsv               # Genome size, GC%, coding density, annotation rate
│       ├── *-AMAPmastertable.tsv       # Coordinates plus eggNOG plus Mycobrowser, one row per gene
│       ├── *-AMAP.COG.tsv              # COG category counts
│       └── *-AMAP.MycoBrowser.tsv      # Mycobrowser category counts
│
├── BAKTA/
│   └── {sampleId}/                     # The uncorrected Bakta annotation, for comparison
│       ├── *.faa  *.gbff  *.gff3  *.tsv
│       ├── *.hypotheticals.faa  *.hypotheticals.tsv
│       ├── *.png                       # Genome map
│       └── *.txt                       # Bakta summary statistics
│
├── CDS_CASES/
│   ├── summary/
│   │   └── *-summary.tsv               # Case counts, all samples
│   └── {sampleId}/
│       ├── *-case[1-4].tsv             # One table per case
│       ├── *-rvmissing.tsv             # Reference genes Bakta did not find
│       ├── *.html                      # Interactive overlap plots
│       └── *_corrected.tsv             # Final corrected coordinates
│
├── DUPLICATIONS/
│   ├── {sampleId}/
│   │   ├── duplications*               # Raw Diamond output, within-genome
│   │   ├── duplications_*_groups.tsv   # Duplicated genes grouped into families
│   │   └── *.dup                       # Duplicated protein sequences
│   ├── homologous_communities.tsv      # Communities across all samples, with network metrics
│   └── homologous_communities.html     # Interactive network
│
├── PANGENE/
│   └── {sampleId}/
│       ├── *.gfa                       # Pangenome graph
│       ├── *.bubbles.txt               # Structural variation in the graph
│       ├── *.Rtab                      # Gene presence/absence matrix
│       └── *.stat                      # Graph statistics
│
├── PATHOTYPR/
│   ├── typing/*_summary.tsv            # Lineage assignment
│   ├── fastas/*.fasta                  # Genomes renamed by lineage
│   └── renamed_fastas.tsv              # Original ID to lineage-based name
│
├── pipeline_info/
│   ├── execution_report_*.html         # Nextflow execution report
│   ├── execution_timeline_*.html       # Timeline
│   ├── pipeline_dag_*.html             # Dependency graph
│   ├── execution_trace_*.txt           # Per-process resource usage
│   ├── software_versions.yml / .tsv    # Exact tool versions
│   ├── params.yml                      # Parameters used for this run
│   ├── validation.log                  # Pre-flight validation results
│   └── run_summary.txt                 # Status, duration, sample counts
│
└── final_report.html                   # The full interactive report
```

## The key files

### `*-AMAP.gff3` and `*-AMAP.gbk`

The corrected annotation, which is the point of the pipeline. The GFF3 is
spec-compliant: gene features parent their CDS features, attributes are
URL-escaped, `##sequence-region` pragmas are emitted per contig, and eggNOG
annotation rides along as `Ontology_term`, `eC_number`, `COG_category`,
`KEGG_ko`, `KEGG_Pathway` and `PFAMs`. The sequences are appended after
`##FASTA`.

### `*-AMAP.faa` versus `*-AMAP.clean.faa`

`*-AMAP.faa` mirrors the annotation exactly, pseudogenes included, so gene
counts match the GFF3.

`*-AMAP.clean.faa` drops pseudogenes and any sequence with an internal stop
codon, and strips the terminal stop. This is the file for orthology,
phylogenomics, pangenome tools and dN/dS, all of which expect valid
translations. It has fewer sequences than the GFF3 has genes, by design.

### `*-AMAPmastertable.tsv`

One row per gene, joining the corrected coordinates, the Bakta annotation, the
eggNOG assignment (COG category, description, preferred name, GO, EC, KEGG,
PFAMs, CAZy) and, for MTBC, the Mycobrowser functional category. This is the
table to load when you want to ask a question about the annotation.

### `*-AMAP.qc.tsv`

One row per sample: genome size, contig count, N50, GC%, total CDS, pseudogene
count, hypothetical protein percentage, mean and median gene length, coding
density, and the eggNOG annotation rate. Useful for spotting the sample that
went wrong.

### `CDS_CASES/`

The audit trail for the correction. Each case table lists the genes that fell
into it, so you can see exactly what was changed and why. The
[correction algorithm](../how-it-works/correction-algorithm.md) explains the
cases; `*-summary.tsv` aggregates the counts across samples and feeds the
report.

### `final_report.html`

Self-contained, so it travels by email. Covers per-sample quality, missing
reference genes by lineage, PE/PPE family distribution, KEGG pathways,
duplication communities and network topology, the pangenome (core, softcore,
shell, cloud, frequency spectrum, presence/absence heatmaps), COG and
Mycobrowser distributions, case statistics and the associated tests, plus a
searchable master-table explorer.

## Effect of the skip options

`--skip_typing`
: No `PATHOTYPR/` directory. Samples keep their original IDs instead of being
  renamed by lineage.

`--skip_functional_annotation`
: The run stops after correction. No `AMAP/`, `DUPLICATIONS/`, `PANGENE/` or
  `final_report.html`. You still get `BAKTA/` and `CDS_CASES/`, including
  `*_corrected.tsv` with the corrected coordinates.

`--skip_mycobrowser`
: `*-AMAP.MycoBrowser.tsv` is written empty and the Mycobrowser sections drop
  out of the report. Everything else is unaffected.

## Next steps

- [How it works](../how-it-works/index.md): what produced all this
- [Duplication metrics](../how-it-works/duplications.md): reading the network numbers
