# Outputs

Everything lands under `--outdir`, `results/` by default, in one subdirectory
per pipeline stage. `{lineage}_{sampleId}` below is the sample's name as
lineage typing leaves it, `L4_G26689` for sample `G26689`, or the sample ID
alone under `--skip_typing`.

## Where to look first

| What you want | File |
| :--- | :--- |
| The whole analysis, after a run with eggNOG ([the full report](report.md)) | `results/final_report.html` |
| Each genome and each gene, question by question, in every run ([the report by question](amap-report.md)) | `results/amap_report.html` |
| Corrected annotation for a sample, after a run with eggNOG | `results/AMAP/{lineage}_{sampleId}/*-AMAP.gff3` |
| Proteins for downstream tools, after a run with eggNOG | `results/AMAP/{lineage}_{sampleId}/*-AMAP.clean.faa` |
| Everything known about every gene, after a run with eggNOG | `results/AMAP/{lineage}_{sampleId}/*-AMAPmastertable.tsv` |
| The corrected genes, their proteins and sequences, in every run | `results/AMAP/{lineage}_{sampleId}/*.AMAP.tsv`, `.faa`, `.fna` |
| Tool versions for your methods section | `results/pipeline_info/software_versions.yml` |

## The full tree

<!-- figure: results-map -->

```text
results/
├── AMAP/
│   └── {lineage}_{sampleId}/
│       ├── *.AMAP.tsv                  # Every run: the corrected genes, one row each
│       ├── *.AMAP.faa                  # Every run: their proteins
│       ├── *.AMAP.fna                  # Every run: their nucleotide sequences
│       ├── *-AMAP.faa                  # Corrected proteins, pseudogenes included
│       ├── *-AMAP.clean.faa            # No pseudogenes, no internal stops
│       ├── *-AMAP.fna                  # Corrected nucleotide sequences
│       ├── *-AMAP.tsv                  # The finished annotation, every feature
│       ├── *-AMAP.gff3                 # Corrected annotation, GFF3
│       ├── *-AMAP.gbk                  # Corrected annotation, GenBank
│       ├── *-AMAP.qc.tsv               # Size, GC%, coding density, annotation rate
│       ├── *-AMAPmastertable.tsv       # Every gene, with eggNOG and Mycobrowser
│       ├── *-AMAP.COG.tsv              # COG category counts
│       ├── *-AMAP.MycoBrowser.tsv      # Mycobrowser category counts
│       └── *-AMAP.noRef.tsv            # The genes that name no reference gene
│
├── BAKTA/
│   └── {lineage}_{sampleId}/           # Bakta's own annotation, uncorrected
│       ├── *.gff3  *.gbff  *.embl  *.tsv  *.json
│       ├── *.faa  *.ffn  *.fna
│       ├── *.hypotheticals.faa  *.hypotheticals.tsv
│       ├── *.inference.tsv
│       ├── *.png  *.svg                # Genome map
│       └── *.txt  *.log                # Bakta summary statistics and log
│
├── CDS_CASES/
│   ├── summary/
│   │   └── *-summary.tsv               # Case counts, all samples
│   └── {lineage}_{sampleId}/
│       ├── *-case[1-4].tsv             # One table per case
│       ├── *-rvmissing.tsv             # Reference genes Bakta missed, put back
│       ├── *.html                      # Interactive overlap plots
│       └── *_corrected.tsv             # As classified, before the stop codon check
│
├── DUPLICATIONS/
│   ├── {lineage}_{sampleId}/
│   │   ├── duplications_*              # Diamond hits within the genome, sorted
│   │   ├── duplications_*_groups.tsv   # Duplicated genes grouped into families
│   │   ├── duplications_*_network.png  # The families drawn
│   │   └── *_interactive.html          # The same, interactive
│   ├── homologous_communities.tsv      # Communities across samples, with metrics
│   └── homologous_communities.html     # Interactive network
│
├── PANGENE/
│   ├── pangenome.gfa                   # Pangenome graph of the whole run
│   ├── pangenome.bubbles.txt           # Structural variation in the graph
│   ├── pangenome.Rtab                  # Gene presence/absence, a column per sample
│   └── pangenome.stat                  # Graph statistics
│
├── PATHOTYPR/
│   ├── typing/*_summary.tsv            # Lineage calls, one row per contig with markers
│   └── fastas/*.fasta                  # Each genome as {lineage}_{sampleId}.fasta
│
├── pipeline_info/
│   ├── execution_report_*.html         # Nextflow execution report
│   ├── execution_timeline_*.html       # Timeline
│   ├── pipeline_dag_*.html             # Dependency graph
│   ├── execution_trace_*.txt           # Per-process resource usage
│   ├── software_versions.yml / .tsv    # Exact tool versions
│   ├── params.yml                      # Parameters used for this run
│   ├── validation.log                  # The checks' verdict, when they pass
│   └── run_summary.txt                 # Status, duration, task counts and the command
│
├── final_report.html                   # The full report, after eggNOG
└── amap_report.html                    # The report by question, in every run
```

The `.AMAP.` files, with a dot, are the correction's own output, written in
every run, `--skip_functional_annotation` included. The `-AMAP.` files, with a
hyphen, are the finished annotation, written after eggNOG. `*.AMAP.faa` and
`*-AMAP.faa` hold the same proteins, and the two `.fna` files the same
sequences.

## The key files

### `*-AMAP.gff3` and `*-AMAP.gbk`

The corrected annotation, which is the point of the pipeline. The GFF3 is
spec-compliant: gene features parent their CDS features, attributes are
URL-escaped, `##sequence-region` pragmas are emitted per contig, and eggNOG
annotation rides along as `Ontology_term`, `eC_number`, `COG_category`,
`KEGG_ko`, `KEGG_Pathway` and `PFAMs`. The sequences are appended after
`##FASTA`.

Two kinds of gene carry marks, written the way NCBI writes them:

- **Pseudogenes** have `pseudo=true` in the GFF3 and `/pseudo` in the GenBank,
  with the reason in `Note` and `/note`: internal stop codons, a frameshift or
  lost stop that makes the frame run on past the end of the reference protein,
  no stop at all, or a stop the reference itself lacks. A `/pseudo` CDS has no
  `/translation`.
- **Genes a contig end cut short**, in a draft assembly, have `partial=true`,
  with `start_range=.,START` or `end_range=END,.` on the open end, and `<` or
  `>` on that end of the GenBank location. They are not pseudogenes. Each
  piece of a cut gene is there under the gene's name, the first in the
  assembly with the gene's identifier and the others with `_2`, `_3`.

The same three facts, `pseudogene`, `partial` (`5'`, `3'` or both) and `note`,
are the last columns of the master table and of `*-AMAP.tsv`.

<!-- figure: excerpt -->

A gene across the origin of a circular sequence is written the way Bakta writes
it. Bakta takes a sequence for a circular replicon when its FASTA description
says `complete`, as NCBI's "complete genome" and "complete sequence" do, or
`circular=true`, as Unicycler's does. In `*-AMAP.tsv` and the master table such
a gene has its start after its end (stcE on plasmid pO157 of *E. coli* O157:H7
Sakai is 92,563 to 2,502); in the GenBank it is `join(92563..92721,1..2502)` in a
record marked circular; and in the GFF3 it ends past the sequence's last base,
by the sequence's length (92,563 to 95,223 on the 92,721 bp plasmid), on a
`region` feature marked `Is_circular=true`, as the GFF3 specification and
Bakta's own GFF3 have it.

### `*-AMAP.faa` versus `*-AMAP.clean.faa`

`*-AMAP.faa` mirrors the annotation exactly, pseudogenes included, so gene
counts match the GFF3.

`*-AMAP.clean.faa` drops the genes AMAP calls pseudogenes, those with
`pseudo=true` in the GFF3, and any sequence with an internal stop codon, and
strips the terminal stop. This is the file for orthology, phylogenomics,
pangenome tools and dN/dS, all of which expect valid translations. It has fewer
sequences than the GFF3 has genes, by design.

### `*-AMAPmastertable.tsv`

One row per gene, joining the corrected coordinates, the Bakta annotation, the
eggNOG assignment (COG category, description, preferred name, GO, EC, KEGG,
PFAMs, CAZy) and, for the *Mycobacterium tuberculosis* complex (MTBC), the Mycobrowser functional category. This is the
table to load when you want to ask a question about the annotation.

### `*-AMAP.qc.tsv`

One row per sample: genome size, contig count, N50, GC%, total CDS, how many of
them AMAP calls pseudogenes (those with `pseudo=true` in the GFF3), hypothetical
protein percentage, mean and median gene length, coding density, and the eggNOG
annotation rate. Useful for spotting the sample that went wrong.

### `CDS_CASES/`

The audit trail for the correction. Each case table lists the genes that fell
into it, so you can see exactly what was changed and why. The
[correction algorithm](../how-it-works/correction-algorithm.md) explains the
cases; `*-summary.tsv` aggregates the counts across samples and feeds the
report. Its last three columns say which sample the counts are for: `Sequence`
is the name its files carry (`L4_G26689`), `Lineage` what typing called it
(`Untyped` under `--skip_typing`) and `Sample` its ID in the sample sheet
(`G26689`).

#### What a healthy summary looks like

Two real runs, those of the [tutorial](../getting-started/tutorial.md), to
calibrate expectations. On the left, *M. tuberculosis* CDC1551 against the
bundled H37Rv reference. On the right, *E. coli* O157:H7 Sakai against K-12
MG1655, using a custom reference.

| Case | CDC1551 against H37Rv | Sakai against K-12 |
| :--- | ---: | ---: |
| case2, perfect match | 2,813 | 3,324 |
| case3, imperfect match | 24 | 31 |
| case4, conflicting overlap, resolved | 1,100 | 507 |
| case1_kept, genuinely novel | 11 | 1,387 |
| rv_missing, rescued | 97 | 133 |
| total_cds before correction | 4,077 | 5,367 |

<!-- figure: case-table -->

The one number that differs by two orders of magnitude is case 1, and that is
biology: CDC1551 is the same species as the reference so almost nothing is
genuinely new, while O157:H7 carries prophages and the LEE island that K-12
does not have.

What should worry you is the opposite: **case 2 at or near zero**. Every gene
disagreeing with the reference by a few bases means something systematic is
wrong with the reference, not with 4,000 individual genes. The usual cause is a
reference proteome whose coordinates do not describe the sequence its proteins
were translated from.

### `PANGENE/pangenome.Rtab`

The gene presence/absence matrix of the whole run, the one both reports draw
the pangenome from: one row per gene, one column per sample, 1 where pangene
found the gene in that sample's genome. Every genome is searched with the same
proteins, the reference proteome and the genes no reference gene matches, so a
gene with a 1 in every column is one all the samples share. A gene is named
after its reference gene (`Rv1506c`) or, when no reference gene matches it,
after one sample's copy of it (`L4_CDC1551_01396`, the ID in that sample's
master table).

### `final_report.html`

Self-contained, so it travels by email. Covers per-sample quality, missing
reference genes by lineage, PE/PPE family distribution, KEGG pathways,
duplication communities and network topology, the pangenome (core, softcore,
shell, cloud, frequency spectrum, presence/absence heatmaps), COG and
Mycobrowser distributions, case statistics and the associated tests, plus a
searchable master-table explorer.

### `amap_report.html`

One self-contained file, organised by the questions a run raises: how each
genome annotated, where genes are missing or broken on the reference, one gene
in every genome, what the correction changed, the pseudogenes by lineage, what
lies outside the reference, and the functional profile. Written with or without
the functional annotation. [The report by question](amap-report.md) goes
through it.

## Effect of the skip options

- **`--skip_typing`**: No `PATHOTYPR/` directory. Samples keep their original IDs instead of being
  renamed by lineage.
- **`--skip_functional_annotation`**: The run stops after correction. No `DUPLICATIONS/`, `PANGENE/` or
  `final_report.html`, and `AMAP/` holds only the corrected table and sequences,
  `*.AMAP.tsv`, `*.AMAP.faa` and `*.AMAP.fna`, without the GFF3, the GenBank or
  the master table. You still get `BAKTA/` and `CDS_CASES/`, including
  `*_corrected.tsv` with the classification's coordinates, and `amap_report.html`,
  without its eggNOG sections.
- **`--skip_mycobrowser`**: `*-AMAP.MycoBrowser.tsv` is written empty and the Mycobrowser sections drop
  out of the report. Everything else is unaffected.

## Next steps

- [Overview](../how-it-works/index.md): what produced all this
- [Duplication metrics](../how-it-works/duplications.md): reading the network numbers
