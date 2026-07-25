---
title: Tutorial
---

# Tutorial: from an empty folder to an annotated genome

Every command spelled out, nothing assumed. It uses a real genome you can
download in one line, so you can follow along and compare your numbers against
the ones printed here. Those numbers come from an actual run, not from memory.

The genome is *Mycobacterium tuberculosis* CDC1551, corrected against the H37Rv
reference that ships inside the AMAP container. At the end there is a second
pass showing the same thing for a non-MTBC organism.

Budget about half an hour of wall clock, most of it spent downloading the
container and building the eggNOG database on first use.

## 1. What you need

Two things, and nothing else. Every tool the pipeline calls lives inside the
AMAP image.

```bash
conda create -n nextflow-env
conda activate nextflow-env
conda install -c bioconda nextflow
```

Nextflow needs Java 11 or newer. Check what you have:

```bash
java -version
nextflow -version
```

!!! warning "Java 8 is not enough"

    Nextflow will refuse to start on Java 8 with a message about an
    unsupported class file version. Install a newer JDK, for example
    `conda install -c conda-forge openjdk=17`.

Then a container engine. Docker on a laptop, Singularity or Apptainer on a
cluster:

```bash
docker --version        # or: singularity --version
```

Full details, including running without a container, are on the
[installation page](installation.md).

## 2. Check the plumbing before touching real data

The test profile annotates the reference genome that ships inside the
container, so it needs no data of your own:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker,test
```

The first run pulls the AMAP image. It is around 3.4 GB compressed and 12 GB
on disk, because it carries the Bakta database and the H37Rv reference files.
That is a one-off; later runs reuse it.

If this finishes without an error, the whole chain works.

## 3. Get a genome

Make a working directory and fetch an assembly:

```bash
mkdir amap-tutorial && cd amap-tutorial

NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF
ASM=000/008/585/GCF_000008585.1_ASM858v1/GCF_000008585.1_ASM858v1_genomic

curl -sL -o cdc1551.fna.gz "$NCBI/$ASM.fna.gz"
gunzip cdc1551.fna.gz
```

Have a look at what you got:

```bash
grep '^>' cdc1551.fna
```

```text
>NC_002755.2 Mycobacterium tuberculosis CDC1551, complete sequence
```

One closed contig, 4,403,837 bp. AMAP takes assembled genomes, so reads are not
an input here.

## 4. Write the sample sheet

A tab-separated file with exactly two columns:

```bash
printf 'sampleId\tfasta\nCDC1551\t%s/cdc1551.fna\n' "$PWD" > samples.tsv
cat samples.tsv
```

```text
sampleId	fasta
CDC1551	/path/to/amap-tutorial/cdc1551.fna
```

!!! danger "It has to be a real tab"

    Spaces will not do, and this is the single most common reason a sample
    sheet is rejected. Check with `cat -t samples.tsv`: real tabs show up as
    `^I`. (`cat -A` does the same on Linux but is not accepted by the `cat`
    that ships with macOS.)

Keep sample IDs to letters, digits and underscores. Several scripts cut the
name at the first dot, so `G26689.1` becomes `G26689`.

## 5. Validate before committing to a long run

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --validate_only
```

```text
--- FASTA Validation (1 samples) ---

============================================================
Results: 1 OK, 0 warnings, 0 failed
Total: 0 error(s), 0 warning(s)
============================================================

  ✓ PASSED — all 1 samples look good
```

This costs seconds and catches the mistakes that would otherwise surface twenty
minutes in. [Preparing inputs](inputs.md) lists every check and its threshold.

## 6. Run it

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker --tsv samples.tsv
```

What happens, in order: lineage typing, Bakta annotation, alignment of the
H37Rv proteome onto your genome, case classification, correction, eggNOG
functional annotation, duplication analysis, pangenome, report.

The first run also downloads and builds the eggNOG database, which is several
gigabytes and dominates the wall clock. Point `--eggnog_db` at somewhere
persistent so the next run reuses it:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --eggnog_db ~/databases/eggnog
```

If the run dies partway, `-resume` picks up from the last completed step.

## 7. Look at what came out

```bash
ls results/
```

```text
AMAP/  BAKTA/  CDS_CASES/  DUPLICATIONS/  PANGENE/  PATHOTYPR/
pipeline_info/  final_report.html
```

Open `results/final_report.html` first. It is self-contained, so it travels by
email.

### Did the correction do anything?

```bash
column -t -s$'\t' results/CDS_CASES/summary/*-summary.tsv
```

```text
Type           Number
case2          2816     coordinates matched the reference exactly
case3            24     coordinates matched, alignment length differed
case4          1097     overlapped the reference but disagreed, resolved
case1_initial   141     Bakta called it, the reference has nothing there
rv_missing       97     the reference has a gene there, Bakta missed it
total_cds      4077     what Bakta called before correction
case1_removed   130     "novel" genes that were not really novel
case1_kept       11     kept as genuinely absent from H37Rv
```

Read that as: of 4077 Bakta calls, 2816 already agreed with the reference,
1097 had to be adjudicated, and 97 reference genes Bakta had missed were
rescued. Only 11 of the 141 apparently novel genes survived scrutiny, which is
what you would expect for a strain of the same species as the reference.

### Sanity check the sample

```bash
column -t -s$'\t' results/AMAP/*/*-AMAP.qc.tsv
```

```text
genome_size  n_contigs  gc_pct  total_cds  pseudogene_count  coding_density_pct
4403837      1          65.61   4020       13                91.5
```

65.6% GC and 91.5% coding density are what a healthy MTBC genome looks like. A
number well outside that is worth investigating before trusting the rest.

### The files you will actually use

| What | Where |
| :--- | :--- |
| Corrected annotation | `results/AMAP/CDC1551/CDC1551-AMAP.gff3` |
| Same, GenBank | `results/AMAP/CDC1551/CDC1551-AMAP.gbk` |
| Proteins for downstream tools | `results/AMAP/CDC1551/CDC1551-AMAP.clean.faa` |
| Everything known per gene | `results/AMAP/CDC1551/CDC1551-AMAPmastertable.tsv` |
| Tool versions for your methods | `results/pipeline_info/software_versions.yml` |

!!! note "Two protein files, and they differ on purpose"

    `*-AMAP.faa` has one entry per annotated gene, pseudogenes included, so it
    matches the GFF3 gene for gene. `*-AMAP.clean.faa` drops pseudogenes and
    anything with an internal stop codon, which is what orthology,
    phylogenomics and dN/dS tools expect. It is legitimately shorter.

Every directory is documented under [outputs](../user-guide/outputs.md).

## 8. More than one genome

Add rows. The pipeline parallelises across samples, and the pangenome and the
cross-genome duplication analysis only become meaningful with several:

```tsv
sampleId	fasta
CDC1551	/data/cdc1551.fna
Erdman	/data/erdman.fna
Beijing_1	/data/beijing_1.fna
```

Past a couple of dozen, move to a cluster:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile singularity,slurm \
    --tsv samples.tsv --eggnog_db /shared/databases/eggnog
```

[Deployment on HPC](../deployment/hpc.md) covers scratch, shared filesystems
and where things should live.

## 9. The same thing for another organism

Everything above assumed MTBC. For anything else you supply the trusted
reference yourself, and turn off the two MTBC-specific steps.

Fetch a reference and a sample. Here, *E. coli* O157:H7 Sakai corrected against
K-12 MG1655:

```bash
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF
K12=000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic
SAKAI=000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic

# the reference: annotation plus the genome it describes
curl -sL -o ref.gbk.gz   "$NCBI/$K12.gbff.gz"
curl -sL -o ref.fna.gz   "$NCBI/$K12.fna.gz"

# the sample
curl -sL -o sakai.fna.gz "$NCBI/$SAKAI.fna.gz"

gunzip ref.gbk.gz ref.fna.gz sakai.fna.gz
printf 'sampleId\tfasta\nSakai\t%s/sakai.fna\n' "$PWD" > ecoli.tsv
```

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv ecoli.tsv \
    --skip_typing --skip_mycobrowser \
    --genus Escherichia --species "Escherichia coli" \
    --ref_annotation ref.gbk --ref_genome ref.fna
```

`--skip_typing` because lineage typing is specific to tuberculosis, and
`--skip_mycobrowser` because those functional categories are an MTBC resource.

Before annotating anything, AMAP extracts the proteome from your reference,
annotates that same reference genome with Bakta, and takes the difference to
work out which genes Bakta invents. For MTBC that list was curated by hand;
here it is derived automatically. On this pair it found 142.

The case counts come out differently from the MTBC run, and the difference is
biology rather than a problem:

```text
case2          3325
case3            34
case4           537
case1_kept     1356    genes O157:H7 has and K-12 does not
rv_missing      102
```

Those 1356 are the prophages and the LEE island that O157:H7 carries and K-12
lacks. Against a reference of your own species you should expect a large case 1
when the strains genuinely differ, and a small one when they do not.

!!! warning "The annotation and the genome must be the same assembly"

    Coordinates from the annotation are used directly against the FASTA. An
    annotation from a different assembly produces nonsense: proteins that do
    not translate, and a nearly empty reference proteome. AMAP prints how many
    proteins it extracted, so check the number looks like a genome's worth.

[Custom reference](../user-guide/custom-reference.md) goes into what makes a
good reference and what changes when you use one.

## When something goes wrong

Nextflow prints the work directory of the failed task. Everything you need is
in it:

```text
work/ab/cd1234.../
├── .command.sh    the exact script that ran
├── .command.err   stderr on its own, almost always the answer
└── .exitcode
```

[Troubleshooting](../deployment/troubleshooting.md) is organised by the error
you are looking at.

## Where to go next

- [Parameters](../user-guide/parameters.md): every flag and its default
- [How it works](../how-it-works/index.md): what the correction actually does, and the thresholds behind it
- [Outputs](../user-guide/outputs.md): every file the run produces
