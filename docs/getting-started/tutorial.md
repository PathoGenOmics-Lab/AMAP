# Tutorial: from an empty folder to an annotated genome

Every command spelled out, nothing assumed. It uses a real genome you can
download in one line, so you can follow along and compare your numbers against
the ones printed here. Those numbers come from actual runs of these commands,
not from memory, made on 18 October 2026 with AMAP at commit `947a99d`. The
runs skipped eggNOG, with `--skip_functional_annotation`. None of the numbers
depends on it, and the one table shown here that is written after it was made
by hand from the same run, as [the sample check](#sanity-check-the-sample)
explains.

The genome is *Mycobacterium tuberculosis* CDC1551, corrected against the H37Rv
reference that ships inside the AMAP container. At the end there is a second
pass showing the same thing for a non-MTBC organism.

The correction itself takes minutes on one genome. What dominates the first run
is setup: pulling a 12 GB image and building the eggNOG database, both of which
happen once and are reused afterwards. Expect the first run to take
considerably longer than every run after it, and do not read a long silence
during the database build as a failure.

## 1. What you need

Nextflow and a container engine. Every tool the pipeline itself calls lives
inside the AMAP image, so there is nothing else to install.

The easiest way to get Nextflow is Conda. If you do not have it,
[Miniforge](https://github.com/conda-forge/miniforge) is the one to install:

```bash
conda create -n nextflow-env -c conda-forge -c bioconda nextflow
conda activate nextflow-env
```

Nextflow needs Java 17 or newer. Check what you have:

```bash
java -version
nextflow -version
```

> [!WARNING]
> **Java 8 is not enough**
>
> Nextflow will refuse to start on Java 8 with a message about an
> unsupported class file version. Install a newer JDK, for example
> `conda install -c conda-forge openjdk=17`.

Then a container engine. Docker on a laptop, Singularity or Apptainer on a
cluster:

```bash
docker --version        # or: singularity --version
```

Full details, including running without a container, are on the
[installation page](installation.md).

> [!WARNING]
> **On an Apple Silicon Mac, read this first**
>
> The AMAP image is `linux/amd64` only, so on an M-series Mac every step runs
> under emulation. The test profile in step 2 is fine. A real genome is slow
> enough that a working run can look like a hung one, which is the part worth
> knowing before you start wondering what you broke.
>
> It is a reasonable way to learn what the pipeline does and to check your
> sample sheet. For results you intend to use, move to an amd64 machine or a
> cluster.

### The machine

Worth checking before you start, because this is what stops a run rather than
anything you will have got wrong.

Disk goes first. The container is 12 GB unpacked. The eggNOG database is built
on first use and cached: AMAP skips the 60 GB Diamond core that eggNOG-mapper
would otherwise fetch and builds a bacteria-only index of about 3 GB instead,
on top of the metadata it downloads. Nextflow keeps every intermediate under
`work/`, which for one genome is a few GB and grows with the sample sheet.

Memory decides whether the run finishes at all. Two steps are heavy:

| Step | CPUs | Memory requested |
| :--- | ---: | ---: |
| eggNOG annotation | 24 | 50 GB |
| eggNOG database build | 16 | 45 GB |
| Pangenome | 16 | 16 GB |
| Duplication analysis | 6 to 13 | 10 GB |
| Report | 1 | 8 GB |
| Bakta annotation | 12 | 6 GB |
| Everything else | 1 to 8 | 6 GB or less |

Nextflow will not start a task that asks for more memory than the machine has,
so on a laptop these requests are the thing that fails, not your data. Lower
the ceiling and every request is clamped to it:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --max_cpus 4 --max_memory 12.GB
```

eggNOG will be slower and may still run out on a genuinely small machine. The
alternative is to skip functional annotation altogether:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --skip_functional_annotation
```

You still get the corrected genes, which is what AMAP exists to produce, as a
table with their proteins and sequences, and `amap_report.html`. What you lose
is eggNOG categories, the GFF3 and GenBank files written after it, the
pangenome, the duplication analysis and `final_report.html`. The heaviest step
left asks for 6 GB.

> [!CAUTION]
> **Docker Desktop has its own memory limit, and it is the one that counts**
>
> On macOS and Windows, Docker runs your containers inside a virtual machine
> with a memory ceiling of its own, often 8 GB by default no matter how much
> the machine has. Nextflow sizes its requests against the host, so it will
> happily start a 12 GB task that the VM then kills. You get exit code 137
> and nothing that mentions Docker.
>
> Raising `--max_memory` makes this worse, not better. Raise the VM instead,
> in Docker Desktop under Settings, Resources, Memory, and keep
> `--max_memory` at or below what you set there.
>
> Linux has no such layer. Containers use host memory directly and this does
> not apply.

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

> [!CAUTION]
> **It has to be a real tab**
>
> Spaces will not do, and this is the single most common reason a sample
> sheet is rejected. Check with `cat -t samples.tsv`: real tabs show up as
> `^I`. (`cat -A` does the same on Linux but is not accepted by the `cat`
> that ships with macOS.)

Keep sample IDs to letters, digits and underscores. Several scripts cut the
name at the first dot, so `G26689.1` becomes `G26689`.

## 5. Validate before committing to a long run

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --validate_only
```

The checks write their verdict to `results/pipeline_info/validation.log`, which
ends:

```bash
tail -n 8 results/pipeline_info/validation.log
```

```text
--- FASTA Validation (1 samples) ---

============================================================
Results: 1 OK, 0 warnings, 0 failed
Total: 0 error(s), 0 warning(s)
============================================================

  ✓ PASSED — all 1 samples look good
```

This costs seconds and catches the mistakes that would otherwise surface
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
ls -F results/
```

```text
AMAP/  BAKTA/  CDS_CASES/  DUPLICATIONS/  PANGENE/  PATHOTYPR/
amap_report.html  final_report.html  pipeline_info/
```

A run without eggNOG, like the ones behind this page, has no `DUPLICATIONS/`,
`PANGENE/` or `final_report.html`.

Open `results/final_report.html` first. It is self-contained, so it travels by
email. `results/amap_report.html` goes through the same run question by
question, and a run without eggNOG still gets it.

### Did the correction do anything?

```bash
column -t -s$'\t' results/CDS_CASES/summary/*-summary.tsv
```

```text
Type           Number
case2          2813     coordinates matched the reference exactly
case3            24     coordinates matched, alignment length differed
case4          1100     overlapped the reference but disagreed, resolved
case1_initial   141     Bakta called it, the reference has nothing there
rv_missing       97     the reference has a gene there, Bakta missed it
total_cds      4077     what Bakta called before correction
case1_removed   130     "novel" genes that were not really novel
case1_kept       11     kept as genuinely absent from H37Rv
```

<!-- figure: case-counts -->

Read that as: of 4077 Bakta calls, 2813 already agreed with the reference,
1100 had to be adjudicated, and 97 reference genes Bakta had missed were
rescued. Only 11 of the 141 apparently novel genes survived scrutiny, which is
what you would expect for a strain of the same species as the reference.

### Sanity check the sample

```bash
column -t -s$'\t' results/AMAP/*/*-AMAP.qc.tsv
```

Six of its eleven columns:

```text
genome_size  n_contigs  gc_pct  total_cds  pseudogene_count  coding_density_pct
4403837      1          65.61   4014       123               91.3
```

65.6% GC and 91.3% coding density are what a healthy MTBC genome looks like. A
number well outside that is worth investigating before trusting the rest.

This table is written after eggNOG, so a run without it has none. The line
above was made by running that step's own script, GET_OUTPUTS, on the
corrected genes of the run behind this page, with an eggNOG table that has no
rows. None of the six columns depends on eggNOG; the one that does,
`eggnog_annotation_pct`, is left out.

### The files you will actually use

Lineage typing puts the lineage in front of the sample ID, so CDC1551 comes
out as `L4_CDC1551`:

| What | Where |
| :--- | :--- |
| Corrected annotation | `results/AMAP/L4_CDC1551/L4_CDC1551-AMAP.gff3` |
| Same, GenBank | `results/AMAP/L4_CDC1551/L4_CDC1551-AMAP.gbk` |
| Proteins for downstream tools | `results/AMAP/L4_CDC1551/L4_CDC1551-AMAP.clean.faa` |
| Everything known per gene | `results/AMAP/L4_CDC1551/L4_CDC1551-AMAPmastertable.tsv` |
| Tool versions for your methods | `results/pipeline_info/software_versions.yml` |

> [!NOTE]
> **Two protein files, and they differ on purpose**
>
> `*-AMAP.faa` has one entry per annotated gene, pseudogenes included, so it
> matches the GFF3 gene for gene. `*-AMAP.clean.faa` drops pseudogenes and
> anything with an internal stop codon, which is what orthology,
> phylogenomics and dN/dS tools expect. It is legitimately shorter.

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
case2          3324
case3            31
case4           507
case1_kept     1387    genes O157:H7 has and K-12 does not
rv_missing      133
```

Those 1387 are the prophages and the LEE island that O157:H7 carries and K-12
lacks. Against a reference of your own species you should expect a large case 1
when the strains genuinely differ, and a small one when they do not.

> [!WARNING]
> **The annotation and the genome must be the same assembly**
>
> Coordinates from the annotation are used directly against the FASTA. An
> annotation from a different assembly produces nonsense: proteins that do
> not translate, and a nearly empty reference proteome. AMAP writes how many
> proteins it extracted into the work directory, where
> `grep -h 'Curated proteins' work/*/*/.command.out` finds it (4318 on this
> pair), so check the number looks like a genome's worth.

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
