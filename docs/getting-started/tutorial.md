# Tutorial: from an empty folder to an annotated genome

Every command spelled out, nothing assumed. It uses a real genome you can
download in one line, so you can follow along and compare your numbers against
the ones printed here. They come from actual runs of these commands, made on
18 October 2026 with AMAP at commit `947a99d`. For exactly these numbers, add
`-r 947a99dee49d2a2429e014296ac35ba08e911108` to each `nextflow run`
([pinning the pipeline](../user-guide/running.md#pinning-the-pipeline) says
why the whole hash).

The runs skipped eggNOG, with `--skip_functional_annotation`, and no number
here depends on it. The one table eggNOG's step writes, in
[the sample check](#sanity-check-the-sample), was made by hand from the same
run.

The genome is *Mycobacterium tuberculosis* CDC1551, corrected against the H37Rv
reference that ships inside the AMAP container. At the end there is a second
pass showing the same thing for a non-MTBC organism.

The correction itself takes minutes on one genome: 2 min 16 s here, without
eggNOG. What dominates the first full run is setup: pulling a 12 GB image and
building the eggNOG database, both of which happen once and are reused
afterwards. Expect the first run to take considerably longer than every run
after it, and do not read a long silence during the database build as a
failure.

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
cluster. With Docker, ask it what your containers can have:

```bash
docker info --format 'Docker can give containers {{.NCPU}} CPUs and {{.MemTotal}} bytes'
```

```text
Docker can give containers 14 CPUs and 8217014272 bytes
```

That is the machine behind this page: an Apple M4 Pro, with Docker Desktop
given 8 GB of its 24. The command fails when the Docker daemon is not
running, which `docker --version` does not catch. Under 8 GB? Raise it in
Docker Desktop, Settings, Resources. On a cluster, `singularity --version`.

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
on first use and cached, and it is the largest thing a run keeps: `eggnog.db`
alone unpacks to 41.4 GB from a 6.8 GB download. For its bacteria-only Diamond
index AMAP also downloads the 9.5 GB eggNOG proteome, and deletes it once the
index is built. Nextflow keeps every intermediate under `work/`, which came to
97 MB for the test in step 2 and 248 MB for all the runs of this page, none
with eggNOG, and grows with the sample sheet.

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
the ceiling to what `docker info` reports, and every request is clamped to it:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --max_cpus 4 --max_memory 8.GB
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

The test profile annotates the H37Rv genome that ships inside the container,
so it needs no data of your own, and `--skip_functional_annotation` keeps it to
the correction, with no eggNOG database to build:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker,test \
    --skip_functional_annotation --outdir amap-test
```

The first run pulls the AMAP image. It is around 3.4 GB compressed and 12 GB
on disk, because it carries the Bakta database and the H37Rv reference files.
That is a one-off; later runs reuse it. With the image already pulled, the test
took 2 min 47 s on the machine above, and wrote `amap-test/` and `work/` in
the folder it ran from.

If the whole chain works, it ends like this, and then names the work
directory and the command:

```text
=======================================
✓ COMPLETED
=======================================
Pipeline    : PathoGenOmics-Lab/AMAP v1.0.1
Nextflow    : 26.10.2
Profile     : docker,test
Duration    : 2m 47s
Succeeded   : 10
Failed      : 0
Ignored     : 0
Output      : amap-test
```

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
functional annotation, duplication analysis, pangenome, the two reports.

The first run also downloads and builds the eggNOG database, `eggnog.db` alone
41.4 GB, and that dominates the wall clock. Point `--eggnog_db` at somewhere
persistent so the next run reuses it:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --eggnog_db ~/databases/eggnog
```

If the run dies partway, `-resume` picks up from the last completed step.

### How long it took

Without eggNOG, on the machine above, the CDC1551 run took 2 min 16 s from
start to finish. Each step, from the run's execution trace in
`results/pipeline_info/`, its `realtime` and `peak_rss`:

| Step | Time | Peak memory |
| :--- | ---: | ---: |
| `PATHOTYPR`, lineage typing | 1.1 s | 6.1 MB |
| `BAKTA_ANNOTATOR` | 1 min 46 s | 3.5 GB |
| `MINIPROT` | 2 s | 324 MB |
| `TYPE_CASES`, the classification | 14.2 s | 219.2 MB |
| `FIX_CDS`, the correction | 3.7 s | 165.2 MB |
| `AMAP_REPORT` | 944 ms | 35.6 MB |

The run's other four tasks, which set up the reference and record the
parameters, the versions and the validation, took under 4 s each. Bakta is the
step to plan for: the *E. coli* run of step 9 took 5 min 47 s, of which Bakta
on Sakai took 2 min 28 s at 5.7 GB, and `PREPARE_REFERENCE`, which runs Bakta on
the reference itself, 2 min 36 s.

## 7. Look at what came out

```bash
ls -1F results/
```

```text
AMAP/
BAKTA/
CDS_CASES/
PATHOTYPR/
amap_report.html
pipeline_info/
```

That is a run without eggNOG, like the ones behind this page. A full run adds
`DUPLICATIONS/`, `PANGENE/` and `final_report.html`.

Open `results/amap_report.html` first: every run writes it, and it goes through
the run question by question. After a full run, `results/final_report.html`
has the whole analysis. Both are self-contained, so they travel by email.

### Did the correction do anything?

Every file of the sample carries its name as lineage typing leaves it,
`L4_CDC1551`. The case summary's first two columns are the counts:

```bash
cut -f1,2 results/CDS_CASES/summary/L4_CDC1551-summary.tsv | column -t
```

```text
Type           Number
case2          2813
case3          24
case4          1100
case1_initial  141
rv_missing     97
total_cds      4077
case1_removed  130
case1_kept     11
```

The other three, `Sequence`, `Lineage` and `Sample`, say whose counts they are,
which matters once the folder holds several samples. What each row means:

| Case | Count | What it means |
| :--- | ---: | :--- |
| `case2` | 2813 | Coordinates matched the reference exactly |
| `case3` | 24 | Coordinates matched, the alignment's length differed |
| `case4` | 1100 | Overlapped the reference but disagreed, and was resolved |
| `case1_initial` | 141 | Bakta called it, the reference has nothing there |
| `rv_missing` | 97 | The reference has a gene there, Bakta missed it |
| `total_cds` | 4077 | What Bakta called before the correction |
| `case1_removed` | 130 | "Novel" genes that were not really novel |
| `case1_kept` | 11 | Kept as genuinely absent from H37Rv |

<!-- figure: case-counts -->

Read that as: of 4077 Bakta calls, 2813 already agreed with the reference,
1100 had to be adjudicated, and 97 reference genes Bakta had missed were
rescued. Only 11 of the 141 apparently novel genes survived scrutiny, which is
what you would expect for a strain of the same species as the reference.

Cases 2, 3, 4 and the initial case 1 add up to 4,078, one more than the calls.
One call, `L4_CDC1551_03643`, has the same ends as two reference proteins:
Rv3478 at its own length, which makes it case 2, and Rv1361c, a PPE protein
three residues longer, which makes it case 3 too.

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

Past about twenty genomes, move to a cluster:

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
biology rather than a problem. Without lineage typing the sample keeps its
name:

```bash
cut -f1,2 results/CDS_CASES/summary/Sakai-summary.tsv | column -t
```

```text
Type           Number
case2          3324
case3          31
case4          507
case1_initial  1506
rv_missing     133
total_cds      5367
case1_removed  122
case1_kept     1387
```

`case1_kept` counts the genes O157:H7 has and K-12 does not. Those 1387 are
the prophages and the LEE island that O157:H7 carries and K-12 lacks. Against a reference of your own species you should expect a large case 1
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
