# Quickstart

This walks from an installed Nextflow to an open report. If Nextflow or a
container engine is still missing, start at [installation](installation.md).

## 1. Check the plumbing

The test profile annotates the H37Rv genome that ships inside the image, so it
needs no data of your own, and `--skip_functional_annotation` keeps it to the
correction, with no eggNOG database to build:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker,test \
    --skip_functional_annotation --outdir amap-test
```

It writes `amap-test/` and `work/` where you run it, and ends with a
`✓ COMPLETED` block. On an Apple M4 Pro, under emulation and with the image
already pulled, it took 2 min 47 s; the first run also pulls the image, about
3.4 GB compressed and 12 GB on disk, once.

## 2. Write a sample sheet

AMAP reads a tab-separated file with two columns, `sampleId` and `fasta`:

```tsv title="samples.tsv"
sampleId	fasta
G26689	/data/genomes/G26689.fasta
G26690	/data/genomes/G26690.fasta
```

The separator must be a real tab. Paths can be absolute, or relative to where
you launch Nextflow. See [preparing inputs](inputs.md) for the details and the
checks AMAP runs over them.

## 3. Check the inputs

Cheap, and it catches the mistakes that would otherwise surface twenty minutes
into a run:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --validate_only
```

When the check passes, its verdict and any warnings are in
`results/pipeline_info/validation.log`, and the terminal says only that the
step ran. When it fails, the run stops and the terminal prints the whole report.

## 4. Run it

Pick what you have, and copy the command:

<!-- figure: builder -->

Or start from one of these:

**Docker**

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker --tsv samples.tsv
```

**Singularity**

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile singularity --tsv samples.tsv
```

**SLURM cluster**

```bash
nextflow run PathoGenOmics-Lab/AMAP \
    -profile singularity,slurm --tsv samples.tsv
```

Without a scheduler profile every task is capped at the machine's CPUs and
16 GB, and with the scheduler profiles at 48 CPUs and 128 GB, whichever order
the profiles come in; a larger request is lowered to the cap. Nextflow refuses
to start a task asking for more memory than the machine has, so on a smaller
machine, lower the ceiling to what `docker info` says your containers can have,
and leave eggNOG out: it loads its whole database into memory, about 45 GB,
whatever the ceiling.

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --max_cpus 4 --max_memory 8.GB \
    --skip_functional_annotation
```

The [tutorial](tutorial.md#the-machine) has the per-step figures and what to do
when the machine cannot reach them.

## 5. Read the results

Everything lands under `results/` unless you set `--outdir`. Each sample's
folder is named `{lineage}_{sampleId}`: lineage typing puts the lineage in
front, so `G26689` comes out as `L4_G26689`, and under `--skip_typing` it is the
sample ID alone ([outputs](../user-guide/outputs.md) has the details). The
things to open first:

| What | Where |
| :--- | :--- |
| Each genome and each gene, question by question, in every run | `results/amap_report.html` |
| The full report, after a run with eggNOG | `results/final_report.html` |
| Corrected annotation, per sample, after a run with eggNOG | `results/AMAP/{lineage}_{sampleId}/*-AMAP.gff3` |
| Protein sequences for downstream tools, after a run with eggNOG | `results/AMAP/{lineage}_{sampleId}/*-AMAP.clean.faa` |
| The corrected genes, their proteins and sequences, in every run | `results/AMAP/{lineage}_{sampleId}/*.AMAP.tsv`, `.faa`, `.fna` |

> [!NOTE]
> **Which protein file do I want?**
>
> `*-AMAP.faa` holds every gene, pseudogenes included, so it mirrors the
> annotation exactly. `*-AMAP.clean.faa` drops pseudogenes and sequences with
> internal stop codons, which is what orthology, phylogenomics and dN/dS
> tools expect. Reach for the clean one unless you specifically need the
> complete set.

The [outputs page](../user-guide/outputs.md) covers every directory.

## 6. Resume an interrupted run

Nextflow caches completed tasks. If a run dies, restart it with `-resume` and it
picks up where it stopped:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker --tsv samples.tsv -resume
```

Keep the same parameters. Changing an input invalidates the cache for every task
downstream of it, which is what you want, but it does mean the run redoes that
work.

## Common variations

**A non-MTBC organism.** Lineage typing and Mycobrowser are specific to
tuberculosis, so turn them off and supply your own reference:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv \
    --skip_typing --skip_mycobrowser \
    --genus Escherichia --species "Escherichia coli" \
    --ref_annotation my_reference.gbk --ref_genome my_reference.fasta
```

See [custom reference](../user-guide/custom-reference.md).

**Just the corrected genes.** If you do not need functional annotation,
duplications or the pangenome, stop after the correction step. It is much
faster, and leaves the corrected genes as a table with their proteins and
sequences; the GFF3 and GenBank files come after eggNOG, so there are none. There
is no `final_report.html`; `amap_report.html` is still written, without its
eggNOG sections:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --skip_functional_annotation
```

## Next steps

- [Running AMAP](../user-guide/running.md): profiles, schedulers and larger runs
- [Parameters](../user-guide/parameters.md): the full reference
- [How it works](../how-it-works/index.md): what the correction actually does
