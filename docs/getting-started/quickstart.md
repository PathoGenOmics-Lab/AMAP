# Quickstart

This walks from an installed Nextflow to an open report. If Nextflow or a
container engine is still missing, start at [installation](installation.md).

## 1. Write a sample sheet

AMAP reads a tab-separated file with two columns, `sampleId` and `fasta`:

```tsv title="samples.tsv"
sampleId	fasta
G26689	/data/genomes/G26689.fasta
G26690	/data/genomes/G26690.fasta
```

The separator must be a real tab. Paths can be absolute, or relative to where
you launch Nextflow. See [preparing inputs](inputs.md) for the details and the
checks AMAP runs over them.

## 2. Check the inputs

Cheap, and it catches the mistakes that would otherwise surface twenty minutes
into a run:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --validate_only
```

## 3. Run it

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

Resource requests are capped at 24 CPUs and 50 GB by default, which is what the
eggNOG annotation step asks for. Nextflow refuses to start a task requesting
more memory than the machine has, so on anything smaller than a server, lower
the ceiling:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --max_cpus 4 --max_memory 12.GB
```

The [tutorial](tutorial.md#the-machine) has the per-step figures and what to do
when the machine cannot reach them.

## 4. Read the results

Everything lands under `results/` unless you set `--outdir`. The three things to
open first:

| What | Where |
| :--- | :--- |
| The full interactive report | `results/final_report.html` |
| Corrected annotation, per sample | `results/AMAP/{sampleId}/*-AMAP.gff3` |
| Protein sequences for downstream tools | `results/AMAP/{sampleId}/*-AMAP.clean.faa` |

> [!NOTE]
> **Which protein file do I want?**
>
> `*-AMAP.faa` holds every gene, pseudogenes included, so it mirrors the
> annotation exactly. `*-AMAP.clean.faa` drops pseudogenes and sequences with
> internal stop codons, which is what orthology, phylogenomics and dN/dS
> tools expect. Reach for the clean one unless you specifically need the
> complete set.

The [outputs page](../user-guide/outputs.md) covers every directory.

## 5. Resume an interrupted run

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

**Just the corrected annotation.** If you do not need functional annotation,
duplications, the pangenome or the report, stop after the correction step. It is
much faster:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --skip_functional_annotation
```

## Next steps

- [Running AMAP](../user-guide/running.md): profiles, schedulers and larger runs
- [Parameters](../user-guide/parameters.md): the full reference
- [How it works](../how-it-works/index.md): what the correction actually does
