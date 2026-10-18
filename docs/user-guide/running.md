# Running AMAP

Nextflow runs the pipeline straight from GitHub, so there is nothing to clone
and nothing to keep up to date by hand.

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker --tsv samples.tsv
```

## Profiles

A profile is a bundle of settings. You normally combine one **container**
profile with one **executor** profile.

| Profile | Kind | What it does |
| :--- | :--- | :--- |
| `docker` | container | Runs every process in the AMAP Docker image. |
| `singularity` | container | Same image via Singularity or Apptainer. |
| `standard` | container | No container. Every tool must already be on `PATH`. |
| `slurm` | executor | Submits each process as a SLURM job. |
| `sge` | executor | Submits to Sun Grid Engine or Univa Grid Engine. |
| `pbs` | executor | Submits to PBS or Torque. |
| `lsf` | executor | Submits to IBM Spectrum LSF. |
| `test` | data | Small bundled dataset, for checking the installation. |
| `highmem` | resources | Raises the ceilings to 48 CPUs, 256 GB, 120 h. |

Combine them with commas:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile singularity,slurm --tsv samples.tsv
```

Write the container profile first. From Nextflow 25.10 on, where two profiles
set the same thing the later one wins, and the container profiles cap every task
at 16 GB, which suits a laptop: `-profile slurm,singularity` would hold a cluster
run to that, where `-profile singularity,slurm` gets the cluster's 128 GB.
Without an executor profile everything runs on the local machine.

## On a laptop or workstation

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker --tsv samples.tsv
```

CPUs are auto-detected and memory is capped at 16 GB. Raise both if the machine
allows:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --max_cpus 8 --max_memory 32.GB
```

A handful of genomes is fine locally. Beyond about twenty, a cluster is worth
the setup.

## On an HPC cluster

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile singularity,slurm --tsv samples.tsv
```

Each process becomes its own scheduler job, with up to 200 running at once, so
Bakta, Miniprot and eggNOG parallelise across all your samples. Defaults reach
48 CPUs and 128 GB per job, and every process carries a time limit so nothing
sits in the queue forever.

To match your site's queues and accounting:

```bash
nextflow run PathoGenOmics-Lab/AMAP \
    -profile singularity,slurm \
    --tsv samples.tsv \
    --max_cpus 64 --max_memory 256.GB \
    -process.queue long \
    -process.clusterOptions '--account=mygroup'
```

[Deployment on HPC](../deployment/hpc.md) goes into the details that actually
bite: shared filesystems, scratch, and where the eggNOG database should live.

> [!TIP]
> **Large datasets**
>
> Past a hundred genomes, use a scheduler. The eggNOG step in particular runs
> once over all samples at a time rather than once per sample, precisely so a
> shared filesystem is not hammered by hundreds of jobs reading the same
> Diamond database.

## Non-MTBC organisms

Lineage typing and Mycobrowser categories are specific to tuberculosis. Turn
them off and bring your own reference:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv \
    --skip_typing --skip_mycobrowser \
    --genus Escherichia --species "Escherichia coli" \
    --ref_annotation my_reference.gbk --ref_genome my_reference.fasta
```

See [custom reference](custom-reference.md).

## Resuming

Nextflow caches every completed task, so an interrupted run picks up where it
stopped:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile singularity,slurm \
    --tsv samples.tsv -resume
```

Keep the parameters identical. Changing an input invalidates that task and
everything downstream, which is correct behaviour but does mean the work is
redone. The cache lives in the work directory, so do not delete it between a run
and its resume.

## Running only part of the pipeline

`--skip_functional_annotation` stops after the correction step. You get the
corrected genes as a table with their proteins and sequences (`*.AMAP.tsv`,
`.faa`, `.fna`) and the case tables, and skip eggNOG, duplications, the
pangenome and the report. The GFF3 and GenBank files are written after eggNOG,
so they are skipped too. Much faster when the corrected genes are all you are
after:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --skip_functional_annotation
```

## What a run produces

Beyond the results themselves, `results/pipeline_info/` holds the execution
report, the timeline, the DAG, a trace with per-process resource usage,
`software_versions.yml` with the exact version of every tool, and `params.yml`
with the parameters used. Those last two are what you want when writing up
methods.

## Next steps

- [Outputs](outputs.md): what every file is
- [Troubleshooting](../deployment/troubleshooting.md): when a run fails
