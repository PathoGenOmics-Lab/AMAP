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

The order does not matter: each profile sets its own kind of thing. The ceilings
on what any one task may ask for come from which profiles are on, not from their
order: 48 CPUs, 256 GB and 120 h with `highmem`, 4 CPUs, 8 GB and 1 h with
`test`, 48 CPUs, 128 GB and 72 h with an executor profile, and otherwise the
machine's CPUs, 16 GB and 72 h. Without an executor profile everything runs on
the local machine.

## On a laptop or workstation

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker --tsv samples.tsv
```

CPUs are auto-detected and memory is capped at 16 GB. On a laptop under Docker
Desktop, the VM's limit, often 8 GB, is the one that counts: lower the cap to
what `docker info` reports, and skip eggNOG, which asks for 50 GB:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --max_cpus 4 --max_memory 8.GB \
    --skip_functional_annotation
```

On a workstation with more, raise it:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --max_cpus 16 --max_memory 64.GB
```

[What each step needs](../getting-started/tutorial.md#the-machine).

A handful of genomes is fine locally. Past about twenty genomes, a cluster is
worth the setup.

## On an HPC cluster

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile singularity,slurm --tsv samples.tsv
```

Each task becomes its own scheduler job, up to 200 at once (the checks, the
versions and the parameters run on the launching node), so Bakta, miniprot and
the correction run side by side across your samples. eggNOG is one job for all
of them. Defaults reach 48 CPUs and 128 GB per job, and every process has a
time limit, so a stuck job ends instead of holding its slot.

To match your site's queues and accounting:

```bash
nextflow run PathoGenOmics-Lab/AMAP \
    -profile singularity,slurm \
    --tsv samples.tsv \
    --max_cpus 64 --max_memory 256.GB \
    -process.queue long \
    -process.clusterOptions '--account=mygroup'
```

[HPC clusters](../deployment/hpc.md) goes into the details that actually
bite: shared filesystems, scratch, and where the eggNOG database should live.

> [!TIP]
> **Large datasets**
>
> Past about twenty genomes, use a scheduler. The eggNOG step in particular
> runs once over all the samples rather than once per sample, precisely
> so a shared filesystem is not hammered by hundreds of jobs reading the same
> Diamond database.

## Non-MTBC organisms

Lineage typing and Mycobrowser categories are specific to the *Mycobacterium tuberculosis* complex (MTBC). Turn
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
`.faa`, `.fna`) and the case tables, and skip eggNOG, the duplications and the
pangenome: no `final_report.html`; `amap_report.html` is still written, without
its eggNOG sections. The GFF3 and GenBank files are written after eggNOG, so
they are skipped too. Much faster when the corrected genes are all you are
after:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --skip_functional_annotation
```

## Pinning the pipeline

These pages describe AMAP as it is on `main`, which moves. Without `-r`,
Nextflow runs the copy of AMAP it downloaded first and keeps it until you
update it, so two machines can run different code under the same command. For
an analysis you will publish, pin both the code and the image:

```bash
nextflow run PathoGenOmics-Lab/AMAP \
    -r 947a99dee49d2a2429e014296ac35ba08e911108 \
    -profile docker --tsv samples.tsv \
    -process.container 'paururo/amap@sha256:08bd17311ab323317ceecba6e5bc639a29f4c54c9d95f9a7c48ba80a9fc640ee'
```

`-r` takes a branch, a tag or a commit. Give a commit as its full hash:
Nextflow 26.10 refuses a short one, with "Remote does not have 947a99d
available for fetch", where 24.04 takes it. The image is pinned by its digest,
which never moves, where a tag such as `latest` can; this one is the image
behind the runs on these pages, and
`docker image inspect paururo/amap:latest --format '{{index .RepoDigests 0}}'`
prints the one you have.

To move to the newest code, `nextflow pull PathoGenOmics-Lab/AMAP`. To see
which code a run used, open `pipeline_info/execution_report_*.html`: it names
the revision and its commit hash, as in `revision main (commit hash 1a1b3959...)`.

## What a run produces

Beyond the results themselves, `results/pipeline_info/` holds the execution
report, the timeline, the DAG, a trace with per-process resource usage,
`software_versions.yml` with the exact version of every tool, and `params.yml`
with the parameters used. Those last two are what you want when writing up
methods.

## Next steps

- [Outputs](outputs.md): what every file is
- [Troubleshooting](../deployment/troubleshooting.md): when a run fails
