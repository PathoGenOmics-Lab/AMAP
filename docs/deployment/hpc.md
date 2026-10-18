# Running on a cluster

Past a couple of dozen genomes, a scheduler is worth the setup. Combine a
container profile with an executor profile:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile singularity,slurm --tsv samples.tsv
```

<!-- figure: hpc -->

## The executor profiles

| Profile | Scheduler | Submission rate |
| :--- | :--- | :--- |
| `slurm` | SLURM | 100 per minute |
| `sge` | Sun Grid Engine, Univa | 20 per minute |
| `pbs` | PBS, Torque | 20 per minute |
| `lsf` | IBM Spectrum LSF | 20 per minute |

All four queue up to 200 jobs at once, default to 48 CPUs and 128 GB per job,
and enable node-local scratch. SGE additionally requests memory through
`-l h_vmem`, and uses the `smp` parallel environment.

Tune the concurrency to your site's policy if 200 is too aggressive:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile singularity,slurm \
    --tsv samples.tsv -process.queueSize 50
```

## Fitting your site

```bash
nextflow run PathoGenOmics-Lab/AMAP \
    -profile singularity,slurm \
    --tsv samples.tsv \
    --max_cpus 64 --max_memory 256.GB \
    -process.queue long \
    -process.clusterOptions '--account=mygroup'
```

`-process.queue` sets the partition or queue, and `-process.clusterOptions`
passes anything else straight through, which is where accounting flags, QoS and
node constraints go.

On SLURM, Bakta's jobs, the longest a run sends, also ask for `--qos=medium`.
The `slurm` profile sets that for Bakta alone, and for Bakta it takes the place
of `-process.clusterOptions`, so give Bakta your site's options in a config
passed with `-c`:

```groovy
process {
    withName: 'BAKTA_ANNOTATOR' {
        clusterOptions = '--qos=medium --account=mygroup'
    }
}
```

## Scratch

`--scratch` is on by default in every HPC profile, and it matters. Without it,
every process reads and writes directly on the shared filesystem, and Bakta and
eggNOG in particular generate enough metadata traffic to slow a cluster down for
everyone. With it, work happens on node-local disk and only the results are
copied back.

Turn it off only if your nodes have no local disk:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile singularity,slurm \
    --tsv samples.tsv --scratch false
```

## Where to put things

**The Singularity cache.** Set it once, somewhere shared:

```bash
export NXF_SINGULARITY_CACHEDIR=/shared/singularity_cache
```

Otherwise the `singularity` profile caches under `--outdir`, and a large run can
have many jobs converting the same image at once.

**The eggNOG database.** Large, built once, reused forever. Point it at
persistent shared storage rather than somewhere per-run:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile singularity,slurm \
    --tsv samples.tsv --eggnog_db /shared/databases/eggnog
```

**The work directory.** It gets large, tens of gigabytes for a big run. Put it
on a filesystem with room and a generous quota:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile singularity,slurm \
    --tsv samples.tsv -work-dir /scratch/$USER/amap_work
```

Keep it until you are sure the run is finished, since `-resume` needs it.

## Launching

Nextflow's own process must survive the whole run, which will outlive an SSH
session. Either run it under `tmux` or `screen`, or submit it as a small job of
its own:

```bash title="run_amap.sh"
#!/bin/bash
#SBATCH --job-name=amap
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=amap_%j.log

module load singularity
export NXF_SINGULARITY_CACHEDIR=/shared/singularity_cache

nextflow run PathoGenOmics-Lab/AMAP \
    -profile singularity,slurm \
    --tsv samples.tsv \
    --eggnog_db /shared/databases/eggnog \
    -work-dir /scratch/$USER/amap_work \
    -resume
```

The launcher itself needs very little: it submits and waits while the real work
happens in the jobs it spawns. Give it plenty of wall time, though, since it has
to outlast everything it submits.

## Memory and retries

A process killed by the out-of-memory killer is resubmitted with more memory,
scaling with the attempt number up to `--max_memory`. That handles the occasional
unusually large genome without you having to size everything for the worst case.

eggNOG is the memory-hungry step. It runs once over all samples rather than once
per sample, which keeps a shared filesystem from being hammered but does mean
one large job. The `highmem` profile exists for this:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile singularity,slurm,highmem \
    --tsv samples.tsv
```

## Watching a run

`squeue -u $USER` shows the jobs. When it finishes, `pipeline_info/` holds the
execution report, the timeline, and a trace with per-process CPU, memory and
runtime, which is what to look at before adjusting resources for the next run.

## Next steps

- [Troubleshooting](troubleshooting.md): when something fails
- [Containers](containers.md): what the image holds
