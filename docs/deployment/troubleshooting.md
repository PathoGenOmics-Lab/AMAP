---
title: Troubleshooting
---

# Troubleshooting

## Reading a failure

When a process fails, Nextflow prints the work directory for the failed task.
That directory has everything you need:

```text
work/ab/cd1234.../
├── .command.sh      # the exact script that ran
├── .command.log     # everything it printed
├── .command.err     # stderr on its own
└── .exitcode        # what it exited with
```

`.command.err` is almost always where the answer is. Re-running `.command.sh`
inside the container reproduces the failure by hand.

`results/pipeline_info/execution_trace_*.txt` shows every task with its status,
CPU, memory and runtime, which is the fastest way to see whether a failure was
resources or logic.

## Before the pipeline starts

**"Unexpected input: '('" in `nextflow.config`**
: Nextflow 25.10 and later read pipelines with a stricter parser, which older
  copies of AMAP do not pass. Nextflow keeps the copy it downloaded first, so
  fetch the current one with `nextflow pull PathoGenOmics-Lab/AMAP` and run
  again. To run an older AMAP release as it is, `export NXF_SYNTAX_PARSER=v1`
  first.

**"Input TSV file not found"**
: The path in `--tsv` is wrong, or relative to somewhere other than where you
  launched Nextflow.

**"Invalid TSV format"**
: The header must be exactly `sampleId` and `fasta`, separated by a real tab.
  Spaces are the usual culprit. `cat -t samples.tsv` shows tabs as `^I`, and
  works on both Linux and macOS.

**"FASTA file not found"**
: A path in the sample sheet does not resolve. Absolute paths avoid the
  ambiguity.

**Validation failures**
: `results/pipeline_info/validation.log` says which check failed and for which
  sample. Thresholds are in [preparing inputs](../getting-started/inputs.md).

## Container problems

**The image will not pull.** Check network access from the compute nodes, which
often differs from the login node. Pre-pulling on the login node and setting
`NXF_SINGULARITY_CACHEDIR` sidesteps it:

```bash
export NXF_SINGULARITY_CACHEDIR=/shared/singularity_cache
singularity pull docker://paururo/amap:latest
```

**Slow on Apple Silicon.** The image is amd64 only, so it emulates. Fine for the
test profile, painful beyond that. Use an amd64 machine.

**Docker permission denied.** Your user is not in the `docker` group. Either add
it, or use Singularity.

## Memory

**A process is killed with exit 137 or 140.** Out of memory. AMAP retries with
more, scaling by attempt, so a single occurrence in the log that later succeeds
is normal. Persistent failures mean raising the ceiling:

```bash
--max_memory 256.GB
```

On a cluster, that is the whole answer. On a macOS or Windows laptop it is the
wrong one: Docker Desktop runs containers inside a VM with its own memory
ceiling, commonly 8 GB, and that ceiling is what killed the process. Nextflow
sized the request against the host and never saw it. Raising `--max_memory`
only makes the request larger than a limit it was already over. Raise the VM in
Docker Desktop under Settings, Resources, Memory, then keep `--max_memory` at
or below it.

**eggNOG runs out of memory.** It is the hungriest step, requesting 50 GB, and
it runs once over all samples together. On a machine that cannot reach that,
`--skip_functional_annotation` drops it along with the pangenome and the
duplication analysis and still gives you the corrected annotation. On a cluster,
use the `highmem` profile, or point it at a large-memory partition:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile singularity,slurm,highmem \
    --tsv samples.tsv
```

## Filesystem

**"No space left on device".** Usually the work directory, which grows to tens
of gigabytes. Move it somewhere with room:

```bash
-work-dir /scratch/$USER/amap_work
```

Also check your `$TMPDIR` quota, since Diamond and Bakta both use it heavily.

**Everything is slow on a shared filesystem.** Confirm `--scratch` is on. It
defaults to true in the HPC profiles, but a `standard` or `docker` profile on a
cluster will not have it.

**Cache misses on resume.** NFS timestamp drift across nodes invalidates
Nextflow's cache. The `slurm` profile already sets `process.cache = 'lenient'`
for this. If you are using a different executor, add it:

```bash
-process.cache lenient
```

## The eggNOG database

**It downloads every run.** `--eggnog_db` is pointing somewhere per-run. Give it
persistent shared storage.

**It rebuilds itself unexpectedly.** A `bacteria.dmnd` under 100 MB is treated as
a truncated download and rebuilt. Usually means an earlier download was
interrupted.

## Results that look wrong

**Far fewer genes than expected.** Check the assembly. A fragmented assembly
loses genes across contig boundaries, and they turn up in `*-rvmissing.tsv`
rather than being corrected. `*-AMAP.qc.tsv` has the contig count and N50.

**"Genes lost to (Start, End, Strand, Contig) collisions" in the log.** Two of
the classification's sources listed a gene at exactly the same place, and one
entry was kept. Genes on different contigs never collide, since the contig is
part of the key.

**`*-AMAP.clean.faa` has fewer sequences than the GFF3 has genes.** Expected.
The clean proteome drops pseudogenes and internal-stop sequences. Use
`*-AMAP.faa` if you need one protein per annotated gene.

**Almost everything is case 1.** Your samples are too far from the reference for
the alignments to find anything, so there is nothing to correct against. Use a
closer reference. See [custom reference](../user-guide/custom-reference.md).

**The Mycobrowser sections are empty.** Expected with `--skip_mycobrowser`, or
when the Mycobrowser file is missing. The pipeline warns and continues.

## Still stuck

Open an issue at
[github.com/PathoGenOmics-Lab/AMAP/issues](https://github.com/PathoGenOmics-Lab/AMAP/issues)
with:

- the command you ran
- the failing task's `.command.err`
- `results/pipeline_info/software_versions.yml`
- `results/pipeline_info/params.yml`

The last two pin down the exact versions and parameters, which is usually what
the answer turns on.
