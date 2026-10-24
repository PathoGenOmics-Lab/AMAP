# Troubleshooting

Each problem below is a heading, with the message the pipeline or the system
prints in quotes where there is one, so you can search for what you see.

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

<!-- figure: failure -->

## Before the pipeline starts

### "Unexpected input: '('" in `nextflow.config`

Nextflow 25.10 and later read pipelines with a stricter parser, which older
copies of AMAP do not pass. Nextflow keeps the copy it downloaded first, so
fetch the current one with `nextflow pull PathoGenOmics-Lab/AMAP` and run
again. To run an older AMAP release as it is, `export NXF_SYNTAX_PARSER=v1`
first.

### "Input TSV file not found"

The path in `--tsv` is wrong, or relative to somewhere other than where you
launched Nextflow.

### "Invalid TSV format"

The header must be exactly `sampleId` and `fasta`, separated by a real tab.
Spaces are the usual culprit. `cat -t samples.tsv` shows tabs as `^I`, and
works on both Linux and macOS.

### "FASTA file not found"

A path in the sample sheet does not resolve. Absolute paths avoid the
ambiguity.

### "not found, neither on the host nor inside the container"

A reference file the run reads is in neither place: a typo in its path, or one
of the defaults, which ship in the container, on a run without one. Give the
path to your own copy with its
[reference argument](../user-guide/parameters.md#reference-and-database-arguments).

### Validation failures

A check that fails stops the run, and Nextflow prints the validator's whole
report, which check failed and for which sample, under `Command output`,
followed by the task's `Work dir`. Only a check that passes writes
`results/pipeline_info/validation.log`; after a failure the log is in that work
folder. With `--strict_validation`, a warning fails the same way. Thresholds
are in [preparing inputs](../getting-started/inputs.md).

## Container problems

### "Cannot connect to the Docker daemon"

Docker is installed but not running: start Docker Desktop, or the Docker
service on Linux. `docker info` answers once it is up, where `docker --version`
answers either way.

### "permission denied while trying to connect to the Docker daemon socket"

Your user is not in the `docker` group. Either add it, or use Singularity.

### The image will not pull

Check network access from the compute nodes, which often differs from the login
node. Pull the image once on the login node, into a folder every node sees, and
have the runs look there. The `singularity` profile keeps its images in
`<outdir>/singularity_cache`, and that setting wins over
`NXF_SINGULARITY_CACHEDIR`, so the variable alone changes nothing; set the
folder in a config passed with `-c` instead:

```groovy title="shared.config"
singularity.cacheDir = '/shared/singularity_cache'
```

```bash
singularity pull /shared/singularity_cache/paururo-amap-latest.img \
    docker://paururo/amap:latest
nextflow run PathoGenOmics-Lab/AMAP -profile singularity,slurm \
    -c shared.config --tsv samples.tsv
```

`paururo-amap-latest.img` is the name Nextflow gives that image in its cache,
so the run finds it there and pulls nothing.

### Slow on Apple Silicon

The image is amd64 only, so it emulates. Fine for the test profile, painful
beyond that. Use an amd64 machine.

## Memory

### Killed with exit 137, 140 or 143

Exit 137 is the out-of-memory killer, 140 a time limit and 143 a SIGTERM, the
way a scheduler usually ends a job at its limits. AMAP retries a task that
exits with any of them, asking for more memory on each attempt, so a single
occurrence in the log that later succeeds is normal. Persistent failures mean raising the ceiling, with
`--max_memory 256.GB` for instance.

On a cluster, that is the whole answer. On a macOS or Windows laptop it is the
wrong one: Docker Desktop runs containers inside a VM with its own memory
ceiling, commonly 8 GB, and that ceiling is what killed the process. Nextflow
sized the request against the host and never saw it. Raising `--max_memory`
only makes the request larger than a limit it was already over. Raise the VM in
Docker Desktop under Settings, Resources, Memory, then keep `--max_memory` at
or below it.

### eggNOG runs out of memory

It is the hungriest step, requesting 50 GB, and it runs once over all samples
together. On a machine that cannot reach that, `--skip_functional_annotation`
drops it along with the pangenome, the duplication analysis and
`final_report.html`, and still gives you the corrected genes and
`amap_report.html`. On a cluster, use the `highmem` profile, or point it at a
large-memory partition:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile singularity,slurm,highmem \
    --tsv samples.tsv
```

## Filesystem

### "No space left on device"

Usually the work directory, which grows to tens of gigabytes on a large run.
Move it somewhere with room by adding `-work-dir /scratch/$USER/amap_work` to
the command. Also check your `$TMPDIR` quota, since Diamond and Bakta both use
it heavily, and the folder `--eggnog_db` names: `eggnog.db` alone is 41.4 GB.

### Everything is slow on a shared filesystem

Confirm `--scratch` is on. It defaults to true in the HPC profiles, but a
`standard` or `docker` profile on a cluster will not have it.

### Cache misses on resume

NFS timestamp drift across nodes invalidates Nextflow's cache. The `slurm`
profile already sets `process.cache = 'lenient'` for this. With another
executor, add `-process.cache lenient` to the command.

## The eggNOG database

### It downloads every run

`--eggnog_db` is pointing somewhere per-run. Give it persistent shared storage.

### It rebuilds itself unexpectedly

A `bacteria.dmnd` under 100 MB is treated as a truncated download and rebuilt.
Usually means an earlier download was interrupted.

### "has eggnog_db.ok, so EGGNOG_DB will not run again"

`eggnog_db.ok` tells Nextflow the database is built, so it never runs
`EGGNOG_DB` again, but what it needs is not all there. Delete the files the
message names and the next run builds it again.

## Results that look wrong

### Far fewer genes than expected

Check the assembly first: `*-AMAP.qc.tsv` has the contig count and the N50. A
gene that a contig end cuts is not lost: it comes back in pieces, each under the
gene's name and marked partial at its open end: `partial` in the master table,
the Partial state in `amap_report.html`. Only a piece under 8 residues is too
short to place: in CDC1551 cut into 40 contigs, 22 of the 25 smaller pieces
were placed, and the other three were five or six codons long
([pieces of cut genes](../how-it-works/correction-algorithm.md#pieces-of-cut-genes)).

`*-rvmissing.tsv` is a different list: the reference genes Bakta missed and AMAP
put back.

### "genes lost to (Start, End, Strand, Contig) collisions"

Two of the classification's sources listed a gene at exactly the same place,
and one entry was kept. Genes on different contigs never collide, since the
contig is part of the key.

### A sample is missing from the results ("Ignored : 1")

A task run per sample that fails for any reason but its limits (exit 137, 140
or 143) is tried twice more and then set aside, so one bad genome does not stop
the others. The run still ends `✓ COMPLETED`, its summary counts the task
under `Ignored`, and a warning after it says how many failed. The genome's
results stop at that step: set aside at the correction or before it, it has no
folder under `AMAP/`, and if it was the correction itself, `amap_report.html`
names it at the top. To find which task it was and why:

```bash
grep -w FAILED results/pipeline_info/execution_trace_*.txt
```

The trace's second column is the task's hash, the start of its folder under
`work/`, and that folder's `.command.err` says why. A dot in a sample ID is one
cause: see [preparing inputs](../getting-started/inputs.md#the-sample-sheet).

### The clean proteome has fewer sequences than the GFF3 has genes

Expected. `*-AMAP.clean.faa` drops pseudogenes and internal-stop sequences. Use
`*-AMAP.faa` if you need one protein per annotated gene.

### Almost everything is case 1

Your samples are too far from the reference for the alignments to find
anything, so there is nothing to correct against. Use a closer reference. See
[custom reference](../user-guide/custom-reference.md).

### The Mycobrowser sections are empty

Expected with `--skip_mycobrowser`, or when the Mycobrowser file is missing.
The pipeline warns and continues.

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
