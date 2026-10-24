# Installation

AMAP needs two things on your machine: **Nextflow**, which runs the pipeline,
and a **container engine**, which supplies every tool the pipeline calls. You do
not install Bakta, miniprot, eggNOG-mapper or any of the rest yourself; they all
live inside the AMAP image.

You also do not need to clone this repository. Nextflow pulls the pipeline
straight from GitHub.

<!-- figure: image -->

## 1. Nextflow

Nextflow needs Java 17 or newer, and Conda brings both:

```bash
conda create -n nextflow-env -c conda-forge -c bioconda nextflow
conda activate nextflow-env
```

If you do not have Conda, [Miniforge](https://github.com/conda-forge/miniforge)
is the one to install: it comes set up for conda-forge, the channel Bioconda
builds on. Any other Nextflow installation works equally well, see the
[Nextflow documentation](https://www.nextflow.io/docs/latest/getstarted.html).

AMAP requires Nextflow 24.04 or newer, which is enforced by the pipeline
manifest.

## 2. A container engine

Pick whichever your machine allows.

**Docker**

Best on a laptop or workstation where you have root. Follow the
[Docker installation guide](https://docs.docker.com/engine/install/) for
your platform.

**Singularity / Apptainer**

The usual choice on an HPC cluster, where Docker is normally not permitted.
Follow the
[Apptainer installation guide](https://apptainer.org/docs/user/latest/quick_start.html).
Many clusters already provide it as a module:

```bash
module load singularity   # or: module load apptainer
```

## 3. Check the installation

```bash
nextflow -version
docker info --format 'Docker can give containers {{.NCPU}} CPUs and {{.MemTotal}} bytes'
```

On a machine with Docker Desktop running, the second line prints something like
this:

```text
Docker can give containers 14 CPUs and 8217014272 bytes
```

Unlike `docker --version`, it fails when the Docker daemon is not running, and
it prints the memory your containers can really have. Under 8 GB? Raise it in
Docker Desktop, Settings, Resources. With Singularity or Apptainer,
`singularity --version` is the check.

## 4. Run the test profile

Before pointing AMAP at real data, confirm the whole chain works. The test
profile annotates the H37Rv genome bundled in the image, so it needs no data of
your own, and `--skip_functional_annotation` keeps it to the correction:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker,test \
    --skip_functional_annotation --outdir amap-test
```

It downloads nothing beyond the pipeline and its image, writes `amap-test/` and
Nextflow's `work/` in the folder you run it from, and ends with a `✓ COMPLETED`
block. On an Apple M4 Pro, under emulation and with the image already pulled,
it took 2 min 47 s.

The first run downloads the AMAP image, about 3.4 GB compressed and 12 GB on
disk, because it bundles the Bakta database and the H37Rv reference files.
Expect that to take a while. Later runs reuse the cached image.

### The eggNOG database

Without `--skip_functional_annotation`, the test also builds the eggNOG
database, the one thing the image does not hold. It is the largest thing AMAP
keeps: `eggnog.db` alone unpacks to 41.4 GB from a 6.8 GB download, and the
build also downloads the 9.5 GB eggNOG proteome, makes the bacteria-only Diamond
index from it and deletes it. It is built once, in the folder `--eggnog_db`
names, and every run after reuses it, so point that somewhere with room.

> [!TIP]
> **Check inputs before a long run**
>
> Once you have your own samples, `--validate_only` checks every FASTA and
> database path without executing anything:
>
> ```bash
> nextflow run PathoGenOmics-Lab/AMAP -profile docker \
>     --tsv samples.tsv --validate_only
> ```
>
> See [inputs](inputs.md) for what it looks for.

## Running without a container

Possible but not recommended: you take on installing Bakta, miniprot,
eggNOG-mapper, Diamond, pangene, pathotypr, R and the report packages yourself,
and pointing AMAP at each reference file with the
[reference arguments](../user-guide/parameters.md#reference-and-database-arguments).

The repository ships an `environment.yml` that reproduces the software stack
under Conda:

```bash
conda env create -f environment.yml
conda activate amap
```

You will still need to supply the Bakta database and the reference files, which
the container would otherwise provide. See
[containers](../deployment/containers.md) for what the image contains.

## Next steps

- [Quickstart](quickstart.md): your first real run
- [Preparing inputs](inputs.md): the sample sheet and what makes a genome pass validation
