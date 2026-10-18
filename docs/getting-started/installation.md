# Installation

AMAP needs two things on your machine: **Nextflow**, which runs the pipeline,
and a **container engine**, which supplies every tool the pipeline calls. You do
not install Bakta, Miniprot, eggNOG-mapper or any of the rest yourself; they all
live inside the AMAP image.

You also do not need to clone this repository. Nextflow pulls the pipeline
straight from GitHub.

## 1. Nextflow

Nextflow needs Java 11 or newer. The simplest route is Conda:

```bash
conda create -n nextflow-env
conda activate nextflow-env
conda install -c bioconda nextflow
```

If you do not have Conda, [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
is a lightweight installer. Any other Nextflow installation works equally well,
see the [Nextflow documentation](https://www.nextflow.io/docs/latest/getstarted.html).

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
docker --version        # or: singularity --version
```

## 4. Run the test profile

Before pointing AMAP at real data, confirm the whole chain works. The test
profile annotates the reference genome bundled in the container, so it needs
no data of your own and nothing to download beyond the image:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker,test
```

The first run downloads the AMAP image, which is several gigabytes because it
bundles the Bakta database and the H37Rv reference files. Expect that to take a
while. Later runs reuse the cached image.

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

Possible but not recommended: you take on installing Bakta, Miniprot,
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
