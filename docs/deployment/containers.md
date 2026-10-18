# Containers

Every tool AMAP calls lives in one image, `paururo/amap`, published on
[Docker Hub](https://hub.docker.com/r/paururo/amap). Both the `docker` and
`singularity` profiles use it, so nothing needs installing beyond a container
engine.

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker --tsv samples.tsv
nextflow run PathoGenOmics-Lab/AMAP -profile singularity --tsv samples.tsv
```

Nextflow pulls the image on first use and caches it. Singularity converts it to
a `.sif` on the way, which takes a few minutes once.

## What is inside

**Tools**: Bakta 1.12.0, Miniprot 0.18, Diamond 2.1.24, eggNOG-mapper 2.1.13,
pangene 1.1, pathotypr 1.0.0, gfatools.

**Runtimes**: Python 3.9 with pandas, NumPy, Biopython, matplotlib, plotly,
NetworkX, intervaltree, pyvis and python-louvain. R 4.4 with the report
packages.

**Bundled data** under `/reference/`, so a standard MTBC run needs no downloads:

| File | What it is |
| :--- | :--- |
| `MTB_ancestor_reference.fasta` | MTBC ancestral reference genome |
| `aa_H37Rv.translated.fas` | H37Rv reference proteome |
| `bad.bakta.translated.fas` | Curated Bakta false positives |
| `new_refH37Rv.fasta` | H37Rv reference genome |
| `sequence.gb` | H37Rv GenBank, for Bakta |
| `Mycobacterium_tuberculosis_H37Rv_txt_v4.txt` | Mycobrowser functional categories |
| `snps_lineages.tsv` | pathotypr lineage markers |

Plus the Bakta light database at `/data/bakta_db`.

The eggNOG database is **not** bundled, because of its size. AMAP downloads and
builds it on first use and caches it at `--eggnog_db`.

## Platform

The image is `linux/amd64` only, and every Dockerfile stage pins it. On Apple
Silicon it runs under emulation, which works but is slow. For anything past a
couple of genomes, use an amd64 machine or a cluster.

## Pinning a version

`latest` tracks the newest release. For a reproducible analysis, pin the
version in `nextflow.config` or on the command line:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv \
    -process.container 'paururo/amap:1.0.1'
```

`results/pipeline_info/software_versions.yml` records the exact version of every
tool a run used, which is what belongs in a methods section.

## Building it yourself

Only needed if you are changing what is in the image.

```bash
git clone https://github.com/PathoGenOmics-Lab/AMAP
cd AMAP
docker build --platform linux/amd64 -f docker/Dockerfile -t amap:dev .
```

The build takes about an hour: it downloads the Bakta database and solves a
conda environment with R. Check the result before using it:

```bash
docker/smoke_test.sh amap:dev
```

That verifies the tools, the Python and R libraries, the bundled reference data,
and runs the repository test suite with the image interpreter. The last part
matters, since a library update can break the pipeline scripts without the image
looking any different.

Publishing is automated: the `Docker image` GitHub Actions workflow builds,
smoke tests and pushes on a `v*` tag or on demand.

> [!WARNING]
> **Reference downloads can expire**
>
> Several reference files are fetched from Google Drive during the build. Those
> links have expired before. If a build fails on a `wget` in the data stage,
> that is why, and the URL needs replacing. The Zenodo-hosted files
> (`MTB_ancestor_reference.fasta`, the pathotypr markers) are not affected.

## Singularity on a cluster

Set a shared cache so every job reuses one converted image instead of each
building its own:

```bash
export NXF_SINGULARITY_CACHEDIR=/shared/singularity_cache
```

Put it in your shell profile or the launch script. Without it, a large run can
have dozens of jobs converting the same image at once.

More cluster specifics in [deployment on HPC](hpc.md).

## Next steps

- [HPC deployment](hpc.md): scratch, shared filesystems, schedulers
- [Troubleshooting](troubleshooting.md): when a run fails
