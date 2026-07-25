---
title: Preparing inputs
---

# Preparing inputs

AMAP takes assembled genomes. It does not do assembly, so reads are not an input
here.

## The sample sheet

A tab-separated file with a header row and exactly two columns:

```tsv title="samples.tsv"
sampleId	fasta
G26689	/data/genomes/G26689.fasta
G26690	/data/genomes/G26690.fasta
L4_G12345	./assemblies/L4_G12345.fa
```

`sampleId`
: The name carried through every output directory and file. Must be unique.

`fasta`
: Path to the assembled genome. Absolute, or relative to the directory you
  launch Nextflow from.

The separator has to be an actual tab. Spaces will not do, and this is the most
common reason a sample sheet is rejected.

!!! warning "Avoid dots in sample IDs"

    Several scripts derive the sample name by cutting at the first dot, so
    `G26689.1` becomes `G26689` and two such samples would collide. Stick to
    letters, digits and underscores.

## Pre-flight validation

`VALIDATE_INPUTS` runs before any computation, on every run. To run only the
checks and stop:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --validate_only
```

Results are written to `results/pipeline_info/validation.log`.

### What fails the run

| Check | Threshold |
| :--- | :--- |
| Genome size | below 100 kb, or above 20 Mb |
| N content | above 10% |
| GC content | below 20% or above 80% |
| Empty sequences | any record with no bases |
| Duplicate contig headers | any repeated identifier |
| Sample sheet | missing `sampleId` or `fasta` column, no data rows, duplicate or empty sample IDs |
| FASTA path | file does not exist or cannot be parsed |

### What only warns

| Check | Threshold |
| :--- | :--- |
| Genome size | below 500 kb, or above 10 Mb |
| N content | above 1% |
| GC content | below 25% or above 75% |
| Assembly fragmentation | more than 100 contigs, and again above 500 |
| Short contigs | any contig under 200 bp in a multi-contig assembly |
| Non-standard characters | anything outside `ATGCNRYSWKMBDHV` |

Warnings do not stop a run. Add `--strict` to the validation step if you want
them treated as failures.

The thresholds are deliberately wide, since AMAP is not restricted to
tuberculosis. A 4.4 Mb MTBC genome at 65% GC sits comfortably inside them, and
so does most of the bacterial world.

## Assembly quality

The correction compares your annotation against a reference proteome, so
assembly quality propagates directly into the result.

**Contiguity matters.** A gene broken across two contigs cannot be recovered.
The more fragmented the assembly, the more reference genes turn up as missing
rather than corrected.

**Ns break coding sequences.** A run of Ns inside a gene stops it translating,
so it is dropped or flagged as a pseudogene. The 1% warning is worth taking
seriously.

!!! warning "Multi-contig assemblies and coordinate collisions"

    The correction step keys genes on `(start, end, strand)` without the contig
    name. On a closed single-contig genome, which is the MTBC case AMAP was
    built around, that is unambiguous. On a draft assembly, two genes at the
    same coordinates on different contigs collide, and one is lost.

    `classify_cases.py` prints a warning when it detects this, and the count
    appears in the run log. In practice the collision needs an exact match of
    start, end and strand across contigs, so it is rare, but on a highly
    fragmented assembly it is worth checking the log for it.

## Reference data

For MTBC you need nothing else. The container bundles the H37Rv proteome, the
Bakta light database, the Mycobrowser table and the pathotypr lineage markers.
See [containers](../deployment/containers.md) for the full list.

For any other organism, supply your own reference with `--ref_annotation` and
`--ref_genome`. See [custom reference](../user-guide/custom-reference.md).

## Next steps

- [Quickstart](quickstart.md): run it
- [Parameters](../user-guide/parameters.md): the full flag reference
