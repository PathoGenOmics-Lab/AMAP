# Preparing inputs

AMAP takes assembled genomes. It does not do assembly, so reads are not an input
here.

## The sample sheet

A tab-separated file with a header row and exactly two columns:

```tsv title="samples.tsv"
sampleId	fasta
G26689	/data/genomes/G26689.fasta
G26690	/data/genomes/G26690.fasta
G12345	./assemblies/G12345.fa
```

<!-- figure: sample-sheet -->

- **`sampleId`**: The name carried through every output directory and file. Must be unique.
  Lineage typing puts the lineage in front, so `G26689` comes out as
  `L4_G26689`.
- **`fasta`**: Path to the assembled genome. Absolute, or relative to the directory you
  launch Nextflow from.

The separator has to be an actual tab. Spaces will not do, and this is the most
common reason a sample sheet is rejected.

> [!WARNING]
> **Avoid dots in sample IDs**
>
> The classification cuts names at the first dot, so for `G26689.1` it writes
> its tables as `L4_G26689` while the pipeline waits for `L4_G26689.1`. The
> step fails, is tried twice more and is set aside, and the run finishes
> without that sample. Validation does not check this: stick to letters,
> digits and underscores.

## Pre-flight validation

`VALIDATE_INPUTS` runs before any computation, on every run except
`-profile test`, which has no sample sheet. To run only the checks and stop:

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv --validate_only
```

A check that passes writes its results to
`results/pipeline_info/validation.log`; one that fails stops the run and prints
them. Two mistakes in the sheet itself stop it even earlier, before any task
starts, with a message of the pipeline's own and no log at all: the drawing of
the sheet above shows both.

<!-- figure: validation-ranges -->

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
| Bakta database | not found, or no `version.json` or `bakta.db` in it |
| eggNOG database | `eggnog_db.ok` in it while the database is incomplete, since Nextflow then never builds it again |
| Reference files | `--ref_annotation`, `--ref_genome` or a [reference argument](../user-guide/parameters.md#reference-and-database-arguments) that the run reads but that is neither on the host nor inside the container |

### What only warns

| Check | Threshold |
| :--- | :--- |
| Genome size | below 500 kb, or above 10 Mb |
| N content | above 1% |
| GC content | below 25% or above 75% |
| Assembly fragmentation | more than 100 contigs, and again above 500 |
| Short contigs | any contig under 200 bp in a multi-contig assembly |
| Non-standard characters | anything outside `ATGCNRYSWKMBDHV` |
| eggNOG database | a `bacteria.dmnd` under 100 MB, which `EGGNOG_DB` takes for a truncated download and builds again |
| Mycobrowser table | `--mycobrowser_ref` neither on the host nor inside the container, so the master table goes without Mycobrowser categories |

Warnings do not stop a run. Add `--strict_validation` to have them stop it too.

An eggNOG database that does not exist yet is neither: the first run that
annotates builds it at `--eggnog_db`.

The thresholds are deliberately wide, since AMAP is not restricted to
tuberculosis. A 4.4 Mb genome of the *Mycobacterium tuberculosis* complex (MTBC) at 65% GC sits comfortably inside them, and
so does most of the bacterial world.

## Assembly quality

The correction compares your annotation against a reference proteome, so
assembly quality propagates directly into the result.

**Contiguity matters.** A gene that a contig end cuts comes back in pieces,
each under the gene's name and marked partial at its open end, the way NCBI
marks it. A piece under 8 residues is too short to place: in CDC1551 cut into
40 contigs, 22 of the 25 smaller pieces were placed, and the other three were
five or six codons long. The more fragmented the assembly, the more partial
genes; [pieces of cut genes](../how-it-works/correction-algorithm.md#pieces-of-cut-genes)
has the rule.

**Ns break coding sequences.** A run of Ns inside a gene stops it translating,
so it is dropped or flagged as a pseudogene. The 1% warning is worth taking
seriously.

> [!NOTE]
> **Contig names are kept**
>
> AMAP runs Bakta with `--keep-contig-headers`, so the sequence names you
> supplied are the ones in `BAKTA/` and in the corrected annotation, and the
> corrected tables say which contig each gene is on. Without that option
> Bakta relabels everything `contig_1`, `contig_2` and so on, and
> `classify_cases.py` stops with an error on a multi-contig assembly rather
> than compare genes across contigs.

> [!NOTE]
> **Closed genomes are read as circular from their header**
>
> Bakta takes a sequence for a complete, circular replicon when the
> description after its name says so: `complete`, as NCBI's "complete genome"
> and "complete sequence" do, or `circular=true`, as Unicycler's does. It then
> joins the two ends, and a gene across the origin comes out whole. Lineage
> typing renames a genome of one sequence after its sample and keeps that
> description. A closed genome whose header says neither is annotated as
> linear, and a gene across its origin is cut in two there; add
> `circular=true` after its name to have it read as circular.

> [!NOTE]
> **Draft assemblies**
>
> Every contig of a draft assembly is numbered from 1, so the same
> coordinates turn up on several contigs. The correction keys every gene on
> its contig as well as its start, end and strand, only ever compares genes
> on the same contig, and translates each one from its own contig. Cut into
> 40 contigs, CDC1551 gets every gene that lies within one contig exactly as
> the closed genome does.

## Reference data

For MTBC you need nothing else. The container bundles the H37Rv proteome, the
Bakta light database, the Mycobrowser table and the pathotypr lineage markers.
See [containers](../deployment/containers.md) for the full list.

For any other organism, supply your own reference with `--ref_annotation` and
`--ref_genome`. See [custom reference](../user-guide/custom-reference.md).

## Next steps

- [Quickstart](quickstart.md): run it
- [Parameters](../user-guide/parameters.md): the full flag reference
