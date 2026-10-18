# Accuracy benchmark

How well the correction algorithm annotates seven genomes whose right answer is
known, and what changes between two versions of it. Unlike the unit tests,
this runs the real inputs the algorithm sees in the pipeline: Bakta and miniprot
on whole genomes.

## Running it

Inside the AMAP container, from the repository root:

```bash
docker run --rm --platform linux/amd64 -v "$PWD:/amap" -w /amap \
    paururo/amap:latest bash tests/benchmark/benchmark.sh
```

It takes three optional arguments: the working directory (default
`benchmark_work/`, ignored by git), a name for the run (default `latest`) and
the name of an earlier run to compare with. A typical cycle while changing the
algorithm:

```bash
# once, on the version to compare against
docker run --rm --platform linux/amd64 -v "$PWD:/amap" -w /amap \
    paururo/amap:latest bash tests/benchmark/benchmark.sh benchmark_work before

# after each change
docker run --rm --platform linux/amd64 -v "$PWD:/amap" -w /amap \
    paururo/amap:latest bash tests/benchmark/benchmark.sh benchmark_work after before
```

The first run downloads four genomes from NCBI, prepares the E. coli
reference and runs Bakta eight times: a quarter of an hour with 12 threads,
longer with fewer. All of that is cached in the working directory, so later runs only
redo `classify_cases.py`, `check_cds.py` and `deduplicate_rvs.py`, with the
arguments `modules/local/cases.nf` gives them, and score the result: a few
minutes for all seven genomes. `THREADS` (default 8) sets the threads for
Bakta, miniprot and minimap2.

## The genomes

| Genome | Assembly | Reference |
| :--- | :--- | :--- |
| H37Rv against itself | AL123456.3, in the container | H37Rv |
| H37Rv, another assembly | CP110619, in the container | H37Rv |
| M. tuberculosis CDC1551 | GCF_000008585.1 | H37Rv |
| M. bovis AF2122/97 | GCF_000195835.3 | H37Rv |
| CDC1551 cut into 40 contigs | from GCF_000008585.1 | H37Rv |
| E. coli K-12 MG1655 against itself | GCF_000005845.2 | K-12, as a custom reference |
| E. coli O157:H7 Sakai | GCF_000008865.2 | K-12, as a custom reference |

The draft is cut at random points with a fixed seed, so every run cuts it the
same way: contigs of very different sizes, about half of them
reverse-complemented, the first named after the sample and the rest
`contig_N`. The E. coli reference goes through what `PREPARE_REFERENCE` does
with `--ref_annotation` and `--ref_genome`: proteins extracted from the GenBank
record, and Bakta run on the reference to find the calls it makes wrongly.

## The truth

None of it comes from AMAP.

- **H37Rv and K-12**: every protein of the reference proteome, whose header
  gives its coordinates (`Rv0001_1-1524`), belongs exactly there, on the strand
  where the genome translates back to it.
- **The other genomes**: those genes are carried across by a whole-genome
  alignment (minimap2, preset `asm5`, primary forward alignments). A gene counts
  only where its sequence is identical in both genomes, length included; the
  rest are marked `changed` or `unmapped` and are not scored. CP110619 also gets
  the ten curated H37Rv entries written in its own coordinates.
- **The draft**: the CDC1551 truth in contig coordinates. Genes that span a cut
  are left out.

So the benchmark checks genes the reference annotates and the sample carries
unchanged. It says nothing about genes that changed, which have no right
answer to compare with.

## The scores

For each genome:

- **exact**: truth genes found at exactly their coordinates and strand, labelled
  with their own name
- **wrong label**: at the right place under another name
- **not exact**: nothing at exactly that place
- **changed genes found with their name**: genes whose sequence differs from
  the reference have no exact answer, but a gene should still sit where the
  alignment puts them, on their strand, covering at least half of that stretch
  and labelled with their name; how many of those AMAP flags as pseudogenes is
  given too
- **ORFs twice**: pairs of genes on one ORF (same contig and strand, the same
  reading frame, overlapping by at least half of the shorter). The references
  have none, except K-12, which nests eleven pairs of real genes this way.
- for the draft, how many proteins are the translation of their gene on the
  contig the table puts it on

A reference can give one gene several proteins from alternative start codons,
all ending at its stop; K-12 does for cheA, infB, clpB and five more. AMAP
annotates such a gene once, so its proteins are scored as one gene, found if
any of them is.

`scores.tsv` in each run directory has the scores as a table. With a run to
compare with, each genome also lists the genes of either run missing from the
other, and whether a gene of the other run sits on the same reading frame there
(the same gene with other boundaries) or nothing does (a gene lost or gained).

## Expected results

On `main` at 9c73f5b (October 2026):

| Genome | exact | changed genes found | genes in the table | ORFs twice |
| :--- | ---: | ---: | ---: | ---: |
| H37Rv against itself | 4026 / 4026 | | 4031 | 0 |
| H37Rv CP110619 | 3986 / 3986 | 49 / 50 | 4036 | 0 |
| CDC1551 | 3097 / 3097 | 855 / 863 | 4014 | 2 |
| M. bovis | 2391 / 2391 | 1506 / 1515 | 3946 | 3 |
| CDC1551, 40 contigs | 3072 / 3072 | | 4015 | 2 |
| E. coli K-12 against itself | 4285 / 4286 | | 4310 | 2 |
| E. coli Sakai | 219 / 220 | 3337 / 3367 | 5320 | 11 |

In the draft every one of the 4015 proteins is the translation of its gene on
its own contig. The first run took 15 minutes with 12 threads on an M-series
Mac, Bakta included; a run on the cache takes under two minutes.

The H37Rv self-test, `tests/selftest/h37rv_selftest.sh`, scores the same genome
against every CDS of the GenBank record instead of the proteome, and expects
4031 of 4031.
