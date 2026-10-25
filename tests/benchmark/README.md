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
redo MINIPROT's second pass, over what contig ends cut off the placements
(it runs AMAP code, `bin/contig_end_pieces.py`), and `classify_cases.py`,
`check_cds.py` and `deduplicate_rvs.py`, with the arguments
`modules/local/miniprot.nf` and `modules/local/cases.nf` give them, and score
the result: a few minutes for all seven genomes. `THREADS` (default 8) sets the threads for
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
  rest are marked `changed`, `deleted` when the stretch they land on is under
  half their length, or `unmapped` when an end does not land at all, and are
  not scored as exact. CP110619 also gets the ten curated H37Rv entries written
  in its own coordinates. Each gene carried across also records how many of its
  bases the alignment pairs, `aligned`: its size in the target, whatever is
  inserted into it. Truth tables made before that column are made again.
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
  alignment puts them, on their strand, covering at least half of the gene
  there and labelled with their name; how many of those AMAP flags as
  pseudogenes is given too. Half of the gene's `aligned` bases, not of the
  stretch: whatever is inserted into a gene, an IS6110, a prophage or a second
  copy of a gene the sample has twice, makes the stretch longer than the gene,
  and half of it can be more than the whole gene. K-12's b1004 lands on 63 kb
  of Sakai, a prophage inside, 597 bp of them its own. Genes marked `deleted`
  are not counted: most of them is gone (PPE50 lands on 87 of its 399 bp in
  CDC1551), and there is nothing to annotate
- **ORFs twice**: pairs of genes on one ORF (same contig and strand, the same
  reading frame, overlapping by at least half of the shorter). The references
  have none, except K-12, which nests eleven pairs of real genes this way.
- for the draft, how many proteins are the translation of their gene on the
  contig the table puts it on, and for the genes the cuts run through, how many
  of their pieces are annotated under the gene's name and marked partial
- **Bakta alone**: exact and changed genes found again, for Bakta's own calls
  of the same genome (the TSV it writes, kept with the inputs) instead of
  AMAP's table. Bakta names a call after the reference where it annotated it
  from the reference proteins AMAP gives it with `--proteins`: its DbXrefs then
  hold `UserProtein:Rv0001`. That is the name scored, so the two are measured
  against the same truth by the same rules

A reference can give one gene several proteins from alternative start codons,
all ending at its stop; K-12 does for cheA, infB, clpB and five more. AMAP
annotates such a gene once, so its proteins are scored as one gene, found if
any of them is.

`scores.tsv` in each run directory has the scores as a table. With a run to
compare with, each genome also lists the genes of either run missing from the
other, and whether a gene of the other run sits on the same reading frame there
(the same gene with other boundaries) or nothing does (a gene lost or gained).

## Expected results

With the changes up to `changed-genes` and `contig-end-pieces` (October 2026):

| Genome | exact | changed genes found | genes in the table | ORFs twice | Bakta alone: exact | Bakta alone: changed genes found |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| H37Rv against itself | 4026 / 4026 | | 4031 | 0 | 2905 | |
| H37Rv CP110619 | 3986 / 3986 | 50 / 50 | 4036 | 0 | 2880 | 29 |
| CDC1551 | 3097 / 3097 | 860 / 862 | 4014 | 2 | 2245 | 705 |
| M. bovis | 2391 / 2391 | 1515 / 1515 | 3947 | 3 | 1742 | 1290 |
| CDC1551, 40 contigs | 3072 / 3072 | | 4050 | 2 | 2231 | |
| E. coli K-12 against itself | 4286 / 4286 | | 4313 | 2 | 3841 | |
| E. coli Sakai | 220 / 220 | 3363 / 3367 | 5334 | 11 | 148 | 3202 |

Bakta alone, scored by the same rules, finds 67 to 73% of the unchanged genes
of each genome at exactly their coordinates and under their name (K-12, 90%),
and AMAP all of them. What Bakta misses it misses by its coordinates: of the
852 unchanged genes of CDC1551 it does not find, 836 have no call at their
coordinates and 16 have one under another name.

CDC1551 has one gene `deleted`, PPE50 (Rv3135). The six changed genes not
found:

- **CDC1551 Rv3018c** (PPE46): an IS6110 splits it into a 5' piece of 306
  residues, identical to H37Rv's, and a 3' piece of 130, and neither is
  annotated under its name. miniprot reports only the placements within 5% of
  a protein's best, and the best of Rv3018c is a full-length paralog, 81%
  identical; the 5' piece shows at `--outs=0.5`.
- **CDC1551 Rv3467**: in CDC1551 the phiRv1 prophage sits in this REP13E12
  copy, not in the one it occupies in H37Rv. The alignment puts Rv3467 on the
  genes at the prophage's left end, 99% identical to Rv3466 and Rv1587c and
  named so, and the copy left whole at H37Rv's phiRv1 site carries Rv3467 among
  its names.
- **Sakai b4518**: the copy at its locus lost one of the protein's two
  22-residue repeats. miniprot does not place it and the search for small
  proteins only allows substitutions; the name goes to a copy elsewhere, 79%
  identical.
- **Sakai b4726, b4671 and b4791**: micro-proteins whose Sakai copy no longer
  encodes them: a stop two codons early after one more change (b4726), a
  frameshift at residue 19 of 33 (b4671), a start codon turned into ATA
  (b4791).

Of the 50 pieces of the 25 genes the draft's cuts run through, 47 are
annotated under the gene's name and marked partial at the cut end. miniprot
places the larger piece of each; the smaller one comes from MINIPROT's second
pass, over what that placement lacks. The other three pieces are five or six
codons long (Rv2643, Rv2881c, Rv3903c), too short to look for.

In the draft every one of the 4050 proteins is the translation of its gene on
its own contig. The first run took 15 minutes with 12 threads on an M-series
Mac, Bakta included; a run on the cache takes under two minutes.

The H37Rv self-test, `tests/selftest/h37rv_selftest.sh`, scores the same genome
against every CDS of the GenBank record instead of the proteome, and expects
4031 of 4031.
