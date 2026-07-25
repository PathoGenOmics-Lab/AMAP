---
title: How it works
---

# How it works

AMAP does not annotate genomes from scratch. It takes an annotation that already
exists, compares it against one that is trusted, and adjudicates the
differences.

The premise is that de novo annotators are good at finding genes and less good
at their exact boundaries. Bakta will call a start codon a few codons late, or
split one gene into two, or invent a short ORF that is not really there. If a
curated reference for a close relative exists, those errors are visible: align
the reference proteins onto the new genome, and where the two disagree, decide
which is more likely right.

On MTBC genomes that adjudication changes roughly 30% of genes.

## The three inputs to a decision

**The Bakta annotation.** Every CDS Bakta called on your genome, with
coordinates and a locus tag.

**The Miniprot alignments.** Every reference protein aligned onto your genome,
with the coordinates where it landed. The reference proteins carry their own
coordinates in their headers (`Rv0001_1-1524`), which lets a later step tell a
real locus from a paralogous placement.

**The false-positive set.** Proteins that a de novo annotator calls on the
reference genome itself but that are not in the curated annotation. Bakta's
characteristic mistakes, in other words. For MTBC this is a curated file; with a
[custom reference](../user-guide/custom-reference.md) AMAP derives it by
annotating the reference and taking the difference.

## The pipeline

Ten steps, from genomes to report:

1. **eggNOG database setup**, skipped if it is already built
2. **Lineage typing** with pathotypr, MTBC only
3. **Annotation** with Bakta
4. **Protein alignment** with Miniprot
5. **Case classification**: every CDS sorted by how it relates to the reference
6. **Correction**: conflicts resolved, missing genes rescued, stop codons checked
7. **Functional annotation** with eggNOG-mapper
8. **Duplication analysis**: within and across genomes
9. **Pangenome** with pangene
10. **Report**

[The workflow page](workflow.md) walks through each with its inputs and outputs.

## The correction itself

Steps 5 and 6 are what makes AMAP more than a wrapper. Every CDS is sorted into
one of five outcomes:

| Case | Meaning |
| :--- | :--- |
| Case 1 | Bakta called it, the reference has nothing there. Potentially novel. |
| Case 2 | Coordinates match the reference exactly. Kept. |
| Case 3 | Coordinates match but alignment length differs. Kept. |
| Case 4 | Overlaps a reference gene, but the coordinates conflict. Must be resolved. |
| Missing | The reference has a gene there, Bakta missed it. Rescued. |

Cases 1 and 4 are where the work is: deciding which novel genes are real, and
which set of coordinates wins a conflict.
[The correction algorithm](correction-algorithm.md) gives the rules and the
thresholds.

## In this section

- [The workflow](workflow.md): the ten steps in detail
- [The correction algorithm](correction-algorithm.md): case classification, conflict resolution, quality control
- [Duplication metrics](duplications.md): what the network numbers mean
