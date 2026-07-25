---
title: The correction algorithm
---

# The correction algorithm

The core of AMAP. It referees between the initial Bakta annotation and the
Miniprot alignments of a trusted reference proteome, classifying every candidate
gene and then applying a set of rules to resolve the disagreements.

Coordinates throughout are 1-based inclusive, so a gene's length is
`end - start + 1`. Every comparison is strand-aware: genes on opposite strands
at the same position are different genes and never conflict.

## 1. Classification

Every Bakta CDS is compared against the Miniprot alignments.

**Case 2, perfect match.** Same start, same end, same strand as an alignment,
and the same length. Correct, kept as is.

**Case 3, imperfect match.** Same coordinates, but the aligned length differs,
typically from indels. Reliable enough, kept.

**Case 4, conflicting overlap.** A Bakta CDS and an alignment overlap on the
same strand, but their coordinates do not match. This is the conflict that has
to be adjudicated.

**Case 1, potentially novel.** A Bakta CDS with no overlap to any alignment.
Either a genuine gene absent from the reference, or something Bakta invented.

**Missing reference genes.** An alignment with no overlapping Bakta CDS. The
reference has a gene there and the annotator missed it. Alignments with no
aligned bases at all are discarded first.

Each set is written to its own table under `CDS_CASES/`, so the whole decision
is auditable.

## 2. Filtering case 1

Novel-gene candidates are scrutinised before being believed. A candidate is
**removed** if either holds:

- it overlaps the [false-positive set](../user-guide/custom-reference.md) by
  **80% or more** of its own length, meaning it is a call the annotator is known
  to make wrongly on the reference genome itself
- it aligns to the reference at **80% or more** match identity, meaning it is
  not novel at all but a shifted version of a gene that is already there

Candidates that pass both, and candidates with no reference alignment to judge
by, are kept as `case1_kept`.

## 3. Resolving case 4

### Single overlaps

One Bakta CDS against one reference alignment. Lengths are compared, and:

- **the reference wins by default**
- the Bakta version wins only if it is longer by **more than 550 bp**

The asymmetry is deliberate. A reference alignment is usually the better
boundary, but a much longer Bakta call normally means the alignment fragmented,
and discarding a full-length gene in favour of a partial alignment would be the
worse error. If the Bakta CDS overlaps more than one alignment, the reference
wins regardless.

### Multiple overlaps

Several Bakta CDS against one alignment, the signature of a gene that has been
split. The best-overlapping Bakta CDS is picked first, then compared to the
alignment on length:

- the reference wins if it is at least as long, or if the Bakta CDS beats it by
  **less than 500 bp**
- otherwise the Bakta CDS wins

## 4. Deduplication on insert

Each winner is checked against the genes already accepted. Where a new
coordinate overlaps an existing one by **70% or more** on the same strand:

- the new gene displaces the existing one if it is longer by **more than 550 bp**
- it also displaces it if the existing gene is shorter by less than that margin
  **and is not a reference gene**, reference genes being identified by the
  coordinates embedded in their identifiers
- otherwise the existing gene stands and the new coordinate is dropped

## 5. Reconciliation

Reference proteins carry their true coordinates in their headers
(`Rv0001_1-1524`). Two things follow.

**Coordinate replacement.** When an accepted gene matches a reference gene by
exact name and their coordinates overlap by more than half the shorter one, with
neither more than twice the other, the reference coordinates replace the
current ones. Name matching is exact after normalisation, so `Rv3508` does not
match `Rv3508c` and `Rv1` does not match `Rv10`.

**Paralog removal.** If a gene's own embedded coordinates overlap its actual
position by less than half, Miniprot placed it on a paralogous region. It is
removed, unless one of these holds:

- another alignment at that exact position does have matching embedded
  coordinates, proving the locus is real
- Bakta independently called a gene at exactly those coordinates
- the protein has more than two alignments across the genome, so it is
  multi-copy and placement away from the reference position is expected, which
  is what keeps PE/PPE and IS elements intact

When a gene is removed and Bakta had called something at roughly that position
with a comparable length, the Bakta call is restored in its place.

## 6. Rescuing missing genes

Reference genes still absent from the set are added back, best alignment first.
A placement is only accepted if it overlaps the protein's embedded coordinates
by at least half, so a gene that only aligns to a paralogous region is left out
rather than placed wrongly. If the best hit fails that test, the next is tried.

## 7. Final deduplication

Where one gene ends up at several coordinates, entries fully contained inside a
different gene are dropped, always keeping at least one entry per gene so
nothing disappears entirely.

## 8. Stop codon check

`check_cds.py` translates every surviving CDS.

**No terminal stop codon.** Coordinates are extended to the next in-frame stop,
searching up to 300 bp in the right direction for the strand. If the extension
introduces internal stops or finds nothing, the original coordinates are kept
and the gene is flagged as a pseudogene.

**Reference genes are exempt.** If the coordinates match the ones embedded in
the identifier, the gene is a known pseudogene in the reference and is left
alone rather than extended.

**Internal stop codons.** Flagged as a pseudogene and kept, with the flag
propagated to the master table and the GFF3. `*-AMAP.clean.faa` excludes them.

**Untranslatable genes are dropped** from the annotation table as well as the
FASTA files, so the two never disagree.

<div class="amap-figure" markdown>
![AMAP correction algorithm](../assets/correction-algorithm.svg)
</div>

## Threshold summary

| Threshold | Where | Value |
| :--- | :--- | :--- |
| False-positive overlap | case 1 filtering | 80% |
| Reference match identity | case 1 filtering | 80% |
| Bakta length advantage | single overlaps | 550 bp |
| Bakta length advantage | multiple overlaps | 500 bp |
| Overlap for deduplication | on insert | 70% |
| Length margin for displacement | on insert | 550 bp |
| Placement overlap | reconciliation and rescue | 50% |
| Multi-copy alignment count | paralog exemption | more than 2 |
| Stop codon search | quality control | 300 bp |

These are constants in `bin/classify/`, tuned on MTBC. Another organism may want
different values, particularly the 550 bp length margins, which assume genes of
roughly bacterial size.

## Next steps

- [The workflow](workflow.md): where this sits in the pipeline
- [Outputs](../user-guide/outputs.md): the case tables that record every decision
