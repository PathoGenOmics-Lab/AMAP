# Duplication metrics

AMAP looks at gene duplication twice: within each genome, and across the whole
set.

## How the families are found

Each corrected proteome is searched against itself with Diamond, excluding
self-hits. A pair is kept when the identical residues of the alignment cover
**80% or more** of either the query or the subject, so a short domain match
between two otherwise unrelated proteins is not counted as duplication. Kept
pairs become links in a graph, and its connected pieces become the
duplication families, listed in
`DUPLICATIONS/{sample}/duplications_{sample}_groups.tsv`.

<!-- figure: dup-family -->

Because the rule takes the larger of the two shares, a fragment joins the
family of the gene it came from: a match that covers almost all of a short
protein counts, however little of the long one it covers.

## How the communities are built

The proteins of every sample's families are then pooled and searched all
against all. Within each genome, the pairs that pass the same rule make a
network, which the Louvain method splits into **communities**, with a fixed
seed so that a rerun splits it the same way. Two communities are
**homologous** when at least one hit joins them, and each homologous pair is
one row of `DUPLICATIONS/homologous_communities.tsv`, with an interactive
network in the matching HTML file.

<!-- figure: dup-communities -->

## The columns

A row describes two communities: the columns ending in `1` are the first, those
ending in `2` the second. The first is the one whose sample comes first by
name, or the one with the lower number when both are in the same sample, and
the rows come in that order too, so a rerun writes the same file. For each:

- **`Genome`, `Community`**: the sample, and the community's number in it.
- **`Size`**: how many proteins it holds.
- **`Elements`**: the eggNOG COG categories of its proteins, with how many
  have each, the largest first. A protein eggNOG gave no category counts as
  `-`, eggNOG's mark for none.
- **`AvgDegree`**: links per protein, on average. In a community of n
  proteins it reaches n - 1 when every protein is linked to every other.
- **`Density`**: the links it has as a share of the pairs it could have. 1
  means every protein matches every other.
- **`AvgClustering`**: for each protein, how often two of its partners are
  linked to each other as well, averaged over the community. 1 when every
  protein is linked to every other, 0 for a chain or a star.
- **`AvgPathLength`**: how many links separate two of its proteins, on
  average, over its largest connected piece.
- **`Nodes`**: the proteins themselves.

## Reading them together

**Density and clustering at 1** mean copies that all match each other: a
recent expansion of near-identical copies. In CDC1551 these are the
transposases of insertion sequences, the ESAT-6-like Esx proteins and PE
proteins.

**Lower density and longer paths** mean copies that have drifted apart, each
still matching some of the others but not all, as M. bovis's copies of the
transposases above have. Or, less happily, an artefact of assembly
fragmentation, where one gene broken across contigs looks like several partial
copies. Cross-checking against the contig count in `*-AMAP.qc.tsv` tells the
two apart.

> [!NOTE]
> **Duplications and assembly quality**
>
> Repetitive families are exactly what short-read assemblers collapse or
> fragment. A genome with markedly fewer duplications than its relatives may
> have collapsed them rather than lost them. The metrics describe the
> assembly you provided, not necessarily the organism.

## Next steps

- [Outputs](../user-guide/outputs.md#the-key-files): where the duplication files live
- [The workflow](workflow.md): where this step sits
