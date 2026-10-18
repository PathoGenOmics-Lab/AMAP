# Duplication metrics

AMAP looks at gene duplication twice: within each genome, and across the whole
set.

## How the communities are built

Each corrected proteome is searched against itself with Diamond, excluding
self-hits. A pair is kept when the aligned identical residues cover **80% or
more** of either the query or the subject, so a short domain match between two
otherwise unrelated proteins is not counted as duplication. Surviving pairs
become edges in a graph, and connected components become duplication families.

All samples are then pooled and searched all-against-all. The resulting graph is
partitioned into **homologous communities** with the Louvain algorithm, and the
metrics below are computed per community. Results land in
`DUPLICATIONS/homologous_communities.tsv`, with an interactive network in the
matching HTML file.

## The metrics

Each appears in the table suffixed `1` and `2`, for the two levels of the
analysis.

- **Degree centrality** (`AverageDegree1-2`): Average number of direct connections per node. High values mean genes in the
  community are similar to many others in it, the signature of a large expanded
  family.
- **Closeness centrality** (`AverageCloseness1-2`): Average distance from a node to every other. High closeness means a compact
  community where members are all mutually similar, rather than a chain of
  progressively diverging copies.
- **Betweenness centrality** (`AverageBetweenness1-2`): How often a node sits on the shortest path between two others. A gene with
  high betweenness bridges otherwise separate subgroups, often an intermediate
  form linking two diverged subfamilies.
- **Eigenvector centrality** (`Eigen1-2`): Weights a node's connections by how well connected its neighbours are.
  Identifies the genes at the core of a family rather than at its edges.
- **Katz centrality** (`Katz1-2`): Like eigenvector centrality but counting distant nodes too, with weight
  decaying by distance. More stable on sparse graphs where eigenvector
  centrality concentrates on one dense cluster.
- **Network density** (`Density1-2`): Realised connections as a fraction of possible ones. Near 1 means every gene
  resembles every other, an unresolved family of near-identical copies. Low
  density means a sparse, diverged group.
- **Clustering coefficient** (`ClusteringCoefficient1-2`): How often two neighbours of a node are themselves connected. High values mean
  well-defined subgroups.
- **Average path length** (`AveragePathLength1-2`): Mean shortest-path distance over all pairs. Short paths mean any member is
  reachable from any other in a few steps.
- **Transitivity** (`transitivity1-2`): The clustering coefficient applied to the whole network: the overall tendency
  to form triangles.
- **Modularity** (`modularity1-2`): How cleanly the network divides into communities. High modularity means dense
  groups with few links between them, so the partition into families is
  well supported. Low modularity means the split is somewhat arbitrary.

## Reading them together

Two patterns come up repeatedly.

**Dense, high-modularity communities** are recent expansions: many near-identical
copies, clearly separated from everything else. In MTBC these are typically
PE/PPE genes and IS elements.

**Sparse, low-modularity communities** are ancient duplications that have
diverged, or, less happily, an artefact of assembly fragmentation where one gene
broken across contigs looks like several partial copies. Cross-checking against
the contig count in `*-AMAP.qc.tsv` distinguishes the two.

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
