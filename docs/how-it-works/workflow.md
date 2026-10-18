# The workflow

The pipeline as a metro map: stations are the processes of `main.nf`, lines
are what travels between them. Pick a station for what it does and where its
files go, or switch an option to see what a run leaves out.

<!-- figure: metro -->

## 1. eggNOG database setup

Checks whether the eggNOG database exists at `--eggnog_db` and builds it if not.
A `bacteria.dmnd` under 100 MB is treated as a truncated download and rebuilt.

The database is large, so point `--eggnog_db` at persistent shared storage on a
cluster. Skipped entirely under `--skip_functional_annotation`.

## 2. Lineage typing

[pathotypr](https://github.com/PathoGenOmics-Lab/pathotypr) assigns each genome
to an MTBC lineage using phylogenetic marker SNPs, and the genomes are renamed
after their sample ID with the lineage in front (`L4_G26689` for sample
`G26689`). That prefix carries through every downstream file and is what lets
the report group results by lineage.

pathotypr types each contig of a draft on its own, and the lineage of the draft
is the one the markers of all its contigs add up to. A genome of one sequence
has that sequence renamed after the sample; the contigs of a draft keep their
names. A genome with no lineage markers at all is named `NA_G26689`.

MTBC-specific. Turn it off with `--skip_typing` for anything else, and sample
IDs stay as you wrote them.

## 3. Annotation with Bakta

[Bakta](https://github.com/oschwengers/bakta) produces the baseline annotation:
CDS, tRNA, rRNA and other features. This is the annotation AMAP will correct,
and it is published unmodified under `BAKTA/` so you can compare.

`--genus` and `--species` are passed through, so setting them correctly improves
the starting point.

## 4. Protein alignment with Miniprot

[Miniprot](https://github.com/lh3/miniprot) aligns the reference proteome onto
each genome, spliced and frameshift-aware. Two alignments are produced: the
reference proteome, and the false-positive set that will filter novel-gene
candidates.

Output is PAF. Coordinates are converted from 0-based to 1-based on read, and
for a protein hitting several positions, the best hit per position is kept so
multi-copy families such as PE/PPE keep all their placements.

In a draft assembly the genes a contig end cuts get a second pass. Miniprot
places the larger piece of such a gene only, so what that placement lacks is
aligned again as a protein of its own, and its placement at the end of another
contig is the rest of the gene
([pieces of cut genes](correction-algorithm.md#pieces-of-cut-genes)).

## 5. Case classification

Every Bakta CDS is compared against the alignments and sorted into cases 1 to 4,
plus the reference genes Bakta missed. Comparisons are strand-aware throughout:
a gene at `(100, 200, +)` and one at `(100, 200, -)` are different genes and are
never merged.

Tables land in `CDS_CASES/`, one per case, plus interactive overlap plots.
[The correction algorithm](correction-algorithm.md) covers the classification.

## 6. Correction

The cases are resolved into one corrected coordinate set: conflicts adjudicated,
novel genes filtered, missing reference genes rescued, duplicates removed.

Then `check_cds.py` translates every surviving CDS. A gene that does not end in
a stop codon has its coordinates extended to the next in-frame stop, up to 300
bp. Genes with internal stops are flagged as pseudogenes rather than discarded.
Anything that cannot be translated at all is dropped, so the annotation table
and the FASTA files always agree.

Outputs are the corrected gene table and its protein and nucleotide files
under `AMAP/` (`*.AMAP.tsv`, `.faa`, `.fna`). The GFF3 and GenBank files follow
once eggNOG has run, in step 7.

## 7. Functional annotation

[eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) assigns COG
categories, KEGG orthologs and pathways, GO terms, EC numbers, PFAM domains and
CAZy families.

This runs **once over all samples together**, not once per sample. With closely
related genomes most proteins are identical across samples, so they are
deduplicated first, annotated once, and the results propagated back. For a few
hundred MTBC genomes that cuts the Diamond search by an order of magnitude and
avoids hundreds of jobs contending for the same database on a shared
filesystem.

The results are merged into `*-AMAPmastertable.tsv` and into the GFF3
attributes.

## 8. Duplication analysis

Each proteome is searched against itself with Diamond to find within-genome
duplications, which are grouped into families by connected components. Then the
proteins of every sample's families are pooled and searched all against all to
build homologous communities across genomes, with network metrics computed per
community.

[Duplication metrics](duplications.md) explains the numbers.

## 9. Pangenome

[pangene](https://github.com/lh3/pangene) builds a pangenome graph from the
corrected annotations and produces a gene presence/absence matrix, which the
report turns into core, softcore, shell and cloud classifications, a frequency
spectrum and presence/absence heatmaps.

Each sample's corrected proteins are aligned with miniprot onto every other
sample's genome: the ones PATHOTYPR renamed or, with `--skip_typing`, the
genomes of the sample sheet, each under its sample ID.

## 10. Report

A single self-contained HTML file built from every output above. Sample quality,
missing reference genes by lineage, PE/PPE families, KEGG pathways, duplication
communities and topology, the pangenome, COG and Mycobrowser distributions, case
statistics with the associated tests, and a searchable table of every gene.

Alongside it, `pipeline_info/` collects the Nextflow execution report, timeline,
DAG and trace, plus `software_versions.yml` and `params.yml` for your methods
section.

## Next steps

- [The correction algorithm](correction-algorithm.md): the rules behind step 6
- [Outputs](../user-guide/outputs.md): every file this produces
