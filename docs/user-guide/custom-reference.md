# Custom reference

AMAP corrects an annotation by refereeing it against a trusted one. For MTBC
that trusted annotation is H37Rv, bundled in the image. Point
`--ref_annotation` and `--ref_genome` at your own and the same machinery works
for any prokaryote.

```bash
nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv \
    --skip_typing --skip_mycobrowser \
    --genus Escherichia --species "Escherichia coli" \
    --ref_annotation my_reference.gbk \
    --ref_genome my_reference.fasta
```

Both must be given together. Setting one without the other stops the run
immediately.

## What counts as a good reference

A **closed, well-curated genome of a close relative**, ideally the type strain
or a reference assembly from RefSeq. The correction assumes the reference genes
are right, so its quality is the ceiling on your result.

- **`--ref_annotation`**: GFF3 (`.gff`, `.gff3`) or GenBank (`.gbk`, `.gb`, `.genbank`, `.gbff`).
  Detected from the extension, with GenBank tried first for anything else.
- **`--ref_genome`**: FASTA of the same genome. It must be the same assembly the annotation
  describes, since coordinates are taken at face value.

## What AMAP does with it

`PREPARE_REFERENCE` runs once, before any sample is processed, and does three
things.

<!-- figure: prepare-reference -->

**1. Extracts the reference proteome.** Every CDS becomes a protein, headed
`GENE_ID_START-END contig=SEQUENCE`, `SEQUENCE` being the name `--ref_genome`
gives the sequence the gene is on.

- **The header matters.** The paralog step reads its coordinates back to
  predict where each protein belongs in the sample, which tells its own locus
  from a paralog's.
- **From GenBank**, the existing `/translation` qualifier is used when present;
  otherwise the sequence is extracted and translated.
- **From GFF3**, the lines of each CDS are joined first, as GenBank's `join()`
  is, and translated with the bacterial code, reading `transl_except` (a
  selenocysteine) where it is given. A GFF3 protein under 8 amino acids is left
  out.
- **The gene ID** comes from `locus_tag`, then `gene`, then `protein_id` for
  GenBank; `locus_tag`, `ID`, `gene`, then `Name` for GFF3. A CDS with none of
  these is skipped.

> [!NOTE]
> **Details**
>
> A CDS across the origin of a circular sequence ends past the sequence's last
> base, the way GFF3 writes it: tagA on plasmid pO157 (92721 bp) of *E. coli*
> O157:H7 Sakai, `join(92527..92721,1..2502)` in its GenBank record, is
> `pO157p01_92527-95223`. On *E. coli* K-12 the two formats give the same 4318
> proteins, and on Sakai the same 5155.

**2. Annotates the reference genome with Bakta.** The same Bakta that will
annotate your samples, run on a genome whose true annotation is known, keeping
the FASTA's names for its sequences (`--keep-contig-headers`).

**3. Works out what Bakta invents.** Any Bakta CDS that does not overlap a
curated gene of its own sequence is a false positive on this genome, and by
extension a call Bakta is prone to make on your samples too. These become the
`bad_bakta` set used to filter novel-gene candidates in
[case 1](../how-it-works/correction-algorithm.md#2-filtering-case-1). Bakta
writes a call across the origin of a sequence whose FASTA description says
complete from its start to its stop, start after end (stcE on pO157 is
92563..2502), and it is compared with the genes on both sides of the origin.
The calls on a sequence the annotation has no gene on are not judged: the
annotation says nothing about them.

For MTBC this set was curated by hand. With a custom reference AMAP derives the
equivalent automatically, which is what makes organism-agnostic mode work
rather than merely run.

## What you must turn off

- **`--skip_typing`**: pathotypr assigns MTBC lineages from tuberculosis-specific marker SNPs.
  Meaningless elsewhere, and it will fail. AMAP warns when `--genus` is not
  Mycobacterium and this is not set.
- **`--skip_mycobrowser`**: Mycobrowser functional categories are an MTBC resource. Auto-enabled when
  `--ref_annotation` is set and `--genus` is not Mycobacterium, but setting it
  explicitly is clearer.

Also set `--genus` and `--species` to match your organism, since Bakta uses
them for naming and reference selection.

## What stays the same

Everything else. Case classification, conflict resolution, missing-gene rescue,
the stop codon check, eggNOG functional annotation, duplication analysis, the
pangenome and the report all work identically. The report simply omits the
lineage and Mycobrowser sections, and the pangenome calls each genome by its
sample ID where a typed run would use its lineage name.

## What to expect from a real run

*E. coli* O157:H7 Sakai corrected against K-12 MG1655, both straight from
RefSeq, in [the tutorial's run](../getting-started/tutorial.md#9-the-same-thing-for-another-organism):

```text
reference proteins extracted from the K-12 GenBank   4318
CDS Bakta called on that same reference genome       4296
Bakta false positives derived from the difference     142

case2                                                3324
case3                                                  31
case4                                                 507
case1_kept                                           1387
rv_missing                                            133
final annotation                                     5335 genes
```

The 142 are the automatically derived equivalent of the hand-curated MTBC
false-positive list, and 1387 case 1 genes is what you would expect between two
*E. coli* strains that genuinely differ.

If your case 2 comes out at or near zero, stop and look at the reference rather
than the sample. Every gene missing an exact match by a constant offset points
at the reference, not at 4000 coincidences.

## Worked example

*E. coli* K-12 MG1655 against the RefSeq reference:

```bash
# Fetch the reference
datasets download genome accession GCF_000005845.2 --include genome,gff3
unzip ncbi_dataset.zip

nextflow run PathoGenOmics-Lab/AMAP -profile docker \
    --tsv samples.tsv \
    --skip_typing --skip_mycobrowser \
    --genus Escherichia --species "Escherichia coli" \
    --ref_annotation ncbi_dataset/data/GCF_000005845.2/genomic.gff \
    --ref_genome ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna
```

## Things that trip people up

> [!WARNING]
> **The annotation and the genome must match**
>
> Coordinates from the annotation are used directly against the FASTA. An
> annotation from a different assembly of the same organism produces
> nonsense: extracted proteins will not translate, and the reference proteome
> will be mostly empty. `PREPARE_REFERENCE` writes the count of extracted
> proteins into its work directory, where
> `grep -h 'Curated proteins' work/*/*/.command.out` finds it, so check it
> looks like a genome's worth.

> [!WARNING]
> **Reference genes need identifiers**
>
> A GFF3 whose CDS features carry no `locus_tag`, `ID`, `gene` or `Name` is
> skipped entirely, and the run fails with "No proteins extracted". This
> happens with minimal GFF3 exports from some annotation tools.

> [!NOTE]
> **Distance to the reference**
>
> The correction expects your samples to be the same species as the
> reference. The further away they are, the more genes fall into case 1
> (nothing in the reference to compare against) and the less AMAP has to say.
> It is a correction tool, not an ortholog finder.

## Next steps

- [How it works](../how-it-works/index.md): what the correction does with the reference
- [Parameters](parameters.md): the full flag reference
