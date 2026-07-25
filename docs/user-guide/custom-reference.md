---
title: Custom reference
---

# Using your own reference

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

`--ref_annotation`
: GFF3 (`.gff`, `.gff3`) or GenBank (`.gbk`, `.gb`, `.genbank`, `.gbff`).
  Detected from the extension, with GenBank tried first for anything else.

`--ref_genome`
: FASTA of the same genome. It must be the same assembly the annotation
  describes, since coordinates are taken at face value.

## What AMAP does with it

`PREPARE_REFERENCE` runs once, before any sample is processed, and does three
things.

**1. Extracts the reference proteome.** Every CDS becomes a protein, headed
`GENE_ID_START-END`. Those embedded coordinates are not decoration: the
reconciliation step reads them back to tell a real gene locus from a paralogous
placement, so the header format matters. From GenBank the existing
`/translation` qualifier is used when present, otherwise the sequence is
extracted and translated. From GFF3 everything is translated from the genome,
and fragments under 10 amino acids are dropped.

The gene ID comes from `locus_tag`, then `gene`, then `protein_id` for GenBank;
`locus_tag`, `ID`, `gene`, then `Name` for GFF3. A CDS with none of these is
skipped.

**2. Annotates the reference genome with Bakta.** The same Bakta that will
annotate your samples, run on a genome whose true annotation is known.

**3. Works out what Bakta invents.** Any Bakta CDS that does not overlap a
curated gene is a false positive on this genome, and by extension a call Bakta
is prone to make on your samples too. These become the `bad_bakta` set used to
filter novel-gene candidates in
[case 1](../how-it-works/correction-algorithm.md#2-filtering-case-1).

For MTBC this set was curated by hand. With a custom reference AMAP derives the
equivalent automatically, which is what makes organism-agnostic mode work
rather than merely run.

## What you must turn off

`--skip_typing`
: pathotypr assigns MTBC lineages from tuberculosis-specific marker SNPs.
  Meaningless elsewhere, and it will fail. AMAP warns when `--genus` is not
  Mycobacterium and this is not set.

`--skip_mycobrowser`
: Mycobrowser functional categories are an MTBC resource. Auto-enabled when
  `--ref_annotation` is set and `--genus` is not Mycobacterium, but setting it
  explicitly is clearer.

Also set `--genus` and `--species` to match your organism, since Bakta uses
them for naming and reference selection.

## What stays the same

Everything else. Case classification, conflict resolution, missing-gene rescue,
the stop codon check, eggNOG functional annotation, duplication analysis, the
pangenome and the report all work identically. The report simply omits the
lineage and Mycobrowser sections.

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

!!! warning "The annotation and the genome must match"

    Coordinates from the annotation are used directly against the FASTA. An
    annotation from a different assembly of the same organism produces
    nonsense: extracted proteins will not translate, and the reference proteome
    will be mostly empty. `PREPARE_REFERENCE` prints the count of extracted
    proteins, so check it looks like a genome's worth.

!!! warning "Reference genes need identifiers"

    A GFF3 whose CDS features carry no `locus_tag`, `ID`, `gene` or `Name` is
    skipped entirely, and the run fails with "No proteins extracted". This
    happens with minimal GFF3 exports from some annotation tools.

!!! note "Distance to the reference"

    The correction expects your samples to be the same species as the
    reference. The further away they are, the more genes fall into case 1
    (nothing in the reference to compare against) and the less AMAP has to say.
    It is a correction tool, not an ortholog finder.

## Next steps

- [How it works](../how-it-works/index.md): what the correction does with the reference
- [Parameters](parameters.md): the full flag reference
