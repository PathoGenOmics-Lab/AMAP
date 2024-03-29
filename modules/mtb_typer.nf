#!/usr/bin/env nextflow
process MTB_TYPER {
    publishDir "${params.output}/mtb_typer", mode: 'copy'
    input:
    tuple val(sampleId), val(fasta)
    output:
    tuple val(sampleId), path("${sampleId}-mtb_typer.txt")
    script:
    """
    echo "Procesando muestra: $sampleId con Fasta: $fasta"
    python3 ${workflow.launchDir}/bin/mtb_typer.py -i $fasta -t "/reference/snps_lineages.tsv" -r "/reference/MTB_ancestor_reference.fasta"
    """
}