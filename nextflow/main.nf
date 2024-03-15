#!/usr/bin/env nextflow

nextflow.enable.dsl=2

Channel
    .fromPath( params.tsv )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple( row.sampleId, file(row.fasta)) }
    .set { sample_fasta }

process mtb_typer {
    publishDir "results/mtb_typer", mode: 'copy'
    input:
    tuple val(sampleId), val(fasta)
    output:
    path "${sampleId}-mtb_typer.txt"
    script:
    """
    echo "Procesando muestra: $sampleId con Fasta: $fasta"
    python3 ${workflow.launchDir}/mtb_typer.py -i $fasta -t "/reference/snps_lineages.tsv" -r "/reference/MTB_ancestor_reference.fasta"
    """
}

process bakta_annotator {
    input:
    path tsvFile

    script:
    """
    read -r fastafile prefix < <(awk -F '\\t' 'NR == 1 {print \$1, \$2}' $tsvFile)
    echo "\$prefix\t\$fastafile"

    bakta --db /data/bakta_db/db-light --skip-trna --output "\$prefix" --prefix "\$prefix" --locus-tag "\$prefix" --proteins /reference/sequence.gb "\$fastafile"
    """
}

workflow {
    mtb_typer(sample_fasta)
        .set { result_mtb_typer }
    columnData = bakta_annotator(result_mtb_typer)


}
