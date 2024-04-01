#!/usr/bin/env nextflow
process BAKTA_ANNOTATOR {
    publishDir "results/bakta_annotation", mode: 'copy'
    input:
    tuple val(sampleId), path(tsvFile)
    output:
    tuple val(sampleId), path("*")
    script:
    """
    read -r fastafile prefix < <(awk -F '\\t' 'NR == 1 {print \$1, \$2}' $tsvFile)
    echo "\$prefix\t\$fastafile"

    bakta --db /data/bakta_db/db-light --skip-trna --output "\$prefix" --prefix "\$prefix" --locus-tag "\$prefix" --proteins /reference/sequence.gb "\$fastafile"
    """
}