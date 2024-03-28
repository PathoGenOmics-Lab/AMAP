#!/usr/bin/env nextflow
process MINIPROT {
    publishDir "results/miniprot_reference", mode: 'copy'
    input:
    tuple val(sampleId), path(tsvFile), path(annotator)
    output:
    tuple val(sampleId), path("${sampleId}.paf"), path("${sampleId}-badbakta.paf"), path("${annotator}-H37Rv.paf")
    script:
    """
    read -r fastafile prefix < <(awk -F '\\t' 'NR == 1 {print \$1, \$2}' $tsvFile)

    miniprot --outs=0.95 -Iut16 -S "\$fastafile" /reference/aa_H37Rv.translated.fas > ${sampleId}.paf
    miniprot --outs=0.95 -Iut16 -S "\$fastafile" /reference/bad.bakta.translated.fas > ${sampleId}-badbakta.paf
    sed -i '/^[^>]/ s/\$/*/g' ${annotator}/${annotator}.faa
    miniprot --outs=0.95 -Iut16 -S /reference/new_refH37Rv.fasta  ${annotator}/${annotator}.faa > ${annotator}-H37Rv.paf
    """
}