#!/usr/bin/env nextflow
process TYPE_CASES {
    publishDir "results/type_cases/summary", pattern: "*summary.tsv", mode: 'copy'
    publishDir "results/type_cases/cases", pattern: "*case*.tsv", mode: 'copy'
    publishDir "results/type_cases/cases", pattern: "*rvmissing.tsv", mode: 'copy'
    input:
    tuple val(sampleId), path(annotator), path(miniprot), path(badbakta), path(cds)
    output:
    path "*"
    script:
    """
    python3 ${workflow.launchDir}/bin/classify_cases.py -b ${annotator}/${annotator}.tsv -p ${miniprot} -bb ${badbakta} -c ${cds}
    """
}