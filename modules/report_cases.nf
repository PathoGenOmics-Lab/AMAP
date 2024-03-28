#!/usr/bin/env nextflow
process REPORT_CASES {
    publishDir "results/type_cases/report_cases", pattern: "*.png", mode: 'copy'
    input:
    path summaries
    output:
    path "cases.png"
    script:
    """
    echo ${summaries} | grep -o "[^ ]*-summary\\.tsv" | xargs cat | grep -v '^Type' > input_summary.txt
    Rscript ${workflow.launchDir}/bin/plot_cases.R input_summary.txt
    """
}