process BAKTA_ANNOTATOR {
    tag "$sampleId"
    cpus 12
    memory '6 GB'
    // This will now publish only the output directory created by Bakta.
    publishDir "${params.outdir}/BAKTA", pattern: "${sampleId}", mode: 'copy'

    input:
    tuple val(sampleId), path(fasta), path(summary_log)

    output:
    // This output captures the entire directory for publishing.
    path(sampleId), emit: bakta_dir
    // These outputs remain to channel specific files to downstream processes.
    tuple val(sampleId), path(fasta), path("${sampleId}/${sampleId}.faa"), path("${sampleId}/${sampleId}.tsv"), emit: for_miniprot
    tuple val(sampleId), path(fasta), path("${sampleId}/${sampleId}.tsv"), path("${sampleId}/${sampleId}.faa"), emit: for_cases
    tuple val(sampleId), path(fasta), path("${sampleId}/${sampleId}.tsv"), emit: for_fix

    script:
    """
    bakta \\
        --db /data/bakta_db/db-light \\
        --skip-trna \\
        --genus "Mycobacterium" \\
        --species "Mycobacterium tuberculosis" \\
        --output ${sampleId} \\
        --prefix ${sampleId} \\
        --locus-tag ${sampleId} \\
        --proteins "${params.bakta_proteins}" \\
        ${fasta}
    """
}
