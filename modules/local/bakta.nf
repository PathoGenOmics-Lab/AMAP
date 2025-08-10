process BAKTA_ANNOTATOR {
    tag "$sampleId"
    cpus 12
    memory '6 GB'
    publishDir "${params.outdir}/BAKTA", pattern: "${sampleId}", mode: 'copy'

    input:
    tuple val(sampleId), path(fasta), path(summary_log), val(genus), val(species)

    output:
    path(sampleId), emit: bakta_dir
    tuple val(sampleId), path(fasta), path("${sampleId}/${sampleId}.faa"), path("${sampleId}/${sampleId}.tsv"), emit: for_miniprot
    tuple val(sampleId), path(fasta), path("${sampleId}/${sampleId}.tsv"), path("${sampleId}/${sampleId}.faa"), emit: for_cases
    tuple val(sampleId), path(fasta), path("${sampleId}/${sampleId}.tsv"), emit: for_fix

    script:
    """
    bakta \\
        --db /data/bakta_db/db-light \\
        --skip-trna \\
        --genus "${genus}" \\
        --species "${species}" \\
        --output ${sampleId} \\
        --prefix ${sampleId} \\
        --locus-tag ${sampleId} \\
        --proteins "${params.bakta_proteins}" \\
        ${fasta}
    """
}
