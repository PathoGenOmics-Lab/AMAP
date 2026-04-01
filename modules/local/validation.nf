/**
 * Pre-flight input validation.
 *
 * Checks FASTA format/quality, database paths, and TSV structure BEFORE
 * any expensive computation starts. Runs locally (never on cluster).
 *
 * Exit codes:
 *   0 = passed (warnings are non-fatal)
 *   1 = warnings in strict mode
 *   2 = critical errors → pipeline aborts
 */
process VALIDATE_INPUTS {
    tag 'validation'
    executor 'local'
    cpus 1
    memory '512.MB'
    errorStrategy 'terminate'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy', pattern: 'validation.log'

    input:
    path tsv_file

    output:
    path 'validation.log'

    script:
    def bakta_arg = params.bakta_db ? "--bakta-db '${params.bakta_db}'" : ""
    def eggnog_arg = (!params.skip_functional_annotation && params.eggnog_db) ? "--eggnog-db '${params.eggnog_db}'" : ""
    """
    python3 ${projectDir}/bin/validate_inputs.py \
        --tsv ${tsv_file} \
        ${bakta_arg} \
        ${eggnog_arg} \
        2>&1 | tee validation.log
    """
}
