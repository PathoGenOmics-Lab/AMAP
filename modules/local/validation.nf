/**
 * Pre-flight input validation.
 *
 * Checks FASTA format/quality, database paths, and TSV structure BEFORE
 * any expensive computation starts. Runs locally (never on cluster).
 *
 * The genomes come in as inputs, not just the sheet that names them. A
 * container sees only what Nextflow mounts for it, which is the folders of a
 * task's inputs, and a relative path in the sheet is relative to where
 * Nextflow was launched, not to the task. So each genome is read where
 * Nextflow staged it, and the sheet's paths only name them in messages. The
 * same goes for a database on the host; one that lives in the container, as
 * the Bakta default does, is checked there.
 *
 * Exit codes:
 *   0 = passed (warnings are non-fatal)
 *   1 = warnings in strict mode
 *   2 = critical errors → pipeline aborts
 */
process VALIDATE_INPUTS {
    tag 'validation'
    // The AMAP code this task runs, in its cache key (see docs/about/contributing.md)
    ext code: { files("${projectDir}/bin/validate_inputs.py") }
    executor 'local'
    cpus 1
    memory '512.MB'
    errorStrategy 'terminate'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy', pattern: 'validation.log'

    input:
    path tsv_file
    // Each in a numbered directory of its own: two genomes in different
    // folders may have the same file name.
    tuple val(sample_ids), path(genomes, stageAs: 'genomes/*/*')
    // The Bakta database when it is on the host, or nothing (see main.nf)
    path bakta_db, stageAs: 'bakta-db/*'
    // eggNOG's by the path EGGNOG_DB builds it at, '' when it is not used,
    // and the database itself if it is there already
    tuple val(eggnog_db_path), path(eggnog_db, stageAs: 'eggnog-db/*')

    output:
    path 'validation.log'

    script:
    // Nextflow hands over one genome as a path, and several as a list
    def staged = genomes instanceof List ? genomes : [genomes]
    // Each sample ID and where its genome is staged, in single quotes so that
    // the shell passes them on as they are
    def quoted = { v -> "'" + v.toString().replace("'", "'\\''") + "'" }
    def genome_rows = [sample_ids, staged].transpose().collect { id, g -> quoted.call(id) + ' ' + quoted.call(g) }.join(' ')
    // Each database by its path, and where it is staged if it was
    def options = []
    if (params.bakta_db) {
        options << '--bakta-db ' + quoted.call(params.bakta_db)
        if (bakta_db) {
            options << '--bakta-db-staged ' + quoted.call(bakta_db)
        }
    }
    if (eggnog_db_path) {
        options << '--eggnog-db ' + quoted.call(eggnog_db_path)
        if (eggnog_db) {
            options << '--eggnog-db-staged ' + quoted.call(eggnog_db)
        }
    }
    if (params.strict_validation) {
        options << '--strict'
    }
    // The table goes in on standard input. Written to a file, it could share
    // the sample sheet's name and go through the link the sheet is staged as,
    // over the sheet itself. printf is a shell builtin, so no limit on the
    // length of a command line applies to a sheet of thousands of samples.
    """
    # AMAP code in this task's cache key: ${task.ext.code.collect { f -> f.name }.join(' ')}
    printf '%s\\t%s\\n' sampleId fasta ${genome_rows} \\
        | python3 ${projectDir}/bin/validate_inputs.py \\
            --tsv ${tsv_file} \\
            --genomes - \\
            ${options.join(' ')} \\
            2>&1 | tee validation.log
    """
}
