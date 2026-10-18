process PATHOTYPR {
    tag "$sampleId"
    // The AMAP code this task runs, in its cache key (see docs/about/contributing.md)
    ext code: { files("${projectDir}/bin/name_by_lineage.py") }

    publishDir "${params.outdir}/PATHOTYPR/typing", pattern: "*_summary.tsv", mode: 'link', overwrite: true
    publishDir "${params.outdir}/PATHOTYPR/fastas", pattern: "*.fasta", mode: 'link', overwrite: true

    input:
    // Staged in a directory of its own: the genome is written out under its
    // new name, and a sample sheet that points sample G26689 at L4_G26689.fasta
    // would otherwise have it written over the input it was read from.
    tuple val(sampleId), path(fasta, stageAs: 'input/*')

    output:
    tuple env('SAMPLE_NAME'), path("*.fasta"), path("*_summary.tsv"), emit: typed_fasta

    script:
    """
    #!/bin/bash
    # AMAP code in this task's cache key: ${task.ext.code.collect { f -> f.name }.join(' ')}
    set -o pipefail
    # NOTE: do NOT use 'set -e' here. If the lineage call fails, the genome
    # is still named, NA_<sampleId>, and the run goes on without a lineage.

    # 1. Run the lineage classifier (may fail for non-MTBC genomes). The
    # subcommand is quoted: Nextflow adds to a task's cache key any bin/ entry
    # its script names as a bare word, and bare, this one named AMAP's own
    # bin/classify package, which this task never runs.
    pathotypr 'classify' \\
        -i ${fasta} \\
        -m "${params.pathotyper_snps}" \\
        -r "${params.pathotyper_ref}" \\
        --output-prefix "${sampleId}" || true

    # 2. Ensure summary file exists (even if empty): Nextflow output glob requires it.
    if [ ! -f "${sampleId}_summary.tsv" ]; then
        echo -e "Sample\\tLineage\\tName" > "${sampleId}_summary.tsv"
    fi

    # 3. Name the genome <lineage>_<sampleId>, or NA_<sampleId> without a
    # lineage. Never from the rows of the summary themselves: there is one per
    # contig with markers, named after the contig, so a draft came out named
    # after whichever of its contigs was listed first.
    SAMPLE_NAME=\$(python3 ${projectDir}/bin/name_by_lineage.py "${sampleId}" "${sampleId}_summary.tsv" "${fasta}") || exit 1
    export SAMPLE_NAME
    """
}
