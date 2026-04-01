process BAKTA_ANNOTATOR {
    tag "$sampleId"
    publishDir "${params.outdir}/BAKTA", pattern: "${sampleId}", mode: 'link', overwrite: true

    input:
    tuple val(sampleId), path(fasta), path(summary_log), val(genus), val(species)

    output:
    path(sampleId), emit: bakta_dir
    tuple val(sampleId), path(fasta), path("${sampleId}/${sampleId}.faa"), path("${sampleId}/${sampleId}.tsv"), emit: for_miniprot
    tuple val(sampleId), path(fasta), path("${sampleId}/${sampleId}.tsv"), path("${sampleId}/${sampleId}.faa"), emit: for_cases
    tuple val(sampleId), path(fasta), path("${sampleId}/${sampleId}.tsv"), emit: for_fix

    script:
    // Use custom reference annotation for --proteins if provided;
    // otherwise fall back to the default H37Rv proteins.
    // Bakta accepts GBK, GFF3, EMBL, and FASTA as --proteins input.
    def proteins_ref = params.ref_annotation ? params.ref_annotation : params.bakta_proteins
    """
    # Pre-check: verify Bakta DB path exists
    if [ ! -d "${params.bakta_db}" ]; then
        echo "ERROR: Bakta database not found at '${params.bakta_db}'" >&2
        echo "Download with: bakta_db download --output /path/to/bakta_db --type light" >&2
        exit 1
    fi

    # Use node-local temp for Bakta's many small intermediate files
    # (avoids NFS metadata storms — Bakta creates ~50 temp files per genome)
    export TMPDIR=\${TMPDIR:-/tmp}/bakta_\$\$
    mkdir -p "\$TMPDIR"
    trap "rm -rf \$TMPDIR" EXIT

    bakta \\
        --db "${params.bakta_db}" \\
        --skip-trna \\
        --threads ${task.cpus} \\
        --genus "${genus}" \\
        --species "${species}" \\
        --output ${sampleId} \\
        --prefix ${sampleId} \\
        --locus-tag ${sampleId} \\
        --proteins "${proteins_ref}" \\
        --tmp-dir "\$TMPDIR" \\
        ${fasta}
    """
}
