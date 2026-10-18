process TYPE_CASES {
    tag "$sampleId"
    // The AMAP code this task runs, in its cache key (see docs/about/contributing.md)
    ext code: { files("${projectDir}/bin/{classify_cases.py,classify/*.py}") }

    publishDir "${params.outdir}/CDS_CASES/summary", pattern: "*-summary.tsv", mode: 'link', overwrite: true
    // Per-sample paths are closures: the strict parser evaluates publishDir once,
    // when the process is defined, so ${sampleId} has to wait for each task.
    publishDir path: { "${params.outdir}/CDS_CASES/${sampleId}" }, pattern: "*-case*.tsv", mode: 'link', overwrite: true
    publishDir path: { "${params.outdir}/CDS_CASES/${sampleId}" }, pattern: "*-rvmissing.tsv", mode: 'link', overwrite: true
    publishDir path: { "${params.outdir}/CDS_CASES/${sampleId}" }, pattern: "*.html", mode: 'link', overwrite: true
    publishDir path: { "${params.outdir}/CDS_CASES/${sampleId}" }, pattern: "*_corrected.tsv", mode: 'link', overwrite: true

    input:
    tuple val(sampleId), path(miniprot_paf), path(badbakta_paf), path(h37rv_paf), path(fasta), path(bakta_tsv), path(bakta_faa)
    path ref_proteome      // the proteome miniprot aligned, to place exactly what it could not

    output:
    tuple val(sampleId), path("${sampleId}_corrected.tsv"), emit: for_fix
    path "*-summary.tsv", emit: summaries
    path "*-case*.tsv"
    path "*-rvmissing.tsv", emit: rvmissing
    path "*.html"

    script:
    """
    # AMAP code in this task's cache key: ${task.ext.code.collect { f -> f.name }.join(' ')}
    echo "▸ [CLASSIFY] Processing ${sampleId}..." >&2
    python3 ${projectDir}/bin/classify_cases.py \\
        -b ${bakta_tsv} \\
        -p ${miniprot_paf} \\
        -bb ${badbakta_paf} \\
        -c ${h37rv_paf} \\
        --genome ${fasta} \\
        --ref-proteome ${ref_proteome} \\
    || {
        echo "ERROR: Classification failed for sample ${sampleId}" >&2
        echo "Input sizes: bakta=\$(wc -l < ${bakta_tsv})L miniprot=\$(wc -l < ${miniprot_paf})L" >&2
        exit 1
    }
    """
}

process FIX_CDS {
    tag "$sampleId"
    // The AMAP code this task runs, in its cache key (see docs/about/contributing.md)
    ext code: { files("${projectDir}/bin/{check_cds.py,deduplicate_rvs.py,classify/*.py}") }

    publishDir path: { "${params.outdir}/AMAP/${sampleId}" }, pattern: "*.AMAP.*", mode: 'link', overwrite: true

    input:
    tuple val(sampleId), path(corrected_tsv), path(fasta), path(bakta_tsv)
    path ref_proteome      // trusted reference proteome (H37Rv default or custom)

    output:
    tuple val(sampleId), path("*.AMAP.faa"), emit: for_eggnog
    tuple val(sampleId), path("*.AMAP.tsv"), path("*.AMAP.faa"), path("*.AMAP.fna"), path(bakta_tsv), path(fasta), emit: for_downstream

    script:
    def baseName = corrected_tsv.baseName.replace('_corrected', '')
    """
    # AMAP code in this task's cache key: ${task.ext.code.collect { f -> f.name }.join(' ')}
    # The reference proteome tells reference pseudogenes from genes the sample
    # cut short, and how much of a protein an alignment left out.
    python3 ${projectDir}/bin/check_cds.py ${fasta} ${corrected_tsv} \\
        --ref-proteome ${ref_proteome}

    # Deduplicate genes with the same Rv assignment, keeping real copies
    python3 ${projectDir}/bin/deduplicate_rvs.py \\
        *.AMAP.tsv \\
        --faa *.AMAP.faa \\
        --fna *.AMAP.fna \\
        --bakta ${bakta_tsv} \\
        --ref-proteome ${ref_proteome}
    """
}
