process TYPE_CASES {
    tag "$sampleId"

    publishDir "${params.outdir}/CDS_CASES/summary", pattern: "*-summary.tsv", mode: 'link', overwrite: true
    publishDir "${params.outdir}/CDS_CASES/${sampleId}", pattern: "*-case*.tsv", mode: 'link', overwrite: true
    publishDir "${params.outdir}/CDS_CASES/${sampleId}", pattern: "*-rvmissing.tsv", mode: 'link', overwrite: true
    publishDir "${params.outdir}/CDS_CASES/${sampleId}", pattern: "*.html", mode: 'link', overwrite: true
    publishDir "${params.outdir}/CDS_CASES/${sampleId}", pattern: "*_corrected.tsv", mode: 'link', overwrite: true

    input:
    tuple val(sampleId), path(miniprot_paf), path(badbakta_paf), path(h37rv_paf), path(fasta), path(bakta_tsv), path(bakta_faa)

    output:
    tuple val(sampleId), path("${sampleId}_corrected.tsv"), emit: for_fix
    path "*-summary.tsv", emit: summaries
    path "*-case*.tsv"
    path "*-rvmissing.tsv", emit: rvmissing
    path "*.html"

    script:
    """
    echo "▸ [CLASSIFY] Processing ${sampleId}..." >&2
    python3 ${projectDir}/bin/classify_cases.py \\
        -b ${bakta_tsv} \\
        -p ${miniprot_paf} \\
        -bb ${badbakta_paf} \\
        -c ${h37rv_paf} \\
    || {
        echo "ERROR: Classification failed for sample ${sampleId}" >&2
        echo "Input sizes: bakta=\$(wc -l < ${bakta_tsv})L miniprot=\$(wc -l < ${miniprot_paf})L" >&2
        exit 1
    }
    """
}

process FIX_CDS {
    tag "$sampleId"

    publishDir "${params.outdir}/AMAP/${sampleId}", pattern: "*.AMAP.*", mode: 'link', overwrite: true

    input:
    tuple val(sampleId), path(corrected_tsv), path(fasta), path(bakta_tsv)

    output:
    tuple val(sampleId), path("*.AMAP.faa"), emit: for_eggnog
    tuple val(sampleId), path("*.AMAP.tsv"), path("*.AMAP.faa"), path("*.AMAP.fna"), path(bakta_tsv), path(fasta), emit: for_downstream

    script:
    def baseName = corrected_tsv.baseName.replace('_corrected', '')
    """
    python3 ${projectDir}/bin/check_cds.py ${fasta} ${corrected_tsv}
    """
}
