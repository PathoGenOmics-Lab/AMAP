process TYPE_CASES {
    tag "$sampleId"
    cpus 1
    memory '1 GB'

    publishDir "${params.outdir}/CDS_CASES/summary", pattern: "*-summary.tsv", mode: 'copy'
    publishDir "${params.outdir}/CDS_CASES/${sampleId}", pattern: "*-case*.tsv", mode: 'copy'
    publishDir "${params.outdir}/CDS_CASES/${sampleId}", pattern: "*-rvmissing.tsv", mode: 'copy'
    publishDir "${params.outdir}/CDS_CASES/${sampleId}", pattern: "*.html", mode: 'copy'
    publishDir "${params.outdir}/CDS_CASES/${sampleId}", pattern: "*.corrected.tsv", mode: 'copy'

    input:
    tuple val(sampleId), path(miniprot_paf), path(badbakta_paf), path(h37rv_paf), path(fasta), path(bakta_tsv), path(bakta_faa)

    output:
    tuple val(sampleId), path("${sampleId}_corrected.tsv"), emit: for_fix
    path "*-summary.tsv", emit: summaries
    path "*-case*.tsv"
    path "*-rvmissing.tsv"
    path "*.html"

    script:
    """
    python3 ${projectDir}/bin/classify_cases.py \\
        -b ${bakta_tsv} \\
        -p ${miniprot_paf} \\
        -bb ${badbakta_paf} \\
        -c ${h37rv_paf}
    """
}

process FIX_CDS {
    tag "$sampleId"
    cpus 1
    memory '1 GB'

    input:
    tuple val(sampleId), path(corrected_tsv), path(fasta), path(bakta_tsv)

    output:
    tuple val(sampleId), path("*.AMAP.faa"), path("*.AMAP.fna"), emit: for_eggnog
    tuple val(sampleId), path("*.AMAP.tsv"), path("*.AMAP.faa"), path("*.AMAP.fna"), path(bakta_tsv), path(fasta), emit: for_downstream

    script:
    def baseName = corrected_tsv.baseName.replace('_corrected', '')
    """
    python3 ${projectDir}/bin/check_cds.py ${fasta} ${corrected_tsv}
    """
}
