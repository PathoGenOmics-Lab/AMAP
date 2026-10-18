/*
 * AMAP_REPORT: the interactive report, results/amap_report.html.
 *
 * bin/report/amap_report.py writes one self-contained HTML file from each
 * genome's own outputs, organised by the questions a run answers rather than by
 * the kind of data. It needs only what the correction leaves, so it is written
 * with or without the functional annotation: the master tables, the pangenome
 * and the duplication communities fill their sections when the run has them,
 * and the report says so where it does not. It sits next to final_report.html,
 * which it does not replace yet.
 *
 * Each kind of input is staged in a folder of its own, by the names the
 * pipeline gave the files, so the script finds every genome's files by its
 * name and no list of paths goes through the command line, however many
 * genomes there are.
 */
process AMAP_REPORT {
    // The AMAP code this task runs, in its cache key (see docs/about/contributing.md)
    ext code: { files("${projectDir}/bin/report/{amap_report.py,amap_report.js,amap_report.css,amap_report.html,amap-logo.png}") }
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path amap_tables, stageAs: 'amap/*'             // <name>.AMAP.tsv, from FIX_CDS
    path case_tables, stageAs: 'cases/*'            // TYPE_CASES' summaries, case tables and rvmissing
    path bakta_dirs, stageAs: 'bakta/*'             // Bakta's folder of each genome
    path master_tables, stageAs: 'master/*'         // with the functional annotation, or nothing
    path ref_proteome, stageAs: 'reference/*'       // its headers place every reference gene
    path pangene, stageAs: 'pangene/*'              // REPORT_INPUTS' merged matrix, or nothing
    path communities, stageAs: 'duplications/*'     // GET_ALLDUP's table, or nothing
    path run_info, stageAs: 'run/*'                 // software_versions.tsv, params.yml, samples.txt
    // The Mycobrowser table on the host, or nothing: then params.mycobrowser_ref
    // is a path inside the container, where the default lives.
    path mycobrowser_ref, stageAs: 'mycobrowser/*'

    output:
    path "amap_report.html"

    script:
    def myco_ref = mycobrowser_ref ?: params.mycobrowser_ref
    // The functional block's outputs, when the run has them
    def functional = [pangene ? "--pangene ${pangene}" : '', communities ? "--communities ${communities}" : ''].findAll { arg -> arg }.join(' ')
    """
    # AMAP code in this task's cache key: ${task.ext.code.collect { f -> f.name }.join(' ')}
    set -euo pipefail

    # Mycobrowser names and describes MTBC genes only; without it the report
    # takes what the genomes' own annotation calls them.
    MYCO_REF="${myco_ref}"
    if [ "${params.skip_mycobrowser}" = "true" ] || [ ! -f "\$MYCO_REF" ]; then
        MYCO_REF=""
    fi

    python3 ${projectDir}/bin/report/amap_report.py \\
        --amap amap \\
        --cases cases \\
        --bakta bakta \\
        --master master \\
        --ref-proteome ${ref_proteome} \\
        --mycobrowser "\$MYCO_REF" ${functional} \\
        --versions run/software_versions.tsv \\
        --params run/params.yml \\
        --expected run/samples.txt \\
        --amap-version "${workflow.manifest.version}" \\
        --out amap_report.html
    """
}
