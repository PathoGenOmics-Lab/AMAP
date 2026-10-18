process REPORT_INPUTS {
    // The AMAP code this task runs, in its cache key (see docs/about/contributing.md)
    ext code: { files("${projectDir}/bin/merge_rtabs.py") }

    input:
    path summaries_annot
    path summaries_cases
    path rvmissing_files
    path mastertable_files
    path pangene_rtab_files
    path qc_stats_files

    output:
    tuple path("input_cog.txt"),
          path("input_myco.txt"),
          path("input_summary.txt"),
          path("input_noref.txt"),
          path("input_rvmissing.txt"),
          path("input_mastertables.txt"),
          path("input_pangene.txt"),
          path("input_qc.txt")

    script:
    """
    # AMAP code in this task's cache key: ${task.ext.code.collect { f -> f.name }.join(' ')}
    # --- SCALABLE FILE COLLECTION ---
    # Uses find + glob instead of echo|tr|grep (which breaks at ~500 files
    # due to ARG_MAX and is O(N²) with string splitting).

    # COG / MycoBrowser / noRef — find by suffix, strip headers, concatenate
    find . -maxdepth 1 -name '*AMAP.COG.tsv' -print0 | sort -z | xargs -0 -I{} tail -n +2 {} > input_cog.txt
    find . -maxdepth 1 -name '*AMAP.MycoBrowser.tsv' -print0 | sort -z | xargs -0 -I{} tail -n +2 {} > input_myco.txt
    find . -maxdepth 1 -name '*AMAP.noRef.tsv' -print0 | sort -z | xargs -0 -I{} tail -n +2 {} > input_noref.txt

    # Case summaries — strip header from each file
    find . -maxdepth 1 -name '*-summary.tsv' -print0 | sort -z | xargs -0 -I{} tail -n +2 {} > input_summary.txt

    # Rv missing — add sample column (extracted from filename)
    touch input_rvmissing.txt
    find . -maxdepth 1 -name '*-rvmissing.tsv' -print0 | sort -z | while IFS= read -r -d '' f; do
        sample=\$(basename "\$f" -rvmissing.tsv | sed 's|^\\./||')
        tail -n +2 "\$f" | awk -v s="\$sample" -F'\\t' '{print \$0 "\\t" s}'
    done >> input_rvmissing.txt

    # Master tables — use first NON-EMPTY file for header, strip rest
    touch input_mastertables.txt
    header_written=false
    find . -maxdepth 1 -name '*AMAPmastertable.tsv' -print0 | sort -z | while IFS= read -r -d '' f; do
        [ ! -s "\$f" ] && continue
        if ! \$header_written; then
            cat "\$f"
            header_written=true
        else
            tail -n +2 "\$f"
        fi
    done >> input_mastertables.txt

    # Pangene Rtab — merge with outer join (handles different gene families)
    RTAB_COUNT=\$(find . -maxdepth 1 -name '*.Rtab' | wc -l)
    if [ "\$RTAB_COUNT" -gt 0 ]; then
        find . -maxdepth 1 -name '*.Rtab' -print0 | sort -z | xargs -0 python3 ${projectDir}/bin/merge_rtabs.py --output input_pangene.txt
    else
        touch input_pangene.txt
    fi

    # QC stats — concatenate per-sample QC TSVs (keep first header only)
    touch input_qc.txt
    header_written=false
    find . -maxdepth 1 -name '*.qc.tsv' -print0 | sort -z | while IFS= read -r -d '' f; do
        [ ! -s "\$f" ] && continue
        sample=\$(basename "\$f" -AMAP.qc.tsv)
        if ! \$header_written; then
            head -1 "\$f" | awk '{print "sample\\t" \$0}'
            header_written=true
        fi
        tail -n +2 "\$f" | awk -v s="\$sample" '{print s "\\t" \$0}'
    done >> input_qc.txt
    """
}

process FINAL_REPORT {
    // The AMAP code this task runs, in its cache key (see docs/about/contributing.md)
    ext code: { files("${projectDir}/bin/report/{AMAP.png,PGO.png,header.html,styles.css,report_AMAP.Rmd}") }
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple path(input_cog), path(input_myco), path(input_summary), path(input_noref),
          path(input_rvmissing), path(input_mastertables), path(input_pangene), path(input_qc)
    tuple path(dup_communities_tsv), path(dup_communities_html)

    output:
    path "final_report.html"

    script:
    """
    # AMAP code in this task's cache key: ${task.ext.code.collect { f -> f.name }.join(' ')}
    cp ${projectDir}/bin/report/AMAP.png .
    cp ${projectDir}/bin/report/PGO.png .
    cp ${projectDir}/bin/report/header.html .
    cp ${projectDir}/bin/report/styles.css .
    cp ${projectDir}/bin/report/report_AMAP.Rmd .
    WORK_DIR=`pwd`
    echo \$WORK_DIR
    echo "▸ [REPORT] Rendering HTML report..." >&2
    Rscript -e "rmarkdown::render('report_AMAP.Rmd', output_file='final_report.html',params=list(version='v${workflow.manifest.version}', cog_file='\$WORK_DIR/${input_cog}', myco_file='\$WORK_DIR/${input_myco}', summary_file='\$WORK_DIR/${input_summary}', empty_file='\$WORK_DIR/${input_noref}', rvmissing_file='\$WORK_DIR/${input_rvmissing}', mastertable_file='\$WORK_DIR/${input_mastertables}', pangene_file='\$WORK_DIR/${input_pangene}', dup_file='\$WORK_DIR/${dup_communities_tsv}', qc_file='\$WORK_DIR/${input_qc}'))" \
    || {
        echo "" >&2
        echo "╔══════════════════════════════════════════════════╗" >&2
        echo "║  ERROR: HTML report generation failed.           ║" >&2
        echo "║                                                  ║" >&2
        echo "║  Common causes:                                  ║" >&2
        echo "║  • Missing R packages (ggplot2, dplyr, plotly)   ║" >&2
        echo "║  • Empty input files (0 samples completed)       ║" >&2
        echo "║  • Incompatible R version (need >= 4.1)          ║" >&2
        echo "║                                                  ║" >&2
        echo "║  The annotation outputs (.faa, .gff, .gbk) are  ║" >&2
        echo "║  still valid — only the summary report failed.   ║" >&2
        echo "╚══════════════════════════════════════════════════╝" >&2
        exit 1
    }
    echo "✓ [REPORT] Done: final_report.html" >&2
    # Keep copied files for -resume compatibility
    """
}
