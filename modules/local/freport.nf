process REPORT_INPUTS {
    input:
    path summaries_annot
    path summaries_cases
    path rvmissing_files
    path mastertable_files
    path pangene_rtab_files

    output:
    tuple path("input_cog.txt"),
          path("input_myco.txt"),
          path("input_summary.txt"),
          path("input_noref.txt"),
          path("input_rvmissing.txt"),
          path("input_mastertables.txt"),
          path("input_pangene.txt")

    script:
    """
    # Helper: safe grep that doesn't fail with set -e when no matches found
    safe_grep() { grep "\$@" || true; }

    # COG: grab every *AMAP.COG.tsv in the bundle, strip header lines exactly
    for f in \$(echo ${summaries_annot} | tr ' ' '\\n' | safe_grep 'AMAP.COG.tsv\$'); do
        tail -n +2 "\$f"
    done > input_cog.txt

    # MycoBrowser: every *AMAP.MycoBrowser.tsv
    for f in \$(echo ${summaries_annot} | tr ' ' '\\n' | safe_grep 'AMAP.MycoBrowser.tsv\$'); do
        tail -n +2 "\$f"
    done > input_myco.txt

    for f in \$(echo ${summaries_annot} | tr ' ' '\\n' | safe_grep 'AMAP.noRef.tsv\$'); do
        tail -n +2 "\$f"
    done > input_noref.txt

    # Case summaries — strip header from each file
    for f in ${summaries_cases}; do
        tail -n +2 "\$f"
    done > input_summary.txt

    # Rv missing — concatenate all rvmissing TSVs (strip headers)
    touch input_rvmissing.txt
    for f in \$(echo ${rvmissing_files} | tr ' ' '\\n' | safe_grep 'rvmissing.tsv\$'); do
        tail -n +2 "\$f" >> input_rvmissing.txt
    done

    # Master tables — concatenate with headers for first file only
    touch input_mastertables.txt
    first=true
    for f in \$(echo ${mastertable_files} | tr ' ' '\\n' | safe_grep 'AMAPmastertable.tsv\$'); do
        if \$first; then
            cat "\$f" >> input_mastertables.txt
            first=false
        else
            tail -n +2 "\$f" >> input_mastertables.txt
        fi
    done

    # Pangene Rtab — concatenate if available
    touch input_pangene.txt
    for f in \$(echo ${pangene_rtab_files} | tr ' ' '\\n' | safe_grep '.Rtab\$'); do
        cat "\$f" >> input_pangene.txt
    done
    """
}

process FINAL_REPORT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple path(input_cog), path(input_myco), path(input_summary), path(input_noref),
          path(input_rvmissing), path(input_mastertables), path(input_pangene)
    tuple path(dup_communities_tsv), path(dup_communities_html)

    output:
    path "final_report.html"

    script:
    """
    cp ${projectDir}/bin/report/AMAP.png .
    cp ${projectDir}/bin/report/header.html .
    cp ${projectDir}/bin/report/styles.css .
    cp ${projectDir}/bin/report/report_AMAP.Rmd .
    WORK_DIR=`pwd`
    echo \$WORK_DIR
    Rscript -e "rmarkdown::render('report_AMAP.Rmd', output_file='final_report.html',params=list(cog_file='\$WORK_DIR/${input_cog}', myco_file='\$WORK_DIR/${input_myco}', summary_file='\$WORK_DIR/${input_summary}', empty_file='\$WORK_DIR/${input_noref}', rvmissing_file='\$WORK_DIR/${input_rvmissing}', mastertable_file='\$WORK_DIR/${input_mastertables}', pangene_file='\$WORK_DIR/${input_pangene}', dup_file='\$WORK_DIR/${dup_communities_tsv}'))"
    # Keep copied files for -resume compatibility
    """
}
