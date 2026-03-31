process REPORT_INPUTS {
    input:
    path summaries_annot
    path summaries_cases

    output:
    tuple path("input_cog.txt"),
          path("input_myco.txt"),
          path("input_summary.txt"),
          path("input_noref.txt")

    script:
    """
    # COG: grab every *AMAP.COG.tsv in the bundle
    cat \$(echo ${summaries_annot} | tr ' ' '\\n' | grep 'AMAP.COG.tsv\$') \\
        | grep -v '^lineage' > input_cog.txt

    # MycoBrowser: every *AMAP.MycoBrowser.tsv
    cat \$(echo ${summaries_annot} | tr ' ' '\\n' | grep 'AMAP.MycoBrowser.tsv\$') \\
        | grep -v '^lineage' > input_myco.txt

    cat \$(echo ${summaries_annot} | tr ' ' '\\n' | grep 'AMAP.noRef.tsv\$') \\
        | grep -v '^genome' > input_noref.txt

    # Case summaries
    cat ${summaries_cases} | grep -v '^Type' > input_summary.txt
    """
}

process FINAL_REPORT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple path(input_cog), path(input_myco), path(input_summary), path(input_noref)

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
    Rscript -e "rmarkdown::render('report_AMAP.Rmd', output_file='final_report.html',params=list(cog_file='\$WORK_DIR/${input_cog}', myco_file='\$WORK_DIR/${input_myco}', summary_file='\$WORK_DIR/${input_summary}', empty_file='\$WORK_DIR/${input_noref}'))"
    
    rm AMAP.png 
    rm header.html 
    rm styles.css 
    rm report_AMAP.Rmd 
    """
}
