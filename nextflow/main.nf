#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\

    ===============================
       A M A P   P I P E L I N E
    ===============================
            _,--',   _._.--._____
     .--.--';_'-.', ";_      _.,-'
    .'--'.  _.'    {`'-;_ .-.>.'
          '-:_      )  / `' '=.
            ) >     {_/,     /~)
            |/               `^ .' 
    ===============================
 `-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`
    `=`,'=/     `=`,'=/     `=`,'=/   
      y==/        y==/        y==/    
    ,=,-<=`.    ,=,-<=`.    ,=,-<=`.  
 ,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,
    fastas: ${params.tsv}
 """

Channel
    .fromPath(params.tsv)
    .splitCsv(header: true, sep: '\t')
    .map { row -> tuple(row.sampleId, file(row.fasta)) }
    .set { sample_fasta }

process MTB_TYPER {
    publishDir "results/mtb_typer", mode: 'copy'
    input:
    tuple val(sampleId), val(fasta)
    output:
    tuple val(sampleId), path("${sampleId}-mtb_typer.txt")
    script:
    """
    echo "Procesando muestra: $sampleId con Fasta: $fasta"
    python3 ${workflow.launchDir}/mtb_typer.py -i $fasta -t "/reference/snps_lineages.tsv" -r "/reference/MTB_ancestor_reference.fasta"
    """
}

process BAKTA_ANNOTATOR {
    publishDir "results/bakta_annotation", mode: 'copy'
    input:
    tuple val(sampleId), path(tsvFile)
    output:
    tuple val(sampleId), path("*")
    script:
    """
    read -r fastafile prefix < <(awk -F '\\t' 'NR == 1 {print \$1, \$2}' $tsvFile)
    echo "\$prefix\t\$fastafile"

    bakta --db /data/bakta_db/db-light --skip-trna --output "\$prefix" --prefix "\$prefix" --locus-tag "\$prefix" --proteins /reference/sequence.gb "\$fastafile"
    """
}

process MINIPROT {
    publishDir "results/miniprot_reference", mode: 'copy'
    input:
    tuple val(sampleId), path(tsvFile), path(annotator)
    output:
    tuple val(sampleId), path("${sampleId}.paf"), path("${sampleId}-badbakta.paf"), path("${annotator}-H37Rv.paf")
    script:
    """
    read -r fastafile prefix < <(awk -F '\\t' 'NR == 1 {print \$1, \$2}' $tsvFile)

    miniprot --outs=0.95 -Iut16 -S "\$fastafile" /reference/aa_H37Rv.translated.fas > ${sampleId}.paf
    miniprot --outs=0.95 -Iut16 -S "\$fastafile" /reference/bad.bakta.translated.fas > ${sampleId}-badbakta.paf
    sed -i '/^[^>]/ s/\$/*/g' ${annotator}/${annotator}.faa
    miniprot --outs=0.95 -Iut16 -S /reference/new_refH37Rv.fasta  ${annotator}/${annotator}.faa > ${annotator}-H37Rv.paf
    """
}

process TYPE_CASES {
    publishDir "results/type_cases/summary", pattern: "*summary.tsv", mode: 'copy'
    publishDir "results/type_cases/cases", pattern: "*case*.tsv", mode: 'copy'
    publishDir "results/type_cases/cases", pattern: "*rvmissing.tsv", mode: 'copy'
    input:
    tuple val(sampleId), path(annotator), path(miniprot), path(badbakta), path(cds)
    output:
    path "*"
    script:
    """
    python3 ${workflow.launchDir}/bin/classify_cases.py -b ${annotator}/${annotator}.tsv -p ${miniprot} -bb ${badbakta} -c ${cds}
    """
}

process REPORT_CASES {
    publishDir "results/type_cases/report_cases", pattern: "*.png", mode: 'copy'
    input:
    path summaries
    output:
    path "cases.png"
    script:
    """
    echo ${summaries} | grep -o "[^ ]*-summary\\.tsv" | xargs cat | grep -v '^Type' > input_summary.txt
    Rscript ${workflow.launchDir}/bin/plot_cases.R input_summary.txt
    """
}

workflow {
    MTB_TYPER(sample_fasta)
        .set { result_mtb_typer }
    BAKTA_ANNOTATOR(result_mtb_typer)
        .set { result_bakta_annotator } 
    MINIPROT(result_mtb_typer.join(result_bakta_annotator))
        .set { result_miniprot_reference }
    TYPE_CASES(result_bakta_annotator.join(result_miniprot_reference))
        .set { result_cases }
    REPORT_CASES(result_cases.collect())
}