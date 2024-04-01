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
    fastas:     ${params.tsv}
    output:     ${params.outdir}
    eggnog_db:  ${params.eggnog_db}
 """

process EGGNOG_DB {
    output:
    file "eggnog_db.log"

    script:
    // Comprueba si params.eggnog_db es una ruta absoluta. Si no, concatena con workflow.launchDir
    def eggnog_db_dir = params.eggnog_db.startsWith("/") ? params.eggnog_db : "${workflow.launchDir}/${params.eggnog_db}"

    """
    mkdir -p ${eggnog_db_dir}
    if [ ! -f ${eggnog_db_dir}/eggnog.db ]; then
        download_eggnog_data.py -y --data_dir ${eggnog_db_dir}
    fi
    echo "eggNOG DB setup complete in ${eggnog_db_dir}" > eggnog_db.log
    """
}

Channel
    .fromPath(params.tsv)
    .splitCsv(header: true, sep: '\t')
    .map { row -> tuple(row.sampleId, file(row.fasta)) }
    .set { sample_fasta }

process MTB_TYPER { //done
    tag "$sampleId"
    publishDir "${params.outdir}/mtb_typer/typing", pattern: "*.", mode: 'copy'
    publishDir "${params.outdir}/mtb_typer/fastas", pattern: "*.fasta", mode: 'copy'
    input:
    tuple val(sampleId), val(fasta)
    output:
    tuple stdout, path("*.fasta"), path("*.")
    script:
    """
    python3 ${workflow.launchDir}/bin/mtb_typer.py -i $fasta -t "/reference/snps_lineages.tsv" -r "/reference/MTB_ancestor_reference.fasta"
    ls *. | tr -d '.' | tr -d '\\n'
    """
}

process BAKTA_ANNOTATOR { //done
    tag "$sampleId"
    publishDir "${params.outdir}/bakta_annotation/", mode: 'copy'
    input:
    tuple val(sampleId), path(fasta), path(nametag)
    output:
    tuple val(sampleId), path(fasta), path("${sampleId}/${sampleId}.faa"), path("${sampleId}/${sampleId}.tsv"), path("${sampleId}/${sampleId}.gff3"), 
    path("${sampleId}/*hypotheticals*"), path("${sampleId}/${sampleId}.png"), path("${sampleId}/${sampleId}.txt"),
    path("${sampleId}/${sampleId}.gbff")
    //tuple val(sampleId), path(fasta), path("${sampleId}/${sampleId}.faa"), path("${sampleId}/${sampleId}.tsv")
    script:
    """
    bakta --db /data/bakta_db/db-light --skip-trna --genus "Mycobacterium" --species "Mycobacterium tuberculosis" --output ${sampleId} --prefix ${sampleId} --locus-tag ${sampleId} --proteins /reference/sequence.gb ${fasta}
    """
}

process MINIPROT { //done
    tag "$sampleId"
    input:
    tuple val(sampleId), path(fasta), path(proteins), path(tsv)
    output:
    tuple val(sampleId),path("${sampleId}.paf"), path("${sampleId}-badbakta.paf"), path("${sampleId}-H37Rv.paf")
    script:
    """
    miniprot --outs=0.95 -Iut16 -S ${fasta} /reference/aa_H37Rv.translated.fas > ${sampleId}.paf
    miniprot --outs=0.95 -Iut16 -S ${fasta} /reference/bad.bakta.translated.fas > ${sampleId}-badbakta.paf
    miniprot --outs=0.95 -Iut16 -S /reference/new_refH37Rv.fasta  ${proteins} > ${sampleId}-H37Rv.paf
    """
}

process TYPE_CASES {
    tag "$sampleId"
    publishDir "${params.outdir}/type_cases/summary", pattern: "*summary.tsv", mode: 'copy'
    publishDir "${params.outdir}/type_cases/${sampleId}", pattern: "*case*.tsv", mode: 'copy'
    publishDir "${params.outdir}/type_cases/${sampleId}", pattern: "*rvmissing.tsv", mode: 'copy'
    publishDir "${params.outdir}/type_cases/${sampleId}", pattern: "*.html", mode: 'copy'
    input:
    tuple val(sampleId), path(miniprot), path(badbakta), path(cds), path(fasta), path(faa), path(tsv)
    output:
    //tuple val(sampleId), path("${sampleId}_corrected.tsv")
    tuple val(sampleId), path("${sampleId}_corrected.tsv"), path("*summary.tsv"),path("*case*.tsv"), path("*rvmissing.tsv"), path("*.html")
    script:
    """
    python3 ${workflow.launchDir}/bin/classify_cases2.py -b ${tsv} -p ${miniprot} -bb ${badbakta} -c ${cds}
    """
}

process FIX_CDS {
    tag "$sampleId"
    input:
    tuple val(sampleId), path(corrected_tsv), path(fasta), path(faa), path(tsv)

    output:
    tuple val(sampleId), path("${sampleId}.AMAP.tsv"), path("${sampleId}.AMAP.faa"), path("${sampleId}.AMAP.fna"), path(tsv), path(fasta)

    script:
    """
    python3 ${workflow.launchDir}/bin/check_cds.py ${fasta} ${corrected_tsv}
    """
}

process ANNOT_EGGNOGMAP {
    tag "$sampleId"
    input:
    tuple val(sampleId), path(fixtsv), path(fixfaa), path(fixfna), path(baktatsv), path(fasta), path(dbMarker)
    output:
    path("${sampleId}.emapper.annotations")
    script:
    def eggnog_db_dir = params.eggnog_db.startsWith("/") ? params.eggnog_db : "${workflow.launchDir}/${params.eggnog_db}"
    """
    emapper.py  -i ${fixfaa} --cpu 0 -m diamond --sensmode fast -o ${sampleId} --data_dir ${eggnog_db_dir}
    """
}

process GET_DUPLICATIONS {
    tag "$sampleId"
    publishDir "${params.outdir}/duplications", mode: 'copy'
    input:
    tuple val(sampleId), path(fixtsv), path(fixfaa), path(fixfna), path(baktatsv), path(fasta)
    output:
    path("duplications_*")
    script:
    """
    makeblastdb -in ${fixfaa} -dbtype prot -out ${sampleId}_faadb
    blastp -query ${fixfaa} -db ${sampleId}_faadb -out duplications_${sampleId} -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen nident"

    python3 ${workflow.launchDir}/bin/get_duplications.py duplications_${sampleId}
    """
}

process GET_OUTPUTS {
    tag "$sampleId"
    publishDir "${params.outdir}/AMAP/", mode: 'copy'
    input:
    tuple val(sampleId), path(fixtsv), path(fixfaa), path(fixfna), path(baktatsv), path(fasta), path(emapper)
    output:
    tuple path("${sampleId}-AMAPmastertable.tsv"),path("${sampleId}-AMAP.tsv"),path("${sampleId}-AMAP.gbk"),path("${sampleId}-AMAP.gff3")
    script:
    """
    sed -i '1d' ${fixtsv}
    python3 ${workflow.launchDir}/bin/generate_tsv_gff_gbk.py ${emapper} ${fixtsv} ${baktatsv} /reference/Mycobacterium_tuberculosis_H37Rv_txt_v4.txt ${fasta}
    """
}

process REPORT_CASES {
    publishDir "${params.outdir}/type_cases/report_cases", pattern: "*.png", mode: 'copy'
    input:
    tuple val(ids),path (summaries)
    output:
    path "cases.png"
    script:
    """
    echo ${summaries} | grep -o "[^ ]*-summary\\.tsv" | xargs cat | grep -v '^Type' > input_summary.txt
    Rscript ${workflow.launchDir}/bin/plot_cases.R input_summary.txt
    """
}

workflow {
    EGGNOG_DB()
        .set { db_done } 
    MTB_TYPER(sample_fasta)
        .set { result_mtb_typer }
    BAKTA_ANNOTATOR(result_mtb_typer)
        .map { sampleId, fasta, proteins, tsv, gff, hyp, png, txt, gb -> tuple(sampleId, fasta, proteins, tsv) }
        .set { result_bakta_annotator } 
    MINIPROT(result_bakta_annotator)
        .set {result_miniprot}
    TYPE_CASES(result_miniprot.join(result_bakta_annotator))
        .map {sampleId, tsv, summary, cases, miss, html -> tuple(sampleId, tsv)}
        .set {result_typecases}
    FIX_CDS(result_typecases.join(result_bakta_annotator))
        .set {result_fix}
    GET_DUPLICATIONS(result_fix)
    ANNOT_EGGNOGMAP(result_fix.combine(db_done))
        .set {eggnogres}
    GET_OUTPUTS(result_fix.combine(eggnogres))
    //REPORT_CASES(result_fix.collect())
    }
    
    //.map { sampleId, ctsv, summary, case, miss, html -> tuple(sampleId, ctsv) }
    //MINIPROT(result_mtb_typer.join(result_bakta_annotator))
    //    .set { result_miniprot_reference }
    //TYPE_CASES(result_bakta_annotator.join(result_miniprot_reference))
    //    .set { result_cases }
    //FIX_CDS(result_cases)
    //    .set { result_cds }
    
    //REPORT_CASES(result_cds.collect())
