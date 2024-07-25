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
    cpus 1
    memory '10 GB'
    output:
    file "eggnog_db.log"

    script:
    // Comprueba si params.eggnog_db es una ruta absoluta. Si no, concatena con projectDir
    def eggnog_db_dir = params.eggnog_db.startsWith("/") ? params.eggnog_db : "${projectDir}/${params.eggnog_db}"

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
    cpus 1
    memory '1 GB'
    publishDir "${params.outdir}/MTB_TYPER/typing", pattern: "*.", mode: 'copy'
    publishDir "${params.outdir}/MTB_TYPER/fastas", pattern: "*.fasta", mode: 'copy'
    input:
    tuple val(sampleId), val(fasta)
    output:
    tuple stdout, path("*.fasta"), path("*.")
    script:
    """
    python3 ${projectDir}/bin/mtb_typer.py -i $fasta -t "/home/ruizro/test_amap/snps_lineages.tsv" -r "/reference/MTB_ancestor_reference.fasta"
    ls *. | tr -d '.' | tr -d '\\n'
    """
}

process BAKTA_ANNOTATOR { //done
    tag "$sampleId"
    cpus 10
    memory '6 GB'
    publishDir "${params.outdir}/BAKTA/", mode: 'copy'
    input:
    tuple val(sampleId), path(fasta), path(nametag)
    output:
    tuple val(sampleId), path(fasta), path("${sampleId}/${sampleId}.faa"), path("${sampleId}/${sampleId}.tsv"), path("${sampleId}/${sampleId}.gff3"), 
    path("${sampleId}/*hypotheticals*"), path("${sampleId}/${sampleId}.png"), path("${sampleId}/${sampleId}.txt"),
    path("${sampleId}/${sampleId}.gbff")
    script:
    """
    bakta --db /data/bakta_db/db-light --skip-trna --genus "Mycobacterium" --species "Mycobacterium tuberculosis" --output ${sampleId} --prefix ${sampleId} --locus-tag ${sampleId} --proteins /reference/sequence.gb ${fasta}
    """
}

process MINIPROT { //done
    tag "$sampleId"
    cpus 1
    memory '500 MB'
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
    cpus 1
    memory '1 GB'
    tag "$sampleId"
    publishDir "${params.outdir}/CDS_CASES/summary", pattern: "*summary.tsv", mode: 'copy'
    publishDir "${params.outdir}/CDS_CASES/${sampleId}", pattern: "*case*.tsv", mode: 'copy'
    publishDir "${params.outdir}/CDS_CASES/${sampleId}", pattern: "*rvmissing.tsv", mode: 'copy'
    publishDir "${params.outdir}/CDS_CASES/${sampleId}", pattern: "*.html", mode: 'copy'
    input:
    tuple val(sampleId), path(miniprot), path(badbakta), path(cds), path(fasta), path(faa), path(tsv)
    output:
    tuple val(sampleId), path("${sampleId}_corrected.tsv"), path("*summary.tsv"),path("*case*.tsv"), path("*rvmissing.tsv"), path("*.html")
    script:
    """
    python3 ${projectDir}/bin/classify_cases2.py -b ${tsv} -p ${miniprot} -bb ${badbakta} -c ${cds}
    """
}

process FIX_CDS {
    cpus 1
    memory '1 GB'
    tag "$sampleId"
    input:
    tuple val(sampleId), path(corrected_tsv), path(summary), path(fasta), path(faa), path(tsv)

    output:
    tuple val(sampleId), path("${sampleId}.AMAP.tsv"), path("${sampleId}.AMAP.faa"), path("${sampleId}.AMAP.fna"), path(tsv), path(fasta)

    script:
    """
    python3 ${projectDir}/bin/check_cds.py ${fasta} ${corrected_tsv}
    """
}

process ANNOT_EGGNOGMAP {
    cpus 32
    memory '10 GB'
    tag "$sampleId"
    input:
    tuple val(sampleId), path(fixtsv), path(fixfaa), path(fixfna), path(baktatsv), path(fasta), path(dbMarker)
    output:
    tuple val(sampleId), path("${sampleId}.emapper.annotations")
    script:
    def eggnog_db_dir = params.eggnog_db.startsWith("/") ? params.eggnog_db : "${projectDir}/${params.eggnog_db}"
    """
    emapper.py  -i ${fixfaa} --cpu 0 -m diamond --sensmode fast -o ${sampleId} --data_dir ${eggnog_db_dir}
    """
}

process GET_DUPLICATIONS {
    tag "$sampleId"
    publishDir "${params.outdir}/DUPLICATIONS", mode: 'copy'
    input:
    tuple val(sampleId), path(fixtsv), path(fixfaa), path(fixfna), path(baktatsv), path(fasta), path(emapper)
    output:
    tuple path("duplications_*"), path("${sampleId}.AMAP.faa.renamed.dup"), path("duplications_${sampleId}_groups.tsv")
    script:
    """
    python3 ${projectDir}/bin/rename_fastas.py ${emapper} ${fixfaa}
    diamond blastp --db ${fixfaa}.renamed -q ${fixfaa}.renamed --outfmt 6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qlen slen nident --no-self-hits --header > duplications_${sampleId}

    python3 ${projectDir}/bin/get_duplicationsv2.py duplications_${sampleId}
    python3 ${projectDir}/bin/dup_ex_multifasta.py duplications_${sampleId}_groups.tsv ${sampleId}.AMAP.faa.renamed
    """
}

process GET_OUTPUTS {
    tag "$sampleId"
    publishDir "${params.outdir}/AMAP/${sampleId}", mode: 'copy'
    input:
    tuple val(sampleId), path(fixtsv), path(fixfaa), path(fixfna), path(baktatsv), path(fasta), path(emapper)
    output:
    tuple val(sampleId),path("${sampleId}-AMAPmastertable.tsv"),path("${sampleId}-AMAP.tsv"),path("${sampleId}-AMAP.gbk"),path("${sampleId}-AMAP.gff3"), path(fixfaa), path(fixfna)
    script:
    """
    sed -i '1d' ${fixtsv}
    python3 ${projectDir}/bin/generate_tsv_gff_gbk.py ${emapper} ${fixtsv} ${baktatsv} /reference/Mycobacterium_tuberculosis_H37Rv_txt_v4.txt ${fasta}
    """
}

process COUNT_ANNOT {
    tag "$sampleId"
    publishDir "${params.outdir}/AMAP/${sampleId}", mode: 'copy'
    input:
    tuple val(sampleId), path(mastertable)
    output:
    path("*.tsv")
    script:
    """
    python3 ${projectDir}/bin/count_annot.py ${mastertable}
    """
}

process REPORT_CASES2 {
    input:
    path(summaries)
    output:
    tuple path("input_cog.txt"), path("input_MycoBrowser.txt")
    script:
    """
    echo ${summaries} | grep -o "[^ ]*.AMAP.COG\\.tsv" | xargs cat | grep -v '^lineage' > input_cog.txt
    echo ${summaries} | grep -o "[^ ]*.AMAP.MycoBrowser\\.tsv" | xargs cat | grep -v '^lineage' > input_MycoBrowser.txt
    """
}
process REPORT_CASES1{
    input:
    path(summaries)
    output:
    path("input_summary.txt")
    script:
    """
    echo ${summaries} | grep -o "[^ ]*-summary\\.tsv" | xargs cat | grep -v '^Type' > input_summary.txt
    """
}
process FINAL_REPORT{
    publishDir "${params.outdir}/", mode: 'copy'
    input:
    tuple path(input_cog), path(input_myco), path(input_summary)
    script:
    """
    echo ${input_cog}
    echo ${input_myco}
    echo ${input_summary}
    WORK_DIR=`pwd`
    echo \$WORK_DIR
    Rscript -e "rmarkdown::render('${projectDir}/bin/check.Rmd', output_file='${projectDir}/${params.outdir}/final_report2.html', params=list(cog_summary='\$WORK_DIR/${input_cog}', myco_summary='\$WORK_DIR/${input_myco}', input_sum='\$WORK_DIR/${input_summary}'))"    
    #Rscript -e "rmarkdown::render('${projectDir}/bin/check_quarto.qmd', output_file='${projectDir}/${params.outdir}/final_report.html', params=list(cog_summary='${input_cog}', myco_summary='${input_myco}', input_sum='${input_summary}'))"    
    #quarto render /home/paula/AMAP/bin/check_quarto.qmd --to html --output 'final_report.html' -P cog_summary:\$WORK_DIR/${input_cog} -P myco_summary:\$WORK_DIR/${input_myco} -P input_sum:\$WORK_DIR/${input_summary}
    """
}

process PASS {
    input:
    tuple val(ids), path(tsv),path(summaries)
    output:
    path(summaries)
    script:
    """
    echo "Summarizing cases"
    """
}

process ALL_DUP {
    cpus 13
    memory '10 GB'
    publishDir "${params.outdir}/finalreport/", mode: 'copy'
    input:
    path(fastas)
    output:
    tuple path("homologous_communities.tsv"), path("homologous_communities.html")
    script:
    """
    cat *.dup > all_fasta.fa
    rm *.dup
    
    diamond cluster -d all_fasta.fa -o all.diamond --approx-id 30 -M 16G
    diamond blastp --db all_fasta.fa -q all_fasta.fa --outfmt 6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qlen slen nident --no-self-hits --header > all.diamond
    python3 ${projectDir}/bin/network_genomes.py all.diamond
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
        .map {sampleId, tsv, summary, cases, miss, html -> tuple(sampleId, tsv, summary)}
        .set {result_typecases}
    FIX_CDS(result_typecases.join(result_bakta_annotator))
        .set {result_fix}
    
    ANNOT_EGGNOGMAP(result_fix.combine(db_done))
        .set {eggnogres}
    GET_OUTPUTS(result_fix.join(eggnogres))
        .map {sampleId, mastertable, tsv, gbk, gff3, fixfaa, fixfna ->tuple(sampleId, mastertable)}
        .set {annotation}
    GET_DUPLICATIONS(result_fix.join(annotation))
        .map {all, fastas, tsv -> fastas}
        .set {res_dup}

    COUNT_ANNOT(annotation)
        .set {count_annotation}
    PASS(result_typecases)
        .set {result_pass}
    REPORT_CASES1(result_pass.collect())
        .set {result_fix1}
    REPORT_CASES2(count_annotation.collect())
        .set {result_fix2}
    FINAL_REPORT(result_fix2.combine(result_fix1))
    ALL_DUP(res_dup.collect())
    }
