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