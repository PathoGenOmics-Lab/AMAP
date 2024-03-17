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