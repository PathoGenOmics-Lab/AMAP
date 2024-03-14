#!/usr/bin/env nextflow

params.input = 'data/*.txt' // Valor por defecto


process myProcess {
    input:
    path file
    output:
    stdout

    script:
    """
    cat ${file}
    bakta --version
    """
}

workflow {
    myProcess(file(params.input)) | view
}
