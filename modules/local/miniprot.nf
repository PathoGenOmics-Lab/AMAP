process MINIPROT {
    tag "$sampleId"
    cpus 1
    memory '500 MB'

    input:
    tuple val(sampleId), path(fasta), path(proteins_faa), path(proteins_tsv)

    output:
    tuple val(sampleId), path("${sampleId}.paf"), path("${sampleId}-badbakta.paf"), path("${sampleId}-H37Rv.paf"), emit: paf_files

    script:
    """
    miniprot --outs=0.95 -Iut16 -S ${fasta} "${params.miniprot_ref_aa}" > ${sampleId}.paf
    miniprot --outs=0.95 -Iut16 -S ${fasta} "${params.miniprot_ref_bad_bakta}" > ${sampleId}-badbakta.paf
    miniprot --outs=0.95 -Iut16 -S "${params.miniprot_ref_h37rv}" ${proteins_faa} > ${sampleId}-H37Rv.paf
    """
}
