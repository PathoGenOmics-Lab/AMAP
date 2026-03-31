process MINIPROT {
    tag "$sampleId"
    cpus 1
    memory '500 MB'

    input:
    tuple val(sampleId), path(fasta), path(proteins_faa), path(proteins_tsv)
    path ref_proteome      // trusted reference proteome (H37Rv default or custom)
    path ref_genome_fasta  // trusted reference genome   (H37Rv default or custom)

    output:
    tuple val(sampleId), path("${sampleId}.paf"), path("${sampleId}-badbakta.paf"), path("${sampleId}-ref.paf"), emit: paf_files

    script:
    """
    # Align reference proteome → sample genome (find matching genes)
    miniprot --outs=0.95 -Iut${task.cpus} -S ${fasta} ${ref_proteome} > ${sampleId}.paf

    # Align 'bad bakta' known false positives → sample genome
    miniprot --outs=0.95 -Iut${task.cpus} -S ${fasta} "${params.miniprot_ref_bad_bakta}" > ${sampleId}-badbakta.paf

    # Reverse alignment: sample proteins → reference genome (for case1 filtering)
    miniprot --outs=0.95 -Iut${task.cpus} -S ${ref_genome_fasta} ${proteins_faa} > ${sampleId}-ref.paf
    """
}
