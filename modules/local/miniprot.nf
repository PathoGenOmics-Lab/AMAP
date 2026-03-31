process MINIPROT {
    tag "$sampleId"
    cpus 1
    memory '500 MB'

    input:
    tuple val(sampleId), path(fasta), path(proteins_faa), path(proteins_tsv)
    path ref_proteome      // trusted reference proteome (H37Rv default or custom)
    path ref_genome_fasta  // trusted reference genome   (H37Rv default or custom)
    path bad_bakta_faa     // Bakta false positives      (H37Rv default or auto-generated)

    output:
    tuple val(sampleId), path("${sampleId}.paf"), path("${sampleId}-badbakta.paf"), path("${sampleId}-ref.paf"), emit: paf_files

    script:
    """
    # Align reference proteome → sample genome (find matching genes)
    miniprot --outs=0.95 -Iut${task.cpus} -S ${fasta} ${ref_proteome} > ${sampleId}.paf

    # Align Bakta false positives → sample genome (for case1 filtering)
    miniprot --outs=0.95 -Iut${task.cpus} -S ${fasta} ${bad_bakta_faa} > ${sampleId}-badbakta.paf

    # Reverse alignment: sample proteins → reference genome
    miniprot --outs=0.95 -Iut${task.cpus} -S ${ref_genome_fasta} ${proteins_faa} > ${sampleId}-ref.paf
    """
}
