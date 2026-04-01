process MINIPROT {
    tag "$sampleId"
    cpus 4
    memory '2 GB'

    input:
    tuple val(sampleId), path(fasta), path(proteins_faa), path(proteins_tsv)
    path ref_proteome      // trusted reference proteome (H37Rv default or custom)
    path ref_genome_fasta  // trusted reference genome   (H37Rv default or custom)
    path bad_bakta_faa     // Bakta false positives      (H37Rv default or auto-generated)

    output:
    tuple val(sampleId), path("${sampleId}.paf"), path("${sampleId}-badbakta.paf"), path("${sampleId}-ref.paf"), emit: paf_files

    script:
    """
    # Run all three alignments in parallel (independent inputs/outputs)
    miniprot --outs=0.95 -Iut${task.cpus} -S ${fasta} ${ref_proteome} > ${sampleId}.paf &
    miniprot --outs=0.95 -Iut${task.cpus} -S ${fasta} ${bad_bakta_faa} > ${sampleId}-badbakta.paf &
    miniprot --outs=0.95 -Iut${task.cpus} -S ${ref_genome_fasta} ${proteins_faa} > ${sampleId}-ref.paf &
    wait
    """
}
