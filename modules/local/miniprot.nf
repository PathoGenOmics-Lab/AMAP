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
    # Track PIDs individually so we can report which one fails
    miniprot --outs=0.95 -Iut${task.cpus} -S ${fasta} ${ref_proteome} > ${sampleId}.paf &
    pid1=\$!
    miniprot --outs=0.95 -Iut${task.cpus} -S ${fasta} ${bad_bakta_faa} > ${sampleId}-badbakta.paf &
    pid2=\$!
    miniprot --outs=0.95 -Iut${task.cpus} -S ${ref_genome_fasta} ${proteins_faa} > ${sampleId}-ref.paf &
    pid3=\$!

    # Wait for each and collect exit codes (don't fail on first error)
    fail=0
    wait \$pid1 || { echo "ERROR: ref proteome alignment failed" >&2; fail=1; }
    wait \$pid2 || { echo "ERROR: bad_bakta alignment failed" >&2; fail=1; }
    wait \$pid3 || { echo "ERROR: reverse alignment failed" >&2; fail=1; }
    exit \$fail
    """
}
