process MINIPROT {
    tag "$sampleId"

    input:
    tuple val(sampleId), path(fasta), path(proteins_faa), path(proteins_tsv)
    path ref_proteome      // trusted reference proteome (H37Rv default or custom)
    path ref_genome_fasta  // trusted reference genome   (H37Rv default or custom)
    path bad_bakta_faa     // Bakta false positives      (H37Rv default or auto-generated)

    output:
    tuple val(sampleId), path("${sampleId}.paf"), path("${sampleId}-badbakta.paf"), path("${sampleId}-ref.paf"), emit: paf_files

    script:
    """
    # Pre-check: verify input files are not empty
    for f in ${ref_proteome} ${bad_bakta_faa} ${proteins_faa}; do
        if [ ! -s "\$f" ] && [ "\$f" != "${bad_bakta_faa}" ]; then
            echo "ERROR: Required input file is empty: \$f" >&2
            exit 1
        fi
    done

    # Run all three alignments in parallel (independent inputs/outputs)
    # Each gets task.cpus/3 threads to stay within the CPU allocation
    THREADS_PER_JOB=\$(( ${task.cpus} / 3 ))
    [ \$THREADS_PER_JOB -lt 1 ] && THREADS_PER_JOB=1

    miniprot --outs=0.95 -Iut\${THREADS_PER_JOB} -S ${fasta} ${ref_proteome} > ${sampleId}.paf &
    pid1=\$!
    miniprot --outs=0.95 -Iut\${THREADS_PER_JOB} -S ${fasta} ${bad_bakta_faa} > ${sampleId}-badbakta.paf &
    pid2=\$!
    miniprot --outs=0.95 -Iut\${THREADS_PER_JOB} -S ${ref_genome_fasta} ${proteins_faa} > ${sampleId}-ref.paf &
    pid3=\$!

    # Wait for each and collect exit codes (don't fail on first error)
    fail=0
    wait \$pid1 || { echo "ERROR: ref proteome alignment failed" >&2; fail=1; }
    wait \$pid2 || { echo "ERROR: bad_bakta alignment failed" >&2; fail=1; }
    wait \$pid3 || { echo "ERROR: reverse alignment failed" >&2; fail=1; }
    exit \$fail
    """
}
