process MINIPROT {
    tag "$sampleId"
    // The AMAP code this task runs, in its cache key (see docs/about/contributing.md)
    ext code: { files("${projectDir}/bin/contig_end_pieces.py") }

    input:
    tuple val(sampleId), path(fasta), path(proteins_faa), path(proteins_tsv)
    path ref_proteome      // trusted reference proteome (H37Rv default or custom)
    path ref_genome_fasta  // trusted reference genome   (H37Rv default or custom)
    path bad_bakta_faa     // Bakta false positives      (H37Rv default or auto-generated)

    output:
    tuple val(sampleId), path("${sampleId}.paf"), path("${sampleId}-badbakta.paf"), path("${sampleId}-ref.paf"), path("${sampleId}-pieces.paf"), emit: paf_files

    script:
    """
    # AMAP code in this task's cache key: ${task.ext.code.collect { f -> f.name }.join(' ')}
    echo "▸ [MINIPROT] Aligning ${sampleId}..." >&2

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
    [ \$fail -eq 0 ] || exit 1

    # A gene a contig end cut is placed by its larger piece alone: --outs=0.95
    # keeps only the placements within 5% of a protein's best. What each such
    # placement lacks is aligned as a protein of its own, whose best placement
    # is the piece at the other contig end; -L 10 and -n 2 let miniprot index
    # and chain the short reading frames an edge leaves there.
    python3 ${projectDir}/bin/contig_end_pieces.py ${sampleId}.paf ${ref_proteome} > ${sampleId}-pieces.faa
    miniprot --outs=0.95 -L 10 -n 2 -Iut${task.cpus} -S ${fasta} ${sampleId}-pieces.faa > ${sampleId}-pieces.paf \\
        || { echo "ERROR: alignment of the pieces contig ends cut off failed" >&2; exit 1; }
    """
}
