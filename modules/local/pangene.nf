/*
 * PANGENE — pangenome graph per sample.
 *
 * Architecture for scalability:
 * Each sample's corrected proteins (fixfaa) are aligned against every OTHER
 * sample's genome (ref_fastas) using miniprot. With N samples, this is
 * N×(N-1) alignments total.
 *
 * The reference FASTAs are provided as a SINGLE directory (via stageAs)
 * to avoid Nextflow staging N-1 individual symlinks per sample — which at
 * 500 samples means 249,500 staging operations.
 *
 * Miniprot jobs run in parallel batches within each process instance,
 * bounded by task.cpus to avoid oversubscription.
 */
process PANGENE {
    tag "$sampleId"

    // A closure: the strict parser evaluates publishDir once, when the process
    // is defined, so ${sampleId} has to wait for each task.
    publishDir path: { "${params.outdir}/PANGENE/${sampleId}" },
               mode: 'link', overwrite: true,
               pattern: "*.{gfa,bubbles.txt,Rtab,stat}"

    input:
    tuple val(sampleId),
          path(fixtsv),
          path(fixfaa),
          path(fixfna),
          path(baktatsv),
          path(fasta),
          path(emapper)
    path ref_fastas, stageAs: 'refs/*'  // staged reference FASTA files from PATHOTYPR

    output:
    tuple val(sampleId),
          path("${sampleId}.gfa"),
          path("${sampleId}.bubbles.txt"),
          path("${sampleId}.Rtab"),
          emit: pangene_out

    script:
    """
    set -euo pipefail
    shopt -s nullglob          # empty glob → empty list

    # 1) create AMAP-style links (overwrite on resume)
    ln -sf ${fixfaa} ${sampleId}-AMAP.faa
    ln -sf ${fixfna} ${sampleId}-AMAP.fna

    # 2) collect reference genomes — exclude self
    refs=()
    for f in refs/*.fasta; do
        [[ "\$(basename \$f)" == "${fasta}" ]] && continue
        refs+=("\$f")
    done

    if [ \${#refs[@]} -eq 0 ]; then
        echo "No reference FASTA files – skipping pangene for ${sampleId}"
        echo "H\tVN:1.0" > ${sampleId}.gfa
        touch ${sampleId}.bubbles.txt ${sampleId}.Rtab ${sampleId}.stat
        exit 0
    fi

    echo "==> ${sampleId}: aligning against \${#refs[@]} references (${task.cpus} CPUs)..."

    # 3) build PAFs — parallel with bounded concurrency
    THREADS_PER_JOB=2
    MAX_JOBS=\$(( ${task.cpus} / THREADS_PER_JOB ))
    [ \$MAX_JOBS -lt 1 ] && MAX_JOBS=1

    pids=()
    running=0
    fail=0
    for ref in "\${refs[@]}"; do
        base=\$(basename "\$ref" .fasta)
        miniprot --outs=0.97 --no-cs -S -Iut\${THREADS_PER_JOB} \\
                 "\$ref" ${fixfaa} \\
                 > "\${base}_${sampleId}.paf" &
        pids+=(\$!)
        running=\$((running + 1))
        if [ \$running -ge \$MAX_JOBS ]; then
            for p in "\${pids[@]}"; do
                wait \$p || { echo "ERROR: miniprot PID \$p failed" >&2; fail=1; }
            done
            pids=()
            running=0
        fi
    done
    for p in "\${pids[@]}"; do
        wait \$p || { echo "ERROR: miniprot PID \$p failed" >&2; fail=1; }
    done
    [ \$fail -eq 0 ] || exit 1

    # 4) run pangene
    pangene -p0 *.paf > ${sampleId}.gfa
    k8 /usr/local/bin/pangene.js call       ${sampleId}.gfa > ${sampleId}.bubbles.txt
    k8 /usr/local/bin/pangene.js gfa2matrix ${sampleId}.gfa > ${sampleId}.Rtab
    gfatools stat ${sampleId}.gfa > ${sampleId}.stat

    rm -f *.paf
    """
}
