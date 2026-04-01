process PANGENE {
    tag "$sampleId"

    publishDir "${params.outdir}/PANGENE/${sampleId}",
               mode: 'copy',
               pattern: "*.{gfa,bubbles.txt,Rtab,stat}"

    input:
    tuple val(sampleId),
          path(fixtsv),
          path(fixfaa),
          path(fixfna),
          path(baktatsv),
          path(fasta),
          path(emapper)
    path ref_fastas  // staged reference FASTA files from PATHOTYPR

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

    # 2) collect reference genomes (staged via Nextflow channel, flat in workdir)
    #    Exclude the sample's own FASTA to avoid self-alignment
    refs=()
    for f in *.fasta; do
        [[ "\$f" == "${fasta}" ]] && continue
        refs+=("\$f")
    done

    if [ \${#refs[@]} -eq 0 ]; then
        echo "No reference FASTA files – skipping pangene for ${sampleId}"
        echo "H\tVN:1.0" > ${sampleId}.gfa
        touch ${sampleId}.bubbles.txt ${sampleId}.Rtab ${sampleId}.stat
        exit 0
    fi

    # 3) build PAFs — run in parallel with background jobs
    #    Each miniprot gets 2 threads; limit concurrent jobs to fit CPU budget
    THREADS_PER_JOB=2
    MAX_JOBS=\$(( ${task.cpus} / THREADS_PER_JOB ))
    [ \$MAX_JOBS -lt 1 ] && MAX_JOBS=1

    pids=()
    running=0
    for ref in "\${refs[@]}"; do
        base=\$(basename "\$ref" .fasta)
        miniprot --outs=0.97 --no-cs -S -Iut\${THREADS_PER_JOB} \\
                 "\$ref" ${fixfaa} \\
                 > "\${base}_${sampleId}.paf" &
        pids+=(\$!)
        running=\$((running + 1))
        # When at capacity, wait for ALL current jobs before next batch
        # (portable — works on bash 3.2+, unlike wait -n which needs 4.3+)
        if [ \$running -ge \$MAX_JOBS ]; then
            wait "\${pids[@]}"
            pids=()
            running=0
        fi
    done
    wait  # wait for remaining jobs

    # 4) run pangene
    pangene -p0  *.paf  > ${sampleId}.gfa
    k8 /usr/local/bin/pangene.js call       ${sampleId}.gfa > ${sampleId}.bubbles.txt
    k8 /usr/local/bin/pangene.js gfa2matrix ${sampleId}.gfa > ${sampleId}.Rtab
    gfatools stat ${sampleId}.gfa > ${sampleId}.stat

    rm -f *.paf
    """
}
