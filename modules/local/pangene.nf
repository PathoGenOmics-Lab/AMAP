process PANGENE {
    tag "$sampleId"
    cpus   16
    memory '16 GB'

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
    refs=( *.fasta )

    if [ \${#refs[@]} -eq 0 ]; then
        echo "No reference FASTA files – skipping pangene for ${sampleId}"
        echo "H\tVN:1.0" > ${sampleId}.gfa
        touch ${sampleId}.bubbles.txt ${sampleId}.Rtab ${sampleId}.stat
        exit 0
    fi

    # 3) build PAFs
    for ref in "\${refs[@]}"; do
        base=\$(basename "\$ref" .fasta)
        miniprot --outs=0.97 --no-cs -S -Iut${task.cpus} \
                 "\$ref"  ${fixfaa} \
                 > "\${base}_${sampleId}.paf"
    done

    # 4) run pangene
    pangene -p0  *.paf  > ${sampleId}.gfa
    k8 /usr/local/bin/pangene.js call       ${sampleId}.gfa > ${sampleId}.bubbles.txt
    k8 /usr/local/bin/pangene.js gfa2matrix ${sampleId}.gfa > ${sampleId}.Rtab
    gfatools stat ${sampleId}.gfa > ${sampleId}.stat

    rm -f *.paf
    """
}
