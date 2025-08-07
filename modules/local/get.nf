process GET_DUPLICATIONS {
    tag "$sampleId"
    cpus   6
    memory '10 GB'

    publishDir "${params.outdir}/DUPLICATIONS/${sampleId}",
           mode: 'copy',
           pattern: "*"                  

    input:
    tuple val(sampleId), path(amap_tsv), path(amap_faa), path(amap_fna), path(bakta_tsv), path(fasta), path(emapper)

    output:
    path("*.dup"), emit: dup_fastas

    script:
    """
    python3 ${projectDir}/bin/rename_fastas.py ${emapper} ${amap_faa}
    diamond blastp \\
        --db ${amap_faa}.renamed \\
        -q ${amap_faa}.renamed \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen nident \\
        --very-sensitive \\
        --no-self-hits \\
        --header \\
        > duplications_${sampleId}

    python3 ${projectDir}/bin/get_duplicationsv2.py duplications_${sampleId}
    python3 ${projectDir}/bin/dup_ex_multifasta.py duplications_${sampleId}_groups.tsv ${sampleId}.AMAP.faa.renamed
    """
}

process GET_OUTPUTS {
    tag "$sampleId"

    /* copy every file whose name contains AMAP */
    publishDir "${params.outdir}/AMAP/${sampleId}",
               mode: 'copy',
               pattern: "*AMAP*"

    input:
    tuple val(sampleId),
          path(amap_tsv),
          path(amap_faa),
          path(amap_fna),
          path(bakta_tsv),
          path(fasta),
          path(emapper)

    output:
    tuple val(sampleId), path("*-AMAPmastertable.tsv"), emit: master_table
    path "*AMAP*"                                           

    script:
    """
    set -euo pipefail

    # remove header once
    if ! grep -q '^lineage' ${amap_tsv}; then
        sed -i '1d' ${amap_tsv}
    fi

    python3 ${projectDir}/bin/generate_tsv_gff_gbk.py \\
        ${emapper} \\
        ${amap_tsv} \\
        ${bakta_tsv} \\
        "${params.mycobrowser_ref}" \\
        ${fasta}
    python3 ${projectDir}/bin/to_panaroo.py --id ${sampleId}

    cat ${amap_faa} > ${sampleId}-AMAP.faa
    cat ${amap_fna} > ${sampleId}-AMAP.fna
    """
}

process GET_COUNT {

    tag "$sampleId"

    /* Copy every summary TSV produced by count_annot.py */
    publishDir "${params.outdir}/AMAP/${sampleId}",
               mode: 'copy',
               pattern: "*-AMAP.*.tsv"

    input:
    tuple val(sampleId),
          path(mastertable)

    output:
    path("*.tsv")

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/count_annot.py ${mastertable}
    """
}

process GET_ALLDUP {
    cpus 13
    memory '10 GB'
    publishDir "${params.outdir}/DUPLICATIONS", mode: 'copy'

    input:
    path dup_fastas

    output:
    tuple path("homologous_communities.tsv"), path("homologous_communities.html")

    script:
    """
    cat *.dup > all_fastas.fa
    rm *.dup

    diamond blastp \\
        --db all_fastas.fa \\
        -q all_fastas.fa \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen nident \\
        --very-sensitive \\
        --no-self-hits \\
        --header \\
        > all_vs_all.diamond

    python3 ${projectDir}/bin/network_genomes.py all_vs_all.diamond
    """
}
