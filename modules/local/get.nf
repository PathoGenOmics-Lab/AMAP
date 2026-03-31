process GET_DUPLICATIONS {
    tag "$sampleId"

    publishDir "${params.outdir}/DUPLICATIONS/${sampleId}",
           mode: 'copy',
           pattern: "duplications*"                  

    input:
    tuple val(sampleId), path(amap_tsv), path(amap_faa), path(amap_fna), path(bakta_tsv), path(fasta), path(emapper)

    output:
    path("*.dup"), emit: dup_fastas
    path("duplications*")

    script:
    """
    python3 ${projectDir}/bin/rename_fastas.py ${emapper} ${amap_faa}

    # Build Diamond DB first, then query against it
    diamond makedb --in ${amap_faa}.renamed --db renamed_db
    diamond blastp \\
        --db renamed_db \\
        -q ${amap_faa}.renamed \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen nident \\
        --very-sensitive \\
        --no-self-hits \\
        --header \\
        > duplications_${sampleId}

    python3 ${projectDir}/bin/get_duplications.py duplications_${sampleId}
    python3 ${projectDir}/bin/dup_ex_multifasta.py duplications_${sampleId}_groups.tsv ${amap_faa}.renamed
    """
}

process GET_OUTPUTS {
    tag "$sampleId"

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

    # Copy input TSV to avoid modifying Nextflow-staged symlink (breaks -resume)
    cp ${amap_tsv} local_amap.tsv

    # Remove header line if present: detect by checking if the first field
    # of line 1 is NOT a valid coordinate (i.e. not purely numeric after the ID)
    if head -1 local_amap.tsv | awk -F'\\t' '{exit (\$3 ~ /^[0-9]+\$/) ? 1 : 0}'; then
        tail -n +2 local_amap.tsv > local_amap_tmp.tsv && mv local_amap_tmp.tsv local_amap.tsv
    fi

    if [ "${params.skip_mycobrowser}" = "true" ] || [ "${params.skip_mycobrowser}" = "false" -a ! -f "${params.mycobrowser_ref}" ]; then
        MYCO_ARG="NONE"
    else
        MYCO_ARG="${params.mycobrowser_ref}"
    fi

    python3 ${projectDir}/bin/generate_tsv_gff_gbk.py \\
        ${emapper} \\
        local_amap.tsv \\
        ${bakta_tsv} \\
        "\$MYCO_ARG" \\
        ${fasta}
    python3 ${projectDir}/bin/to_panaroo.py --id ${sampleId}

    cat ${amap_faa} > ${sampleId}-AMAP.faa
    cat ${amap_fna} > ${sampleId}-AMAP.fna
    """
}

process GET_COUNT {
    tag "$sampleId"

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
    publishDir "${params.outdir}/DUPLICATIONS", mode: 'copy'

    input:
    path dup_fastas

    output:
    tuple path("homologous_communities.tsv"), path("homologous_communities.html")

    script:
    """
    cat *.dup > all_fastas.fa
    rm *.dup

    diamond makedb --in all_fastas.fa --db all_fastas_db
    diamond blastp \\
        --db all_fastas_db \\
        -q all_fastas.fa \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen nident \\
        --very-sensitive \\
        --no-self-hits \\
        --header \\
        > all_vs_all.diamond

    python3 ${projectDir}/bin/network_genomes.py all_vs_all.diamond
    """
}
