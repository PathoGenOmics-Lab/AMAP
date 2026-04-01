process GET_DUPLICATIONS {
    tag "$sampleId"

    publishDir "${params.outdir}/DUPLICATIONS/${sampleId}",
           mode: 'link', overwrite: true,
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
    diamond makedb --in ${amap_faa}.renamed --db renamed_db --threads ${task.cpus}
    diamond blastp \\
        --db renamed_db \\
        -q ${amap_faa}.renamed \\
        --threads ${task.cpus} \\
        --tmpdir \${TMPDIR:-/tmp} \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen nident \\
        --very-sensitive \\
        --no-self-hits \\
        --header \\
        > duplications_${sampleId}.raw

    # Sort for reproducibility (Diamond multi-thread output order is non-deterministic)
    grep '^#' duplications_${sampleId}.raw > duplications_${sampleId} || true
    grep -v '^#' duplications_${sampleId}.raw | sort -k1,1 -k2,2 -k12,12gr >> duplications_${sampleId}
    rm duplications_${sampleId}.raw

    python3 ${projectDir}/bin/get_duplications.py duplications_${sampleId}
    python3 ${projectDir}/bin/dup_ex_multifasta.py duplications_${sampleId}_groups.tsv ${amap_faa}.renamed
    """
}

process GET_OUTPUTS {
    tag "$sampleId"

    publishDir "${params.outdir}/AMAP/${sampleId}",
               mode: 'link', overwrite: true,
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

    # Remove header line if present: check if line 1 starts with known header
    # keywords. Avoid heuristic on column content which can misfire.
    if head -1 local_amap.tsv | grep -qiE '^(id|gene_id|locus|#)'; then
        tail -n +2 local_amap.tsv > local_amap_tmp.tsv && mv local_amap_tmp.tsv local_amap.tsv
    fi

    # Determine Mycobrowser argument at runtime (file may not exist outside container)
    MYCO_REF="${params.mycobrowser_ref}"
    if [ "${params.skip_mycobrowser}" = "true" ]; then
        MYCO_ARG="NONE"
    elif [ -f "\$MYCO_REF" ]; then
        MYCO_ARG="\$MYCO_REF"
    else
        echo "WARNING: Mycobrowser reference not found at \$MYCO_REF — skipping Mycobrowser" >&2
        MYCO_ARG="NONE"
    fi

    python3 ${projectDir}/bin/generate_tsv_gff_gbk.py \\
        ${emapper} \\
        local_amap.tsv \\
        ${bakta_tsv} \\
        "\$MYCO_ARG" \\
        ${fasta}
    python3 ${projectDir}/bin/to_panaroo.py --id ${sampleId}

    # Hard link (same inode, no extra I/O) — avoids broken symlinks with publishDir
    ln -f ${amap_faa} ${sampleId}-AMAP.faa || cp ${amap_faa} ${sampleId}-AMAP.faa
    ln -f ${amap_fna} ${sampleId}-AMAP.fna || cp ${amap_fna} ${sampleId}-AMAP.fna
    """
}

process GET_COUNT {
    tag "$sampleId"

    publishDir "${params.outdir}/AMAP/${sampleId}",
               mode: 'link', overwrite: true,
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
    publishDir "${params.outdir}/DUPLICATIONS", mode: 'link', overwrite: true

    input:
    path dup_fastas

    output:
    tuple path("homologous_communities.tsv"), path("homologous_communities.html")

    script:
    """
    # Concatenate all per-sample duplication FASTAs
    # Use find+xargs instead of glob to avoid ARG_MAX with 500+ files
    find . -maxdepth 1 -name '*.dup' -print0 | xargs -0 cat > all_fastas.fa
    find . -maxdepth 1 -name '*.dup' -delete

    # Handle case where no samples have duplicated genes (empty concatenated FASTA)
    if [ ! -s all_fastas.fa ]; then
        echo "No duplicated genes across any sample — creating empty outputs."
        echo -e "" > all_vs_all.diamond
        echo -e "" > homologous_communities.tsv
        # Minimal valid HTML for the network visualization
        echo "<html><body><p>No homologous communities found</p></body></html>" > homologous_communities.html
        exit 0
    fi

    # Scale Diamond block size with available RAM for faster processing
    MEM_GB=\$(echo "${task.memory}" | grep -oP '\\d+' | head -1)
    BLOCK_SIZE=\$(( MEM_GB / 6 ))
    [ \$BLOCK_SIZE -lt 2 ] && BLOCK_SIZE=2
    [ \$BLOCK_SIZE -gt 8 ] && BLOCK_SIZE=8

    diamond makedb --in all_fastas.fa --db all_fastas_db --threads ${task.cpus}
    diamond blastp \\
        --db all_fastas_db \\
        -q all_fastas.fa \\
        --threads ${task.cpus} \\
        --tmpdir \${TMPDIR:-/tmp} \\
        --block-size \$BLOCK_SIZE \\
        --index-chunks 1 \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen nident \\
        --very-sensitive \\
        --no-self-hits \\
        --header \\
        > all_vs_all.raw

    # Sort for reproducibility (Diamond multi-thread output order is non-deterministic)
    grep '^#' all_vs_all.raw > all_vs_all.diamond || true
    grep -v '^#' all_vs_all.raw | sort -k1,1 -k2,2 -k12,12gr >> all_vs_all.diamond
    rm all_vs_all.raw

    python3 ${projectDir}/bin/network_genomes.py all_vs_all.diamond
    """
}
