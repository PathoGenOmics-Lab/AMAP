process EGGNOG_DB {
    tag 'eggnog_db_setup'

    output:
    path "eggnog_db", emit: db_path

    script:
    """
    set -euo pipefail

    # Resolve the user-provided eggnog_db path
    EGGNOG_DB_DIR="${params.eggnog_db}"
    if [[ ! "\$EGGNOG_DB_DIR" = /* ]]; then
        EGGNOG_DB_DIR="${projectDir}/${params.eggnog_db}"
    fi

    mkdir -p "\$EGGNOG_DB_DIR"

    # ------------------------------------------------------------
    # 0) Quick exit if the bacteria DB already exists
    # ------------------------------------------------------------
    if [ -f "\$EGGNOG_DB_DIR/bacteria.dmnd" ]; then
        echo "==> bacteria.dmnd found – skipping download and build."
        # Use relative symlink when possible for portability with -resume
        if command -v realpath >/dev/null 2>&1; then
            ln -sfn "\$(realpath --relative-to=. "\$EGGNOG_DB_DIR")" eggnog_db
        else
            ln -sfn "\$EGGNOG_DB_DIR" eggnog_db
        fi
        exit 0
    fi

    # ------------------------------------------------------------
    # 1) Download metadata only (skip 60-GB Diamond core)
    # ------------------------------------------------------------
    echo "==> Downloading metadata only (skipping Diamond core)…"
    download_eggnog_data.py -D -y --data_dir "\$EGGNOG_DB_DIR"

    # ------------------------------------------------------------
    # 2) Build the bacteria-only Diamond index (~3 GB)
    # ------------------------------------------------------------
    echo "==> Building bacteria.dmnd (taxid 2)…"
    create_dbs.py -m diamond \
                  --dbname bacteria \
                  --taxa Bacteria \
                  -y \
                  --data_dir "\$EGGNOG_DB_DIR"

    # ------------------------------------------------------------
    # 3) Clean up build-only files
    # ------------------------------------------------------------
    echo "==> Removing build-only FASTA and taxid table…"
    rm -f "\$EGGNOG_DB_DIR"/e5.proteomes.faa  || true
    rm -f "\$EGGNOG_DB_DIR"/e5.taxid_info.tsv || true

    # Symlink into workdir for Nextflow output channel
    if command -v realpath >/dev/null 2>&1; then
        ln -sfn "\$(realpath --relative-to=. "\$EGGNOG_DB_DIR")" eggnog_db
    else
        ln -sfn "\$EGGNOG_DB_DIR" eggnog_db
    fi
    """
}

/*
 * BATCH eggNOG annotation — runs eggNOG-mapper ONCE on all samples combined.
 *
 * WHY: On HPC clusters with shared filesystems (NFS/Lustre), running N
 * parallel eggNOG-mapper jobs that each load the ~3 GB Diamond database
 * causes I/O contention. With 200+ concurrent jobs, the filesystem becomes
 * the bottleneck and annotation is slower than serial.
 *
 * HOW: Concatenate all per-sample FAA files with sample-prefixed headers,
 * run eggNOG-mapper once with all CPUs, then split the output back by sample.
 * The DB is loaded into memory once (--dbmem), and Diamond uses all cores.
 */
process ANNOT_EGGNOGMAP_BATCH {
  tag "batch_${faa_files.size()}_samples"

  input:
  path faa_files    // collected FAA files from all samples
  val  sample_ids   // list of sampleIds in matching order
  path eggnog_db

  output:
  path "split/*.emapper.annotations", emit: split_annotations

  script:
  """
  set -euo pipefail
  mkdir -p tmp scr_tmp split

  # 1. Concatenate all FAAs with sample-prefixed headers
  #    Original header: >L8_G26689_Rv0001;+;0-1524
  #    Tagged header:   >SAMPLE=L8_G26689|||L8_G26689_Rv0001;+;0-1524
  for faa in ${faa_files}; do
      # Extract sampleId from filename: {sampleId}.AMAP.faa
      sid=\$(basename "\$faa" .AMAP.faa)
      awk -v sid="\$sid" '/^>/{sub(/^>/, ""); print ">SAMPLE=" sid "|||" \$0; next} {print}' "\$faa"
  done > all_proteins.faa

  # 2. Run eggNOG-mapper ONCE on the combined file
  emapper.py \\
    -i all_proteins.faa \\
    -m diamond \\
    --dmnd_db eggnog_db/bacteria.dmnd \\
    --cpu ${task.cpus} \\
    --usemem \\
    --dbmem \\
    --sensmode fast \\
    --tax_scope 2 \\
    --target_orthologs one2one \\
    -o batch \\
    --data_dir eggnog_db \\
    --temp_dir tmp \\
    --scratch_dir scr_tmp

  # 3. Split results back into per-sample annotation files
  #    - Strip the SAMPLE=xxx||| prefix from the query column
  #    - Route each line to the correct sample file
  head -n 5 batch.emapper.annotations | grep '^#' > header.txt || true

  # Create per-sample files with headers
  for faa in ${faa_files}; do
      sid=\$(basename "\$faa" .AMAP.faa)
      cp header.txt "split/\${sid}.emapper.annotations"
  done

  # Split annotation lines by sample prefix
  grep -v '^#' batch.emapper.annotations | while IFS=\$'\\t' read -r query rest; do
      # Extract sample ID from SAMPLE=xxx|||original_query
      sid=\$(echo "\$query" | sed 's/^SAMPLE=//; s/|||.*//')
      # Restore original query name
      orig_query=\$(echo "\$query" | sed 's/^SAMPLE=[^|]*|||//')
      echo -e "\${orig_query}\\t\${rest}" >> "split/\${sid}.emapper.annotations"
  done
  """
}

/*
 * Reshape batch output back into per-sample tuples for downstream joins.
 */
process SPLIT_EGGNOGMAP {
  tag "$sampleId"

  input:
  tuple val(sampleId), path(annotation)

  output:
  tuple val(sampleId), path("${sampleId}.emapper.annotations")

  script:
  """
  # Just rename/copy — the file is already split by ANNOT_EGGNOGMAP_BATCH
  cp ${annotation} ${sampleId}.emapper.annotations
  """
}

