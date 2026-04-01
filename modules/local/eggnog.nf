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
  #    Delimiter: @@SMPL@@ (safe — no protein ID contains this)
  #    Original:  >L8_G26689_Rv0001;+;0-1524
  #    Tagged:    >SAMPLE@@SMPL@@L8_G26689@@SMPL@@L8_G26689_Rv0001;+;0-1524
  for faa in ${faa_files}; do
      sid=\$(basename "\$faa" .AMAP.faa)
      awk -v sid="\$sid" '/^>/{sub(/^>/, ""); print ">SAMPLE@@SMPL@@" sid "@@SMPL@@" \$0; next} {print}' "\$faa"
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
  #    All in one awk pass — O(N) instead of O(N×M) subprocesses
  # Capture ALL header lines (any number of #-prefixed lines)
  grep '^#' batch.emapper.annotations > header.txt || true

  # Create per-sample files with headers
  for faa in ${faa_files}; do
      sid=\$(basename "\$faa" .AMAP.faa)
      cp header.txt "split/\${sid}.emapper.annotations"
  done

  # Split annotation lines by sample — single awk pass over the file
  # Handles empty results gracefully (grep -v returns 1 if no lines → || true)
  grep -v '^#' batch.emapper.annotations 2>/dev/null | awk -F'\\t' '
  {
      # Parse query: SAMPLE@@SMPL@@sampleId@@SMPL@@original_query
      query = \$1
      # Extract sampleId (between first and second @@SMPL@@)
      split(query, parts, "@@SMPL@@")
      sid = parts[2]
      orig_query = parts[3]
      # Reconstruct line with original query
      \$1 = orig_query
      print >> ("split/" sid ".emapper.annotations")
  }
  ' OFS='\\t' || true
  """
}



