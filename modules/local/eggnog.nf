process EGGNOG_DB {
    tag 'eggnog_db_setup'

    // storeDir caches the DB outside the work directory.
    // Once built, it's reused across ALL pipeline runs — never rebuilt.
    storeDir "${params.eggnog_db}"

    output:
    path "eggnog_db", emit: db_path

    script:
    """
    set -euo pipefail

    # Resolve the user-provided eggnog_db path
    # Use launchDir (where user ran nextflow) for relative paths, NOT projectDir
    # (which points to the pipeline cache and is wrong for user-supplied paths)
    EGGNOG_DB_DIR="${params.eggnog_db}"
    if [[ ! "\$EGGNOG_DB_DIR" = /* ]]; then
        EGGNOG_DB_DIR="${launchDir}/${params.eggnog_db}"
    fi

    mkdir -p "\$EGGNOG_DB_DIR"

    # ------------------------------------------------------------
    # 0) Quick exit if the bacteria DB already exists
    # ------------------------------------------------------------
    if [ -f "\$EGGNOG_DB_DIR/bacteria.dmnd" ]; then
        # Integrity check: Diamond DB should be >100 MB
        DB_SIZE=\$(stat -f%z "\$EGGNOG_DB_DIR/bacteria.dmnd" 2>/dev/null || stat -c%s "\$EGGNOG_DB_DIR/bacteria.dmnd" 2>/dev/null || echo 0)
        if [ "\$DB_SIZE" -lt 100000000 ]; then
            echo "WARNING: bacteria.dmnd is only \${DB_SIZE} bytes (expected >100 MB) — may be corrupted, rebuilding..." >&2
            rm -f "\$EGGNOG_DB_DIR/bacteria.dmnd"
        else
            echo "==> bacteria.dmnd found (\${DB_SIZE} bytes) – skipping download and build."
            # Use relative symlink when possible for portability with -resume
            if command -v realpath >/dev/null 2>&1; then
                ln -sfn "\$(realpath --relative-to=. "\$EGGNOG_DB_DIR")" eggnog_db
            else
                ln -sfn "\$EGGNOG_DB_DIR" eggnog_db
            fi
            exit 0
        fi
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

  echo "▸ [eggNOG] Preparing batch annotation for ${faa_files.size()} samples..." >&2

  # 1. Concatenate all FAAs with sample-prefixed headers
  #    Delimiter: @@SMPL@@ (safe — no protein ID contains this)
  #    Original:  >L8_G26689_Rv0001;+;0-1524
  #    Tagged:    >SAMPLE@@SMPL@@L8_G26689@@SMPL@@L8_G26689_Rv0001;+;0-1524
  for faa in ${faa_files}; do
      sid=\$(basename "\$faa" .AMAP.faa)
      awk -v sid="\$sid" '/^>/{sub(/^>/, ""); print ">SAMPLE@@SMPL@@" sid "@@SMPL@@" \$0; next} {print}' "\$faa"
  done > all_proteins.faa

  # 2. Run eggNOG-mapper ONCE on the combined file
  # Check if there's anything to annotate (all samples might have 0 CDS)
  N_PROTS=\$(grep -c '^>' all_proteins.faa || echo 0)
  if [ "\$N_PROTS" -eq 0 ]; then
      echo "WARNING: No proteins to annotate (all samples have 0 CDS)" >&2
      # Create minimal output so split step works
      echo "# emapper: no input proteins" > batch.emapper.annotations
  else
      echo "==> Annotating \$N_PROTS proteins from ${faa_files.size()} samples..."

      # Optimize Diamond block size based on available memory.
      # Larger blocks = fewer disk passes = faster. Rule: ~6 GB RAM per block.
      MEM_GB=\$(echo "${task.memory}" | grep -oP '\\d+' | head -1)
      BLOCK_SIZE=\$(( MEM_GB / 8 ))
      [ \$BLOCK_SIZE -lt 2 ] && BLOCK_SIZE=2
      [ \$BLOCK_SIZE -gt 12 ] && BLOCK_SIZE=12

      emapper.py \\
        -i all_proteins.faa \\
        -m diamond \\
        --dmnd_db eggnog_db/bacteria.dmnd \\
        --cpu ${task.cpus} \\
        --usemem \\
        --dbmem \\
        --sensmode fast \\
        --dmnd_block_size \$BLOCK_SIZE \\
        --dmnd_index_chunks 1 \\
        --tax_scope 2 \\
        --target_orthologs one2one \\
        -o batch \\
        --data_dir eggnog_db \\
        --temp_dir tmp \\
        --scratch_dir scr_tmp
  fi

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
  # NOTE: close() after each print >> to avoid hitting ulimit -n (open files)
  # on large runs (500+ samples). Without close(), awk keeps all output files
  # open simultaneously which fails at ~256-1024 samples on most systems.
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
      outfile = "split/" sid ".emapper.annotations"
      print >> outfile
      close(outfile)
  }
  ' OFS='\\t' || true
  """
}



