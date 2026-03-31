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
        ln -sfn "\$EGGNOG_DB_DIR" eggnog_db
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
    ln -sfn "\$EGGNOG_DB_DIR" eggnog_db
    """
}

process ANNOT_EGGNOGMAP {
  tag "$sampleId"

  input:
  tuple val(sampleId), path(fixfaa)
  path eggnog_db

  output:
  tuple val(sampleId), path("${sampleId}.emapper.annotations")

  script:
  """
  set -euo pipefail
  mkdir -p tmp
  mkdir -p scr_tmp

  emapper.py \
    -i ${fixfaa} \
    -m diamond \
    --dmnd_db eggnog_db/bacteria.dmnd \
    --cpu ${task.cpus} \
    --usemem \
    --dbmem \
    --sensmode fast \
    --tax_scope 2 \
    --target_orthologs one2one \
    -o ${sampleId} \
    --data_dir eggnog_db \
    --temp_dir tmp \
    --scratch_dir scr_tmp
  """
}

