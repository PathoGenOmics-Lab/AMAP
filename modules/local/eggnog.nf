process EGGNOG_DB {
    tag 'eggnog_db_setup'
    cpus   16
    memory '45 GB'

    output:
    path "eggnog_db", emit: db_path

    script:
    def eggnog_db_dir = params.eggnog_db.startsWith('/') \
        ? params.eggnog_db \
        : "${projectDir}/${params.eggnog_db}"

    """
    set -euo pipefail
    mkdir -p ${eggnog_db_dir}

    # ------------------------------------------------------------
    # 0) Quick exit if the bacteria DB already exists
    # ------------------------------------------------------------
    if [ -f "${eggnog_db_dir}/bacteria.dmnd" ]; then
        echo "==> bacteria.dmnd found – skipping download and build."
        ln -sfn ${eggnog_db_dir} eggnog_db
        exit 0
    fi

    # ------------------------------------------------------------
    # 1) Download metadata only (skip 60-GB Diamond core)
    # ------------------------------------------------------------
    echo "==> Downloading metadata only (skipping Diamond core)…"
    download_eggnog_data.py -D -y --data_dir ${eggnog_db_dir}

    # ------------------------------------------------------------
    # 2) Build the bacteria-only Diamond index (~3 GB)
    # ------------------------------------------------------------
    echo "==> Building bacteria.dmnd (taxid 2)…"
    create_dbs.py -m diamond \
                  --dbname bacteria \
                  --taxa Bacteria \
                  -y \
                  --data_dir ${eggnog_db_dir}

    # ------------------------------------------------------------
    # 3) Clean up build-only files
    # ------------------------------------------------------------
    echo "==> Removing build-only FASTA and taxid table…"
    rm -f ${eggnog_db_dir}/e5.proteomes.faa  || true
    rm -f ${eggnog_db_dir}/e5.taxid_info.tsv || true
    # (Keep eggnog.taxa.db.traverse.pkl; delete it if you need the extra 150 MB.)

    # Symlink for downstream processes
    ln -sfn ${eggnog_db_dir} eggnog_db
    """
}

process ANNOT_EGGNOGMAP {
  tag "$sampleId"
  cpus 24
  memory '50 GB'

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

