process EGGNOG_DB {
    tag 'eggnog_db_setup'
    cpus   16
    memory '45 GB'

    output:
    path "eggnog_db", emit: db_path

    script:
    def eggnog_db_dir = params.eggnog_db.startsWith('/')
        ? params.eggnog_db
        : "${projectDir}/${params.eggnog_db}"

    """
    set -euo pipefail

    mkdir -p ${eggnog_db_dir}

    echo "==> Downloading metadata only (skipping Diamond core)…"
    download_eggnog_data.py -D -y --data_dir ${eggnog_db_dir}

    # ------------------------------------------------------------------
    # Build bacteria-only Diamond DB once (≈3 GB)
    # ------------------------------------------------------------------
    if [ ! -f "${eggnog_db_dir}/diamond/bacteria.dmnd" ]; then
        echo "==> Building bacteria.dmnd (taxid 2)…"
        create_dbs.py -m diamond \\
                      --dbname bacteria \\
                      --taxa Bacteria \\
                      -y \\
                      --data_dir ${eggnog_db_dir}
    else
        echo "==> bacteria.dmnd already present – skipping build."
    fi

    # ------------------------------------------------------------------
    # Clean up files no longer needed for annotation
    # ------------------------------------------------------------------
    echo "==> Removing build-only FASTA and taxid table…"
    rm -f ${eggnog_db_dir}/e5.proteomes.faa  || true
    rm -f ${eggnog_db_dir}/e5.taxid_info.tsv || true


    # Symlink for downstream processes
    ln -sfn ${eggnog_db_dir} eggnog_db
    """
}
