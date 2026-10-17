/**
 * PREPARE_REFERENCE
 *
 * When the user provides --ref_annotation + --ref_genome:
 * 1. Extracts protein sequences from the curated annotation → ref proteome
 * 2. Annotates the SAME genome with Bakta
 * 3. Compares Bakta vs curated → finds genes Bakta invents (false positives)
 * 4. Outputs: ref proteome, ref genome, and bad_bakta FAA
 *
 * This replaces the hardcoded H37Rv files and the manually-curated
 * bad.bakta.translated.fas.
 */
process PREPARE_REFERENCE {
    tag 'prepare_ref'
    // Resources defined in nextflow.config (PREPARE_REFERENCE)

    input:
    path ref_annotation
    path ref_genome

    output:
    path "custom_ref_ref_proteins.faa", emit: ref_proteome
    path ref_genome, emit: ref_genome_fasta
    path "custom_ref_bad_bakta.faa", emit: bad_bakta

    script:
    """
    set -euo pipefail

    # 1. Extract curated proteins from the user's trusted annotation
    python3 ${projectDir}/bin/extract_ref_proteins.py \
        --ref-annotation ${ref_annotation} \
        --ref-genome ${ref_genome} \
        --prefix custom_ref

    # 2. Annotate the same reference genome with Bakta
    #    (to find what Bakta invents on a known genome)
    bakta \
        --db "${params.bakta_db}" \
        --skip-trna \
        --threads ${task.cpus} \
        --genus "${params.genus}" \
        --species "${params.species}" \
        --output bakta_ref \
        --prefix bakta_ref \
        --locus-tag bakta_ref \
        ${ref_genome}

    # 3. Compare Bakta annotation vs curated → extract false positives
    python3 ${projectDir}/bin/find_bakta_false_positives.py \
        --bakta-tsv bakta_ref/bakta_ref.tsv \
        --curated-faa custom_ref_ref_proteins.faa \
        --bakta-faa bakta_ref/bakta_ref.faa \
        --prefix custom_ref

    echo "==> Reference preparation complete"
    echo "    Curated proteins: \$(grep -c '^>' custom_ref_ref_proteins.faa)"
    echo "    Bakta false positives: \$(grep -c '^>' custom_ref_bad_bakta.faa || echo 0)"
    """
}

/**
 * BUNDLED_REFERENCE
 *
 * The default H37Rv reference lives inside the container, under /reference.
 * Handing those paths to a process as `path` inputs makes Nextflow bind-mount
 * /reference from the host, where it does not exist: Docker Desktop refuses
 * the mount, and Docker on Linux creates an empty directory and mounts it over
 * the container's own, hiding every file in it. Singularity skips the missing
 * path, which is why this only ever showed with Docker. Copying the files out
 * from inside the container gives the rest of the pipeline ordinary files to
 * stage, the same way TEST_GENOME does for the genome.
 */
process BUNDLED_REFERENCE {
    tag 'h37rv'

    output:
    path "ref_proteome.faa", emit: ref_proteome
    path "ref_genome.fasta", emit: ref_genome_fasta
    path "bad_bakta.faa", emit: bad_bakta

    script:
    """
    copy() {
        if [ ! -s "\$1" ]; then
            echo "ERROR: bundled reference file not found: \$1" >&2
            echo "       It ships in the AMAP container; outside it, pass your own with" >&2
            echo "       --ref_annotation and --ref_genome." >&2
            exit 1
        fi
        cp "\$1" "\$2"
    }
    copy "${params.miniprot_ref_aa}" ref_proteome.faa
    copy "${params.miniprot_ref_h37rv}" ref_genome.fasta
    copy "${params.miniprot_ref_bad_bakta}" bad_bakta.faa
    """
}
