/**
 * PREPARE_REFERENCE
 *
 * Extracts protein sequences from a user-provided reference genome + annotation
 * (GFF3 or GenBank). Produces the proteome FAA needed for Miniprot alignment
 * and the reference genome FASTA for reverse alignment.
 *
 * When params.ref_annotation is set, this replaces the hardcoded H37Rv files.
 * When not set, the pipeline uses the bundled H37Rv references (default).
 */
process PREPARE_REFERENCE {
    tag 'prepare_ref'

    input:
    path ref_annotation
    path ref_genome

    output:
    path "custom_ref_ref_proteins.faa", emit: ref_proteome
    path ref_genome, emit: ref_genome_fasta

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/extract_ref_proteins.py \
        --ref-annotation ${ref_annotation} \
        --ref-genome ${ref_genome} \
        --prefix custom_ref
    """
}
