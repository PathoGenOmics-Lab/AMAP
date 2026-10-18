process PATHOTYPR {
    tag "$sampleId"

    publishDir "${params.outdir}/PATHOTYPR/typing", pattern: "*_summary.tsv", mode: 'link', overwrite: true
    publishDir "${params.outdir}/PATHOTYPR/fastas", pattern: "*.fasta", mode: 'link', overwrite: true

    input:
    tuple val(sampleId), path(fasta)

    output:
    tuple env('SAMPLE_NAME'), path("*.fasta"), path("*_summary.tsv"), emit: typed_fasta
    // The name, the sample ID and the lineage, for the steps that write the
    // last two down: a sample ID can hold underscores, so they cannot be read
    // back out of the name.
    tuple env('SAMPLE_NAME'), val(sampleId), env('LINEAGE'), emit: names

    script:
    """
    #!/bin/bash
    set -o pipefail
    # NOTE: do NOT use 'set -e' here. If the lineage call fails,
    # we need to reach the fallback code that sets SAMPLE_NAME.
    # With set -e, the script exits before SAMPLE_NAME is exported,
    # leaving env() empty → downstream sampleId is blank.

    # 1. Run the lineage classifier (may fail for non-MTBC genomes). The
    # subcommand is quoted: Nextflow adds to a task's cache key any bin/ entry
    # its script names as a bare word, and bare, this one named AMAP's own
    # bin/classify package, which this task never runs.
    pathotypr 'classify' \\
        -i ${fasta} \\
        -m "${params.pathotyper_snps}" \\
        -r "${params.pathotyper_ref}" \\
        --output-prefix "${sampleId}" || true

    # 2. Ensure summary file exists (even if empty) — Nextflow output glob requires it.
    if [ ! -f "${sampleId}_summary.tsv" ]; then
        echo -e "Sample\\tLineage\\tName" > "${sampleId}_summary.tsv"
    fi

    # 3. Check if the summary file was created and is not empty.
    if [ -s "${sampleId}_summary.tsv" ]; then
        # 3. Extract the new name: <lineage>_<sampleId>
        newName=\$(awk -F'\\t' 'NR==2 {print \$NF"_"\$1}' "${sampleId}_summary.tsv")
        lineage=\$(awk -F'\\t' 'NR==2 {print \$NF}' "${sampleId}_summary.tsv")

        if [ -n "\$newName" ]; then
            # Rename FASTA header (first contig only)
            awk -v name="\$newName" 'BEGIN{first=1} /^>/{if(first){print ">"name; first=0}else{print}; next} {print}' "${fasta}" > "\${newName}.fasta"
            # Export as env var for Nextflow (no stdout pollution)
            export SAMPLE_NAME="\$newName"
            export LINEAGE="\${lineage:-NA}"
        else
            echo "Error: Could not determine new name from summary file. Using original sampleId." >&2
            awk -v name="NA_${sampleId}" 'BEGIN{first=1} /^>/{if(first){print ">"name; first=0}else{print}; next} {print}' "${fasta}" > "NA_${sampleId}.fasta"
            export SAMPLE_NAME="NA_${sampleId}"
            export LINEAGE=NA
        fi
    else
        echo "Error: Summary file not found or empty. Using original sampleId." >&2
        awk -v name="NA_${sampleId}" 'BEGIN{first=1} /^>/{if(first){print ">"name; first=0}else{print}; next} {print}' "${fasta}" > "NA_${sampleId}.fasta"
        export SAMPLE_NAME="NA_${sampleId}"
        export LINEAGE=NA
    fi
    """
}
