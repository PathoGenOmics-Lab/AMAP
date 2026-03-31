process PATHOTYPR {
    tag "$sampleId"
    cpus 1
    memory '1 GB'

    publishDir "${params.outdir}/PATHOTYPR/typing", pattern: "*_summary.tsv", mode: 'copy'
    publishDir "${params.outdir}/PATHOTYPR/fastas", pattern: "*.fasta", mode: 'copy'

    input:
    tuple val(sampleId), path(fasta)

    output:
    tuple stdout, path("*.fasta"), path("*_summary.tsv"), emit: typed_fasta

    script:
    """
    #!/bin/bash
    set -e
    set -o pipefail

    # 1. Run the lineage classifier.
    # This command generates an output file named `${sampleId}_summary.tsv`.
    pathotypr classify \\
        -i ${fasta} \\
        -m "${params.pathotyper_snps}" \\
        -r "${params.pathotyper_ref}" \\
        --output-prefix "${sampleId}"

    # 2. Check if the summary file was created and is not empty.
    if [ -s "${sampleId}_summary.tsv" ]; then
        # 3. Extract the new name from the summary file.
        #    - `awk` is used to read the second line (NR==2) of the TSV file (-F'\\t').
        #    - A new name is constructed with the format: <last_column>_<first_column>
        #      Example: L9_G48909
        newName=\$(awk -F'\\t' 'NR==2 {print \$NF"_"\$1}' "${sampleId}_summary.tsv")

        # 4. Check if the new name was extracted successfully.
        if [ -n "\$newName" ]; then
            # 5. Modify the FASTA header and save to a new file.
            #    - `sed` replaces the line starting with '>' with the new name.
            # Only rename the FIRST header; preserve other contig headers in multi-contig FASTAs
            awk -v name="\$newName" 'BEGIN{first=1} /^>/{if(first){print ">"name; first=0}else{print}; next} {print}' "${fasta}" > "\${newName}.fasta"

            # 6. Print the new name to standard output for Nextflow to capture.
            echo -n "\$newName"
        else
            echo "Error: Could not determine new name from summary file. Using original sampleId." >&2
            awk -v name="NA_${sampleId}" 'BEGIN{first=1} /^>/{if(first){print ">"name; first=0}else{print}; next} {print}' "${fasta}" > "NA_${sampleId}.fasta"
            echo -n "${sampleId}"
        fi
    else
        echo "Error: Summary file not found or empty. Using original sampleId." >&2
        awk -v name="NA_${sampleId}" 'BEGIN{first=1} /^>/{if(first){print ">"name; first=0}else{print}; next} {print}' "${fasta}" > "NA_${sampleId}.fasta"
        echo -n "${sampleId}"
    fi
    """
}
