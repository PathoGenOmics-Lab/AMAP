/**
 * Capture software versions for reproducibility.
 *
 * Logs the exact version of every tool used in the pipeline.
 * Output is a YAML file suitable for inclusion in publications
 * and a TSV for the HTML report.
 */
process SOFTWARE_VERSIONS {
    tag 'versions'
    executor 'local'
    cpus 1
    memory '256.MB'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    output:
    path 'software_versions.yml', emit: versions_yml
    path 'software_versions.tsv', emit: versions_tsv

    script:
    """
    cat <<-YAML > software_versions.yml
    pipeline:
      name: "${workflow.manifest.name ?: 'AMAP'}"
      version: "${workflow.manifest.version ?: 'dev'}"
      nextflow: "${nextflow.version}"
      nextflow_build: "${nextflow.build}"
    tools:
    YAML

    # Collect tool versions (each with fallback if not installed)
    bakta_ver=\$(bakta --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' || echo 'N/A')
    miniprot_ver=\$(miniprot --version 2>&1 || echo 'N/A')
    diamond_ver=\$(diamond version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' || echo 'N/A')
    emapper_ver=\$(emapper.py --version 2>&1 | head -1 | grep -oP '\\d+\\.\\d+\\.\\d+' || echo 'N/A')
    pangene_ver=\$(pangene --version 2>&1 || echo 'N/A')
    python_ver=\$(python3 --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+')
    pandas_ver=\$(python3 -c "import pandas; print(pandas.__version__)" 2>/dev/null || echo 'N/A')
    biopython_ver=\$(python3 -c "import Bio; print(Bio.__version__)" 2>/dev/null || echo 'N/A')
    intervaltree_ver=\$(python3 -c "import intervaltree; print(intervaltree.__version__)" 2>/dev/null || echo 'N/A')
    r_ver=\$(R --version 2>&1 | head -1 | grep -oP '\\d+\\.\\d+\\.\\d+' || echo 'N/A')
    pathotypr_ver=\$(pathotypr --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' || echo 'N/A')

    cat <<-YAML >> software_versions.yml
      bakta: "\$bakta_ver"
      miniprot: "\$miniprot_ver"
      diamond: "\$diamond_ver"
      eggnog_mapper: "\$emapper_ver"
      pangene: "\$pangene_ver"
      pathotypr: "\$pathotypr_ver"
      python: "\$python_ver"
      pandas: "\$pandas_ver"
      biopython: "\$biopython_ver"
      intervaltree: "\$intervaltree_ver"
      r: "\$r_ver"
    YAML

    # TSV for easy parsing
    echo -e "Tool\\tVersion" > software_versions.tsv
    echo -e "AMAP\\t${workflow.manifest.version ?: 'dev'}" >> software_versions.tsv
    echo -e "Nextflow\\t${nextflow.version}" >> software_versions.tsv
    echo -e "Bakta\\t\$bakta_ver" >> software_versions.tsv
    echo -e "Miniprot\\t\$miniprot_ver" >> software_versions.tsv
    echo -e "Diamond\\t\$diamond_ver" >> software_versions.tsv
    echo -e "eggNOG-mapper\\t\$emapper_ver" >> software_versions.tsv
    echo -e "Pangene\\t\$pangene_ver" >> software_versions.tsv
    echo -e "Pathotypr\\t\$pathotypr_ver" >> software_versions.tsv
    echo -e "Python\\t\$python_ver" >> software_versions.tsv
    echo -e "pandas\\t\$pandas_ver" >> software_versions.tsv
    echo -e "Biopython\\t\$biopython_ver" >> software_versions.tsv
    echo -e "intervaltree\\t\$intervaltree_ver" >> software_versions.tsv
    echo -e "R\\t\$r_ver" >> software_versions.tsv
    """
}

/**
 * Dump all pipeline parameters to a YAML file for full reproducibility.
 * This captures the EXACT configuration used for a given run.
 */
process DUMP_PARAMS {
    tag 'params'
    executor 'local'
    cpus 1
    memory '128.MB'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    output:
    path 'params.yml'

    script:
    // Write all params as YAML
    def paramsList = params.collect { k, v -> "  ${k}: \"${v}\"" }.join('\n')
    """
    cat <<-YAML > params.yml
    # AMAP Pipeline Parameters
    # Run: ${new java.util.Date()}
    # Nextflow: ${nextflow.version}
    # Pipeline: ${workflow.manifest.version ?: 'dev'}
    # Command: ${workflow.commandLine}
    # Profile: ${workflow.profile}
    # Work dir: ${workflow.workDir}
    # Launch dir: ${workflow.launchDir}
    params:
    ${paramsList}
    YAML
    """
}
