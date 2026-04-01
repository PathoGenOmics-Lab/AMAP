#!/usr/bin/env nextflow
/*
========================================================================================
    A M A P   P I P E L I N E
========================================================================================
    Assisted MTBC Annotator Pipeline
----------------------------------------------------------------------------------------
    Publication:
    DOI: https://doi.org/10.5281/zenodo.8024446
----------------------------------------------------------------------------------------
    Usage:
    nextflow run PathoGenOmics-Lab/AMAP -profile singularity --tsv <your_samples.tsv>
----------------------------------------------------------------------------------------
*/

// Enable Nextflow DSL 2
nextflow.enable.dsl=2

// --- MANIFEST ---
manifest {
    name            = 'PathoGenOmics-Lab/AMAP'
    author          = 'Paula Ruiz-Rodriguez'
    description     = 'Assisted Microbial Annotator Pipeline'
    mainScript      = 'main.nf'
    version         = '1.0.0'
    nextflowVersion = '>=23.04.0'
    doi             = 'https://doi.org/10.5281/zenodo.8024446'
    homePage        = 'https://github.com/PathoGenOmics-Lab/AMAP'
    defaultBranch   = 'main'
}

//
// --- PARAMS ---
//
params.tsv = 'example_files/test.tsv'
params.outdir = 'results'
params.eggnog_db = 'eggnog_db'
params.help = false
params.validate_only = false
params.skip_typing = false
params.skip_functional_annotation = false
params.genus = "Mycobacterium"
params.species = "Mycobacterium tuberculosis"
params.bakta_db = "/data/bakta_db/db-light"
params.skip_mycobrowser = false

// Custom reference (set both to use your own reference instead of H37Rv)
params.ref_annotation = null   // GFF3 or GenBank of a trusted reference
params.ref_genome = null       // FASTA of the reference genome

// Default H37Rv reference files (used when ref_annotation/ref_genome are not set)
params.pathotyper_snps = "/reference/snps_lineages.tsv"
params.pathotyper_ref = "/reference/MTB_ancestor_reference.fasta"
params.miniprot_ref_aa = "/reference/aa_H37Rv.translated.fas"
params.miniprot_ref_bad_bakta = "/reference/bad.bakta.translated.fas"
params.miniprot_ref_h37rv = "/reference/new_refH37Rv.fasta"
params.bakta_proteins = "/reference/sequence.gb"
params.mycobrowser_ref = "/reference/Mycobacterium_tuberculosis_H37Rv_txt_v4.txt"

//
// --- LOG AND HELP ---
//
def logo = """\
=======================================
       A M A P   P I P E L I N E
=======================================
           _,--',   _._.--._____
   .--.--';_'-.', ";_         _.,-'
   .'--'. _.'     {`'-;_ .-.>.'
           '-:_       )  / `' '=.
             ) >     {_/,     /~)
             |/               `^ .'
=======================================
`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`
    `=,´=/     `=,´=/     `=,´=/
      y==/         y==/       y==/
    ,=,-<=`.     ,=,-<=`.    ,=,-<=`.
 ,-'-'   `-=_, -'-'   `-=_, -'-'   `-=
"""

if (params.help) {
    log.info logo
    log.info """
    Usage:
    nextflow run ${workflow.manifest.name} -profile <profile> [options]

    Required Arguments:
    --tsv                       Path to the input TSV file.
                                Format: sampleId<TAB>fasta_path

    Optional Arguments:
    --outdir                    Path to the output directory. (Default: "${params.outdir}")
    --skip_typing               Skip the initial PATHOTYPR lineage typing step. (Default: ${params.skip_typing})
    --skip_functional_annotation Skip all functional annotation and downstream analysis steps. (Default: ${params.skip_functional_annotation})
    --genus                     Genus for Bakta annotation. (Default: "${params.genus}")
    --species                   Species for Bakta annotation. (Default: "${params.species}")
    --bakta_db                  Path to the Bakta database directory. (Default: "${params.bakta_db}")
    --skip_mycobrowser          Skip Mycobrowser functional categories (for non-MTBC organisms). (Default: ${params.skip_mycobrowser})
    --validate_only             Run input validation only (no pipeline execution).
    --max_cpus                  Maximum CPUs per process. (Default: ${params.max_cpus})
    --max_memory                Maximum memory per process. (Default: ${params.max_memory})
    --max_time                  Maximum time per process. (Default: ${params.max_time})
    --scratch                   Use node-local scratch for I/O-heavy processes. (Default: ${params.scratch})
    --help                      Display this help message.
    --version                   Print pipeline version and exit.

    Organism-Agnostic Mode:
    --ref_annotation            Path to a trusted reference annotation (GFF3 or GenBank).
    --ref_genome                Path to the trusted reference genome (FASTA).
                                Both must be set together. When used, the pipeline extracts
                                proteins from your reference and uses them instead of H37Rv.
    --skip_mycobrowser          Auto-enabled when genus is not Mycobacterium.
    --skip_typing               Required for non-MTBC organisms (Pathotypr is MTBC-specific).
    
    Example — annotating E. coli:
      nextflow run PathoGenOmics-Lab/AMAP -profile docker \\
        --tsv samples.tsv --genus Escherichia --species "Escherichia coli" \\
        --ref_annotation ref.gbk --ref_genome ref.fasta --skip_typing

    Reference & Database Arguments:
    --eggnog_db                 Path to the pre-formatted eggNOG database directory. (Default: "${params.eggnog_db}")
    --pathotyper_snps           Path to the SNPs TSV file for Pathotyper. (Default: "${params.pathotyper_snps}")
    --pathotyper_ref            Path to the reference FASTA file for Pathotyper. (Default: "${params.pathotyper_ref}")
    --miniprot_ref_aa           Path to the H37Rv amino acid reference for Miniprot. (Default: "${params.miniprot_ref_aa}")
    --miniprot_ref_bad_bakta    Path to the 'bad bakta' reference FASTA for Miniprot. (Default: "${params.miniprot_ref_bad_bakta}")
    --miniprot_ref_h37rv        Path to the H37Rv reference genome FASTA for Miniprot. (Default: "${params.miniprot_ref_h37rv}")
    --bakta_proteins            Path to the reference proteins GBK file for Bakta. (Default: "${params.bakta_proteins}")
    --mycobrowser_ref           Path to the Mycobrowser reference TXT file. (Default: "${params.mycobrowser_ref}")
    """
    exit 0
}

// --- VERSION FLAG ---
if (args.containsKey('version') || args.containsKey('v')) {
    log.info "${manifest.name} v${manifest.version}"
    log.info "Nextflow ${nextflow.version}"
    System.exit(0)
}

log.info logo
log.info """
    ┌─────────────────────────────────────────────────┐
    │  Pipeline    : ${manifest.name} v${manifest.version}
    │  Nextflow    : ${nextflow.version}
    │  Profile     : ${workflow.profile ?: 'standard'}
    ├─────────────────────────────────────────────────┤
    │  Input TSV   : ${params.tsv}
    │  Output      : ${params.outdir}
    │  Reference   : ${params.ref_annotation ? 'Custom (' + params.ref_annotation + ')' : 'H37Rv (built-in)'}
    │  Organism    : ${params.genus} ${params.species}
    ├─────────────────────────────────────────────────┤
    │  Typing      : ${params.skip_typing ? '⏭  skipped' : '✓  Pathotypr'}
    │  Annotation  : ${params.skip_functional_annotation ? '⏭  skipped' : '✓  eggNOG + Bakta'}
    │  Mycobrowser : ${params.skip_mycobrowser ? '⏭  skipped' : '✓  enabled'}
    │  Validate    : ${params.validate_only ? '🔍 only (dry run)' : '✓  pre-flight'}
    ├─────────────────────────────────────────────────┤
    │  Max CPUs    : ${params.max_cpus}
    │  Max Memory  : ${params.max_memory}
    │  Max Time    : ${params.max_time}
    │  Scratch     : ${params.scratch ? '✓  node-local' : '✗  disabled'}
    └─────────────────────────────────────────────────┘
"""

// --- COMPLETION HANDLER ---
workflow.onComplete {
    def duration = workflow.duration
    def succeeded = workflow.stats.succeededCount
    def failed = workflow.stats.failedCount
    def ignored = workflow.stats.ignoredCount
    def status = workflow.success ? "✓ COMPLETED" : "✗ FAILED"

    log.info """
    =======================================
    ${status}
    =======================================
    Pipeline    : ${manifest.name} v${manifest.version}
    Nextflow    : ${nextflow.version}
    Profile     : ${workflow.profile}
    Duration    : ${duration}
    Succeeded   : ${succeeded}
    Failed      : ${failed}
    Ignored     : ${ignored}
    Output      : ${params.outdir}
    Work dir    : ${workflow.workDir}
    Command     : ${workflow.commandLine}
    =======================================
    """.stripIndent()

    if (failed > 0 && ignored > 0) {
        log.warn """
    ${failed} process(es) failed, ${ignored} ignored (non-critical samples).
    
    To investigate:
      1. Check execution report: ${params.outdir}/pipeline_info/execution_report_*.html
      2. Check trace file:       ${params.outdir}/pipeline_info/execution_trace_*.txt
      3. Full log:               .nextflow.log
      
    To retry failed samples only:
      nextflow run PathoGenOmics-Lab/AMAP -resume [same params]
        """.stripIndent()
    } else if (failed > 0) {
        log.error """
    ${failed} process(es) failed.
    
    To investigate:
      1. Check execution report: ${params.outdir}/pipeline_info/execution_report_*.html
      2. Full log:               .nextflow.log
      
    To retry:
      nextflow run PathoGenOmics-Lab/AMAP -resume [same params]
        """.stripIndent()
    }
    
    // Write machine-readable run summary
    def summary_file = new File("${params.outdir}/pipeline_info/run_summary.txt")
    summary_file.parentFile.mkdirs()
    summary_file.text = """\
    pipeline=${manifest.name}
    version=${manifest.version}
    nextflow=${nextflow.version}
    profile=${workflow.profile}
    status=${workflow.success ? 'COMPLETED' : 'FAILED'}
    duration=${duration}
    succeeded=${succeeded}
    failed=${failed}
    ignored=${ignored}
    start=${workflow.start}
    complete=${workflow.complete}
    command=${workflow.commandLine}
    workDir=${workflow.workDir}
    launchDir=${workflow.launchDir}
    projectDir=${workflow.projectDir}
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline failed: ${workflow.errorMessage}"
    log.error "Check .nextflow.log for full stack trace."
}

// Validate custom reference params
if ((params.ref_annotation && !params.ref_genome) || (!params.ref_annotation && params.ref_genome)) {
    error "Both --ref_annotation and --ref_genome must be provided together. You set one but not the other."
}

// Auto-detect organism-agnostic mode:
// When using a custom reference AND genus is not Mycobacterium, auto-enable
// --skip_mycobrowser to avoid referencing MTBC-specific categories.
if (params.ref_annotation && params.genus != 'Mycobacterium' && !params.skip_mycobrowser) {
    log.info "Auto-enabling --skip_mycobrowser (genus '${params.genus}' is not Mycobacterium)"
    params.skip_mycobrowser = true
}

// Suggest --skip_typing when genus is not Mycobacterium
if (params.genus != 'Mycobacterium' && !params.skip_typing) {
    log.warn "Genus '${params.genus}' is not Mycobacterium but --skip_typing is not set. " +
             "Pathotypr lineage typing is MTBC-specific and will likely fail for other organisms. " +
             "Consider adding --skip_typing to your command."
}

// validate_only handled in workflow block below
//
// --- MODULE INCLUDES ---
//
include { PATHOTYPR }           from './modules/local/pathotypr'
include { BAKTA_ANNOTATOR }     from './modules/local/bakta'
include { MINIPROT }            from './modules/local/miniprot'
include { TYPE_CASES }          from './modules/local/cases'
include { FIX_CDS }             from './modules/local/cases'
include { EGGNOG_DB }               from './modules/local/eggnog'
include { ANNOT_EGGNOGMAP_BATCH }   from './modules/local/eggnog'
include { GET_COUNT }           from './modules/local/get'
include { GET_DUPLICATIONS }    from './modules/local/get'
include { GET_OUTPUTS }         from './modules/local/get'
include { GET_ALLDUP }          from './modules/local/get'
include { PANGENE }             from './modules/local/pangene'
include { REPORT_INPUTS }       from './modules/local/freport'
include { FINAL_REPORT }        from './modules/local/freport'
include { PREPARE_REFERENCE }   from './modules/local/reference'
include { VALIDATE_INPUTS }     from './modules/local/validation'
include { SOFTWARE_VERSIONS }   from './modules/local/versions'
include { DUMP_PARAMS }         from './modules/local/versions'
//
// --- INPUT CHANNEL ---
//
Channel
    .fromPath(params.tsv)
    .ifEmpty { error """
    ╔══════════════════════════════════════════════════════════╗
    ║  Input TSV file not found: ${params.tsv}
    ║
    ║  Make sure the path is correct and the file exists.
    ║  Use --tsv to specify the input file.
    ╚══════════════════════════════════════════════════════════╝
    """.stripIndent() }
    .splitCsv(header: true, sep: '\t')
    .map { row ->
        if (!row.containsKey('sampleId') || !row.containsKey('fasta')) {
            error """
    ╔══════════════════════════════════════════════════════════╗
    ║  Invalid TSV format!
    ║
    ║  Expected columns: sampleId  fasta
    ║  Found columns:    ${row.keySet().join('  ')}
    ║
    ║  Your TSV should look like:
    ║    sampleId\\tfasta
    ║    sample1\\t/path/to/sample1.fasta
    ║    sample2\\t/path/to/sample2.fasta
    ╚══════════════════════════════════════════════════════════╝
    """.stripIndent()
        }
        def f = file(row.fasta)
        if (!f.exists()) {
            error """
    ╔══════════════════════════════════════════════════════════╗
    ║  FASTA file not found!
    ║
    ║  Sample:  ${row.sampleId}
    ║  Path:    ${row.fasta}
    ║
    ║  Check that the path is correct (absolute or relative
    ║  to where you launched nextflow).
    ╚══════════════════════════════════════════════════════════╝
    """.stripIndent()
        }
        tuple(row.sampleId, f)
    }
    .set { fasta_ch }

// Count samples for the user
fasta_ch.count().subscribe { n ->
    log.info "    ✓ Loaded ${n} sample(s) from ${params.tsv}"
    if (n == 0) log.warn "    ⚠ No samples found! Check your TSV file."
}

//
// --- WORKFLOW ---
//
workflow {
    // 0. Input validation + reproducibility — always runs first
    VALIDATE_INPUTS(Channel.value(file(params.tsv)))
    SOFTWARE_VERSIONS()
    DUMP_PARAMS()

    // If validate_only mode, stop here
    if (params.validate_only) {
        return
    }

    // 1. Setup eggNOG database if not present
    if (!params.skip_functional_annotation) {
        EGGNOG_DB()
    }

    // 2. Perform MTB lineage typing or skip it
    if (!params.skip_typing) {
        PATHOTYPR(fasta_ch)
        typed_ch = PATHOTYPR.out.typed_fasta
    } else {
        // If skipping typing, create a placeholder file for caching.
        CREATE_EMPTY_FILE()
        typed_ch = fasta_ch.combine(CREATE_EMPTY_FILE.out).map { sampleId, fasta, empty_file ->
            tuple(sampleId, fasta, empty_file)
        }
    }

    // 3. Annotate genomes with Bakta
    BAKTA_ANNOTATOR(typed_ch.combine(Channel.value(params.genus)).combine(Channel.value(params.species)))

    // 4. Prepare reference proteome (custom or default H37Rv)
    if (params.ref_annotation && params.ref_genome) {
        log.info "Using custom reference: ${params.ref_annotation} + ${params.ref_genome}"
        PREPARE_REFERENCE(
            Channel.value(file(params.ref_annotation)),
            Channel.value(file(params.ref_genome))
        )
        ref_proteome_ch     = PREPARE_REFERENCE.out.ref_proteome
        ref_genome_fasta_ch = PREPARE_REFERENCE.out.ref_genome_fasta
        bad_bakta_ch        = PREPARE_REFERENCE.out.bad_bakta
    } else {
        ref_proteome_ch     = Channel.value(file(params.miniprot_ref_aa))
        ref_genome_fasta_ch = Channel.value(file(params.miniprot_ref_h37rv))
        bad_bakta_ch        = Channel.value(file(params.miniprot_ref_bad_bakta))
    }

    // 5. Align proteins to the genome with Miniprot
    MINIPROT(BAKTA_ANNOTATOR.out.for_miniprot, ref_proteome_ch, ref_genome_fasta_ch, bad_bakta_ch)

    // 6. Classify CDS cases based on alignments
    TYPE_CASES(MINIPROT.out.paf_files.join(BAKTA_ANNOTATOR.out.for_cases))

    // 7. Correct CDS based on classified cases
    FIX_CDS(TYPE_CASES.out.for_fix.join(BAKTA_ANNOTATOR.out.for_fix))

    // --- FUNCTIONAL ANNOTATION AND DOWNSTREAM ANALYSIS ---
    // This entire block will be skipped if --skip_functional_annotation is used
    if (!params.skip_functional_annotation) {
        // 7. Annotate corrected proteins with eggNOG-mapper (BATCH mode)
        //    Collects all FAA files, runs eggNOG ONCE, splits results back.
        //    This avoids N parallel jobs competing for the same Diamond DB on
        //    shared filesystems — critical for HPC with NFS/Lustre.
        eggnog_faa_ch = FIX_CDS.out.for_eggnog.map { it[1] }.collect()
        eggnog_ids_ch = FIX_CDS.out.for_eggnog.map { it[0] }.collect()

        ANNOT_EGGNOGMAP_BATCH(eggnog_faa_ch, eggnog_ids_ch, EGGNOG_DB.out.db_path)

        // Flatten batch output back to per-sample tuples: [sampleId, annotation]
        eggnog_results_ch = ANNOT_EGGNOGMAP_BATCH.out.split_annotations
            .flatten()
            .map { f ->
                def sid = f.baseName.replace('.emapper', '')
                tuple(sid, f)
            }

        // Combine results for downstream processing
        ch_for_downstream = FIX_CDS.out.for_downstream.join(eggnog_results_ch)

        // 8. Generate final output files (GFF, GBK, Master Table)
        GET_OUTPUTS(ch_for_downstream)

        // 9. Find and analyze gene duplications
        GET_DUPLICATIONS(ch_for_downstream)
        
        // 10. Runs pangenome for the genes of each genome
        // Collect typed FASTA files from PATHOTYPR step (or empty channel if skipped)
        if (!params.skip_typing) {
            ref_fastas_ch = PATHOTYPR.out.typed_fasta.map { it[1] }.collect()
        } else {
            // Use a collected empty list instead of Channel.empty() so that
            // PANGENE still executes (with 0 reference FASTAs → creates
            // placeholder outputs). Channel.empty() prevents process execution entirely.
            ref_fastas_ch = Channel.value([])
        }
        PANGENE(ch_for_downstream, ref_fastas_ch)

        // 11. Count annotations per category
        GET_COUNT(GET_OUTPUTS.out.master_table)

        // 12. Generate final reports
        // Collect summaries from different processes
        summaries_cases_ch = TYPE_CASES.out.summaries.collect()
        summaries_annot_ch = GET_COUNT.out.collect()

        // 13. Create report inputs
        // Collect rvmissing files and master tables for comprehensive report
        rvmissing_ch = TYPE_CASES.out.rvmissing.collect()
        mastertables_ch = GET_OUTPUTS.out.master_table.map { it[1] }.collect()

        // Collect pangene Rtab files (or empty if skipped)
        pangene_rtab_ch = PANGENE.out.pangene_out.map { it[3] }.collect().ifEmpty( file('NO_PANGENE') )
        qc_stats_ch = GET_OUTPUTS.out.qc_stats.collect()

        REPORT_INPUTS(summaries_annot_ch, summaries_cases_ch, rvmissing_ch, mastertables_ch, pangene_rtab_ch, qc_stats_ch)

        // 14. Analyze duplications across all samples
        GET_ALLDUP(GET_DUPLICATIONS.out.dup_fastas.collect())

        // 15. Render the final HTML report
        FINAL_REPORT(REPORT_INPUTS.out, GET_ALLDUP.out)
    }
}

//
// --- PROCESSES ---
//
/**
 * Creates a single, stable empty file to be used as a placeholder.
 * This ensures that downstream processes that need this file as an
 * input can be cached correctly when using -resume.
 */
process CREATE_EMPTY_FILE {
    // This process will be cached after its first successful run.
    output:
    path 'empty_summary.log'

    script:
    "touch empty_summary.log"
}
