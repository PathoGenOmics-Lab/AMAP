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

// Manifest is defined in nextflow.config

//
// --- PARAMS ---
//
params.tsv = 'samples.tsv'
params.outdir = 'results'
params.eggnog_db = 'eggnog_db'
params.help = false
params.version = false
params.validate_only = false
// Validation warnings stop the run too, not only its errors
params.strict_validation = false
params.skip_typing = false
params.skip_functional_annotation = false
params.genus = "Mycobacterium"
params.species = "Mycobacterium tuberculosis"
params.bakta_db = "/data/bakta_db/db-light"
// Set by the `test` profile: take the genome from the container, not a TSV.
params.test_data = false

// Custom reference (set both to use your own reference instead of H37Rv)
params.ref_annotation = null   // GFF3 or GenBank of a trusted reference
params.ref_genome = null       // FASTA of the reference genome

// Mycobrowser categories describe MTBC genes only, so a custom reference from
// another genus skips them unless --skip_mycobrowser is given either way. This
// is decided here, in the default, because Nextflow ignores any later
// assignment to a param.
params.skip_mycobrowser = params.ref_annotation && params.genus != 'Mycobacterium'

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
// Nextflow 25.10 and later parse scripts with the strict syntax by default,
// which allows no statements outside a workflow or a function. Everything
// below that used to run at the top of this file is therefore a function, and
// the entry workflow calls them in the same order.
def logo() {
    return """\
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
}

def helpMessage() {
    return """
    Usage:
    nextflow run ${workflow.manifest.name} -profile <profile> [options]

    Required Arguments:
    --tsv                       Path to the input TSV file.
                                Format: sampleId<TAB>fasta_path
                                Not needed with -profile test, which annotates
                                the reference genome bundled in the container.

    Optional Arguments:
    --outdir                    Path to the output directory. (Default: "${params.outdir}")
    --skip_typing               Skip the initial PATHOTYPR lineage typing step. (Default: ${params.skip_typing})
    --skip_functional_annotation Skip all functional annotation and downstream analysis steps. (Default: ${params.skip_functional_annotation})
    --genus                     Genus for Bakta annotation. (Default: "${params.genus}")
    --species                   Species for Bakta annotation. (Default: "${params.species}")
    --bakta_db                  Path to the Bakta database directory. (Default: "${params.bakta_db}")
    --skip_mycobrowser          Skip Mycobrowser functional categories (for non-MTBC organisms). (Default: ${params.skip_mycobrowser})
    --validate_only             Run input validation only (no pipeline execution).
    --strict_validation         Stop on validation warnings too, not only on errors. (Default: ${params.strict_validation})
    --max_cpus                  Maximum CPUs per process. (Default: ${params.max_cpus})
    --max_memory                Maximum memory per process. (Default: ${params.max_memory})
    --max_time                  Maximum time per process. (Default: ${params.max_time})
    --scratch                   Use node-local scratch for I/O-heavy processes. (Default: ${params.scratch})
    --help                      Display this help message.
    --version                   Print pipeline and Nextflow versions and exit.

    Organism-Agnostic Mode:
    --ref_annotation            Path to a trusted reference annotation (GFF3 or GenBank).
    --ref_genome                Path to the trusted reference genome (FASTA).
                                Both must be set together. When used, the pipeline extracts
                                proteins from your reference and uses them instead of H37Rv.
    --skip_mycobrowser          Auto-enabled when a custom reference is set and
                                genus is not Mycobacterium. Set it yourself
                                otherwise.
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
}

def runSummary() {
    return """
    ┌─────────────────────────────────────────────────┐
    │  Pipeline    : ${workflow.manifest.name} v${workflow.manifest.version}
    │  Nextflow    : ${nextflow.version}
    │  Profile     : ${workflow.profile ?: 'standard'}
    ├─────────────────────────────────────────────────┤
    │  Input       : ${params.test_data ? 'bundled test genome' : params.tsv}
    │  Output      : ${params.outdir}
    │  Reference   : ${params.ref_annotation ? 'Custom (' + params.ref_annotation + ')' : 'H37Rv (built-in)'}
    │  Organism    : ${params.genus} ${params.species}
    ├─────────────────────────────────────────────────┤
    │  Typing      : ${params.skip_typing ? '⏭  skipped' : '✓  Pathotypr'}
    │  Annotation  : ${params.skip_functional_annotation ? '⏭  skipped' : '✓  eggNOG + Bakta'}
    │  Mycobrowser : ${params.skip_mycobrowser ? '⏭  skipped' : '✓  enabled'}
    │  Validate    : ${params.validate_only ? '🔍 only (dry run)' : '✓  pre-flight'}${params.strict_validation ? ', strict' : ''}
    ├─────────────────────────────────────────────────┤
    │  Max CPUs    : ${params.max_cpus}
    │  Max Memory  : ${params.max_memory}
    │  Max Time    : ${params.max_time}
    │  Scratch     : ${params.scratch ? '✓  node-local' : '✗  disabled'}
    └─────────────────────────────────────────────────┘
"""
}

// --- COMPLETION HANDLER ---
def onCompleteHandler() {
    def duration = workflow.duration
    def succeeded = workflow.stats.succeededCount
    def failed = workflow.stats.failedCount
    def ignored = workflow.stats.ignoredCount
    def status = workflow.success ? "✓ COMPLETED" : "✗ FAILED"

    log.info """
    =======================================
    ${status}
    =======================================
    Pipeline    : ${workflow.manifest.name} v${workflow.manifest.version}
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
    pipeline=${workflow.manifest.name}
    version=${workflow.manifest.version}
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

def onErrorHandler() {
    log.error "Pipeline failed: ${workflow.errorMessage}"
    log.error "Check .nextflow.log for full stack trace."
}

def checkParams() {
    // Validate custom reference params
    if ((params.ref_annotation && !params.ref_genome) || (!params.ref_annotation && params.ref_genome)) {
        error "Both --ref_annotation and --ref_genome must be provided together. You set one but not the other."
    }

    // Organism-agnostic mode: params.skip_mycobrowser is decided with the
    // other params at the top of this file; this only says why.
    if (params.ref_annotation && params.genus != 'Mycobacterium' && params.skip_mycobrowser) {
        log.info "Skipping Mycobrowser categories (genus '${params.genus}' is not Mycobacterium)"
    }

    // Suggest --skip_typing when genus is not Mycobacterium
    if (params.genus != 'Mycobacterium' && !params.skip_typing) {
        log.warn "Genus '${params.genus}' is not Mycobacterium but --skip_typing is not set. " +
                 "Pathotypr lineage typing is MTBC-specific and will likely fail for other organisms. " +
                 "Consider adding --skip_typing to your command."
    }
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
include { PANGENE_GENOME }      from './modules/local/pangene'
include { REPORT_INPUTS }       from './modules/local/freport'
include { FINAL_REPORT }        from './modules/local/freport'
include { AMAP_REPORT }         from './modules/local/amap_report'
include { PREPARE_REFERENCE; BUNDLED_REFERENCE } from './modules/local/reference'
include { VALIDATE_INPUTS }     from './modules/local/validation'
include { SOFTWARE_VERSIONS }   from './modules/local/versions'
include { DUMP_PARAMS }         from './modules/local/versions'
//
// --- INPUT CHANNEL ---
//
// The test profile does not read a sample sheet. Its genome ships inside the
// container, and the check below runs in the Nextflow head process on the host,
// before any container exists, so a /reference path in a TSV would be reported
// as missing. TEST_GENOME copies it out from inside the container instead. That
// keeps a 4 MB genome out of this repository and off the host entirely.
def readSampleSheet() {
    // Checked here rather than with ifEmpty on the channel: fromPath throws on a
    // missing literal path before ifEmpty ever runs, so the message below was
    // never the one anybody saw.
    if (!file(params.tsv).exists()) {
        error """
        ╔══════════════════════════════════════════════════════════╗
        ║  Input TSV file not found: ${params.tsv}
        ║
        ║  Make sure the path is correct and the file exists.
        ║  Use --tsv to specify the input file, or try the pipeline
        ║  with no data of your own:
        ║
        ║    nextflow run PathoGenOmics-Lab/AMAP -profile docker,test
        ╚══════════════════════════════════════════════════════════╝
        """.stripIndent()
    }

    def fasta_ch = Channel
        .fromPath(params.tsv)
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

    // Count samples for the user
    fasta_ch.count().subscribe { n ->
        log.info "    ✓ Loaded ${n} sample(s) from ${params.tsv}"
        if (n == 0) log.warn "    ⚠ No samples found! Check your TSV file."
    }
    return fasta_ch
}

//
// --- REFERENCE FILES ---
//
// A reference file a task reads by path, as the input that task gets: the file,
// staged like any other input so that a container sees it, when it is on the
// host; nothing otherwise, and the task reads the path it was given inside the
// container, where the defaults live. Handed over as a path, a default would
// make Nextflow mount the host's missing /reference over the container's own
// (see BUNDLED_REFERENCE).
def stagedIfOnHost(path) {
    return (path && file(path).exists()) ? channel.value(file(path)) : channel.value([])
}

// The reference files this run reads, each with the param that gives it
def referenceFiles() {
    def files = params.ref_annotation
        ? [['ref_annotation', params.ref_annotation], ['ref_genome', params.ref_genome]]
        : [['bakta_proteins', params.bakta_proteins],
           ['miniprot_ref_aa', params.miniprot_ref_aa],
           ['miniprot_ref_h37rv', params.miniprot_ref_h37rv],
           ['miniprot_ref_bad_bakta', params.miniprot_ref_bad_bakta]]
    if (!params.skip_typing) {
        files += [['pathotyper_snps', params.pathotyper_snps], ['pathotyper_ref', params.pathotyper_ref]]
    }
    // Read by AMAP_REPORT in every run, and by GET_OUTPUTS with functional
    // annotation, so it is checked with or without it
    if (!params.skip_mycobrowser) {
        files += [['mycobrowser_ref', params.mycobrowser_ref]]
    }
    return files
}

//
// --- WORKFLOW ---
//
workflow {
    // What used to run at the top of this file, in the same order: help, the
    // run summary, the handlers, the parameter checks and the sample sheet.
    if (params.help) {
        log.info logo()
        log.info helpMessage()
        exit 0
    }
    if (params.version) {
        // println rather than log.info, which Nextflow 26 prints with no line
        // break between messages and `nextflow -q` does not print at all.
        println "Pipeline    : ${workflow.manifest.name} v${workflow.manifest.version}"
        println "Nextflow    : ${nextflow.version}"
        exit 0
    }
    log.info logo()
    log.info runSummary()
    workflow.onComplete = { onCompleteHandler() }
    workflow.onError = { onErrorHandler() }
    checkParams()
    if (params.test_data) {
        log.info "    ✓ Test profile: using the reference genome bundled in the container"
        fasta_ch = channel.empty()
    } else {
        fasta_ch = readSampleSheet()
    }

    // The databases. One on the host is staged into every task that reads it,
    // like any other input, so that a container sees it; otherwise the tasks
    // get nothing and look for its path inside the container, where the Bakta
    // default lives. Handed over as a path, that default would make Nextflow
    // mount the host's missing /data/bakta_db over the container's own (see
    // BUNDLED_REFERENCE).
    bakta_db_ch = (params.bakta_db && file(params.bakta_db).exists()) ? channel.value(file(params.bakta_db)) : channel.value([])
    // eggNOG's is checked only when it is used, by the path EGGNOG_DB builds it
    // at (a relative one is relative to launchDir), and staged if it is there
    // already: it need not be.
    def eggnog_dir = (!params.skip_functional_annotation && params.eggnog_db) ? file(params.eggnog_db) : null
    eggnog_db_check_ch = channel.value(eggnog_dir ? tuple(eggnog_dir.toString(), eggnog_dir.exists() ? eggnog_dir : []) : tuple('', []))

    // The reference files the tasks read, each staged into them when it is on
    // the host (see stagedIfOnHost). Bakta's proteins are the custom
    // reference's annotation when there is one.
    bakta_proteins_ch  = stagedIfOnHost(params.ref_annotation ?: params.bakta_proteins)
    pathotyper_snps_ch = stagedIfOnHost(params.pathotyper_snps)
    pathotyper_ref_ch  = stagedIfOnHost(params.pathotyper_ref)
    mycobrowser_ref_ch = stagedIfOnHost(params.mycobrowser_ref)
    // and those not on the host, which the tasks read inside the container,
    // for VALIDATE_INPUTS to check that they are there
    refs_in_container_ch = channel.value(referenceFiles().findAll { _param, given -> !(given && file(given).exists()) })

    // 0. Input validation and reproducibility, before anything else. The test
    // profile has no sample sheet to validate, only the reference genome the
    // container ships, so validating params.tsv would fail on a file that was
    // never meant to exist. The genomes go in with the sheet, as the pipeline
    // will annotate them: a container sees only the inputs Nextflow stages.
    if (!params.test_data) {
        VALIDATE_INPUTS(
            Channel.value(file(params.tsv)),
            fasta_ch.toList().map { rows -> tuple(rows.collect { row -> row[0] }, rows.collect { row -> row[1] }) },
            bakta_db_ch,
            eggnog_db_check_ch,
            refs_in_container_ch
        )
    }
    SOFTWARE_VERSIONS()
    DUMP_PARAMS()

    // With --validate_only, stop after the checks above. Nextflow 26.10 and
    // later refuse a return in a workflow body, and exit() would cancel those
    // checks too, so the rest of the pipeline sits inside this if instead.
    if (!params.validate_only) {
        // 1. Setup eggNOG database if not present. Its folder is made here, on
        // the host, and handed to EGGNOG_DB as an input, so that its container
        // sees the folder (see EGGNOG_DB).
        if (!params.skip_functional_annotation) {
            eggnog_dir?.mkdirs()
            if (!eggnog_dir?.isDirectory()) {
                error "Cannot make a folder for the eggNOG database at --eggnog_db '${params.eggnog_db}'"
            }
            EGGNOG_DB(channel.value(eggnog_dir))
        }

        // 1b. Where the genomes come from: the sample sheet, or the container.
        genomes_ch = params.test_data ? TEST_GENOME() : fasta_ch

        // 2. Perform MTB lineage typing or skip it. names_ch says, for every
        // genome, the name its files carry, its sample ID and its lineage, for
        // the steps that write the last two down: typing names a genome
        // <lineage>_<sampleId>, and a sample ID can hold underscores of its
        // own, so neither can be read back out of the name.
        if (!params.skip_typing) {
            PATHOTYPR(genomes_ch, pathotyper_snps_ch, pathotyper_ref_ch)
            typed_ch = PATHOTYPR.out.typed_fasta
            names_ch = PATHOTYPR.out.names
        } else {
            // If skipping typing, create a placeholder file for caching.
            CREATE_EMPTY_FILE()
            typed_ch = genomes_ch.combine(CREATE_EMPTY_FILE.out).map { sampleId, fasta, empty_file ->
                tuple(sampleId, fasta, empty_file)
            }
            names_ch = genomes_ch.map { sampleId, _fasta -> tuple(sampleId, sampleId, 'Untyped') }
        }

        // 3. Annotate genomes with Bakta
        BAKTA_ANNOTATOR(typed_ch.combine(Channel.value(params.genus)).combine(Channel.value(params.species)), bakta_db_ch, bakta_proteins_ch)

        // 4. Prepare reference proteome (custom or default H37Rv)
        if (params.ref_annotation && params.ref_genome) {
            log.info "Using custom reference: ${params.ref_annotation} + ${params.ref_genome}"
            PREPARE_REFERENCE(
                Channel.value(file(params.ref_annotation)),
                Channel.value(file(params.ref_genome)),
                bakta_db_ch
            )
            ref_proteome_ch     = PREPARE_REFERENCE.out.ref_proteome
            ref_genome_fasta_ch = PREPARE_REFERENCE.out.ref_genome_fasta
            bad_bakta_ch        = PREPARE_REFERENCE.out.bad_bakta
        } else {
            // Files given on the host are staged as usual. The defaults live inside
            // the container, and BUNDLED_REFERENCE copies them out, together with
            // any of the three given on the host (see there).
            def refs = [params.miniprot_ref_aa, params.miniprot_ref_h37rv, params.miniprot_ref_bad_bakta]
            if (refs.every { file(it).exists() }) {
                ref_proteome_ch     = Channel.value(file(params.miniprot_ref_aa))
                ref_genome_fasta_ch = Channel.value(file(params.miniprot_ref_h37rv))
                bad_bakta_ch        = Channel.value(file(params.miniprot_ref_bad_bakta))
            } else {
                BUNDLED_REFERENCE(
                    stagedIfOnHost(params.miniprot_ref_aa),
                    stagedIfOnHost(params.miniprot_ref_h37rv),
                    stagedIfOnHost(params.miniprot_ref_bad_bakta)
                )
                ref_proteome_ch     = BUNDLED_REFERENCE.out.ref_proteome.first()
                ref_genome_fasta_ch = BUNDLED_REFERENCE.out.ref_genome_fasta.first()
                bad_bakta_ch        = BUNDLED_REFERENCE.out.bad_bakta.first()
            }
        }

        // 5. Align proteins to the genome with Miniprot
        MINIPROT(BAKTA_ANNOTATOR.out.for_miniprot, ref_proteome_ch, ref_genome_fasta_ch, bad_bakta_ch)

        // 6. Classify CDS cases based on alignments
        TYPE_CASES(MINIPROT.out.paf_files.join(BAKTA_ANNOTATOR.out.for_cases).join(names_ch), ref_proteome_ch)

        // 7. Correct CDS based on classified cases
        FIX_CDS(TYPE_CASES.out.for_fix.join(BAKTA_ANNOTATOR.out.for_fix), ref_proteome_ch)

        // What the functional block below hands the interactive report, if it
        // runs: the master tables, the pangenome and the duplication communities.
        report_master_ch = channel.value([])
        report_pangene_ch = channel.value([])
        report_communities_ch = channel.value([])

        // --- FUNCTIONAL ANNOTATION AND DOWNSTREAM ANALYSIS ---
        // This entire block will be skipped if --skip_functional_annotation is used
        if (!params.skip_functional_annotation) {
            // 7. Annotate corrected proteins with eggNOG-mapper (BATCH mode)
            //    Collects all FAA files, runs eggNOG ONCE, splits results back.
            //    This avoids N parallel jobs competing for the same Diamond DB on
            //    shared filesystems: critical for HPC with NFS/Lustre.
            eggnog_faa_ch = FIX_CDS.out.for_eggnog.map { it[1] }.collect()
            eggnog_ids_ch = FIX_CDS.out.for_eggnog.map { it[0] }.collect()

            // Wait for EGGNOG_DB to complete, then pass the folder it built the database in
            eggnog_db_ch = EGGNOG_DB.out.db_path.map { eggnog_dir }
            ANNOT_EGGNOGMAP_BATCH(eggnog_faa_ch, eggnog_ids_ch, eggnog_db_ch)

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
            GET_OUTPUTS(ch_for_downstream.join(names_ch.map { name, _sample, lineage -> tuple(name, lineage) }), mycobrowser_ref_ch)

            // 9. Find and analyze gene duplications
            GET_DUPLICATIONS(ch_for_downstream)
        
            // 10. Runs pangenome for the genes of each genome
            // Each sample's proteins are aligned onto every other sample's
            // genome, which PANGENE tells apart by file name. The genomes
            // PATHOTYPR renamed have names of their own; without typing, each
            // sample's genome from the sheet is copied to <sampleId>.fasta, and
            // that copy is also the one PANGENE leaves out as the sample's own.
            if (!params.skip_typing) {
                ref_fastas_ch = PATHOTYPR.out.typed_fasta.map { it[1] }.collect()
                pangene_ch = ch_for_downstream
            } else {
                PANGENE_GENOME(genomes_ch)
                ref_fastas_ch = PANGENE_GENOME.out.map { _sampleId, genome -> genome }.collect()
                pangene_ch = ch_for_downstream.join(PANGENE_GENOME.out).map { sampleId, amap_tsv, amap_faa, amap_fna, bakta_tsv, _fasta, emapper, genome ->
                    tuple(sampleId, amap_tsv, amap_faa, amap_fna, bakta_tsv, genome, emapper)
                }
            }
            PANGENE(pangene_ch, ref_fastas_ch)

            // 11. Count annotations per category
            GET_COUNT(GET_OUTPUTS.out.master_table.join(names_ch))

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

            // ifEmpty: a step that never runs, all its genomes having failed,
            // must not keep the interactive report from being written.
            report_master_ch = GET_OUTPUTS.out.master_table.map { _name, table -> table }.collect().ifEmpty([])
            report_pangene_ch = REPORT_INPUTS.out.map { inputs -> inputs[6] }.ifEmpty([])
            report_communities_ch = GET_ALLDUP.out.map { tsv, _html -> tsv }.ifEmpty([])
        }

        // 16. The interactive report, results/amap_report.html. It needs only
        // what the correction leaves, so it is written with or without the
        // functional annotation. The sample IDs of the sheet let it name the
        // genomes that produced nothing; it reads the Mycobrowser table as
        // GET_OUTPUTS does (mycobrowser_ref_ch, see stagedIfOnHost).
        sample_ids_ch = genomes_ch.map { sampleId, _fasta -> sampleId }.collectFile(name: 'samples.txt', newLine: true, sort: true)
        AMAP_REPORT(
            FIX_CDS.out.for_downstream.map { row -> row[1] }.collect().ifEmpty([]),
            TYPE_CASES.out.summaries.mix(TYPE_CASES.out.case_tables, TYPE_CASES.out.rvmissing).collect().ifEmpty([]),
            BAKTA_ANNOTATOR.out.bakta_dir.collect().ifEmpty([]),
            report_master_ch,
            ref_proteome_ch,
            report_pangene_ch,
            report_communities_ch,
            SOFTWARE_VERSIONS.out.versions_tsv.mix(DUMP_PARAMS.out, sample_ids_ch).collect(),
            mycobrowser_ref_ch
        )
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
/**
 * Emits the reference genome that ships in the container as a sample, so
 * `-profile test` needs neither a genome in this repository nor one on the
 * host. Runs inside the container, which is the only place /reference exists.
 */
process TEST_GENOME {
    tag 'test_data'

    output:
    tuple val('H37Rv_test'), path('H37Rv_test.fasta')

    script:
    """
    if [ ! -s "${params.miniprot_ref_h37rv}" ]; then
        echo "ERROR: the test profile expects the bundled reference at ${params.miniprot_ref_h37rv}." >&2
        echo "It ships in the AMAP container, so run this with -profile docker,test" >&2
        echo "or -profile singularity,test rather than without a container." >&2
        exit 1
    fi
    cp "${params.miniprot_ref_h37rv}" H37Rv_test.fasta
    """
}

process CREATE_EMPTY_FILE {
    // This process will be cached after its first successful run.
    output:
    path 'empty_summary.log'

    script:
    "touch empty_summary.log"
}
