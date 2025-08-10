#!/usr/bin/env nextflow
/*
========================================================================================
    A M A P   P I P E L I N E
========================================================================================
    Assisted MTBC Annotator Pipeline
----------------------------------------------------------------------------------------
    Usage:
    nextflow run PathoGenOmics-Lab/AMAP -profile <docker/singularity/slurm> --tsv <your_samples.tsv>
----------------------------------------------------------------------------------------
*/

// Enable Nextflow DSL 2
nextflow.enable.dsl=2

// --- LOGO and HELP MESSAGE ---
def logo = """
=======================================
        A M A P   P I P E L I N E
=======================================
           _,--',   _._.--._____
   .--.--';_'-.', ";_         _.,-'
   .'--'.  _.'     {`'-;_ .-.>.'
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

    Required arguments:
    --tsv                       Path to the input TSV file.
                                Format: sampleId<TAB>fasta_path

    Optional arguments:
    --outdir                    Path to the output directory. (Default: results)
    --skip_typing               Skip the PATHOTYPR lineage typing step. (Default: false)
    --help                      Display this help message.

    All other parameters, such as reference file paths, can be configured in the
    `nextflow.config` file or overridden on the command line.
    """
    exit 0
}

log.info logo
log.info """
    Input TSV       : ${params.tsv}
    Output dir      : ${params.outdir}
    Skip Typing     : ${params.skip_typing}
"""

// --- MODULE INCLUDES ---
include { PATHOTYPR }           from './modules/pathotypr'
include { BAKTA_ANNOTATOR }     from './modules/bakta'
include { MINIPROT }            from './modules/miniprot'
include { TYPE_CASES }          from './modules/cases'
include { FIX_CDS }             from './modules/cases'
include { EGGNOG_DB }           from './modules/eggnog'
include { ANNOT_EGGNOGMAP }     from './modules/eggnog'
include { GET_COUNT }           from './modules/get'
include { GET_DUPLICATIONS }    from './modules/get'
include { GET_OUTPUTS }         from './modules/get'
include { GET_ALLDUP }          from './modules/get'
include { PANGENE }             from './modules/pangene'
include { REPORT_INPUTS }       from './modules/freport'
include { FINAL_REPORT }        from './modules/freport'

// --- INPUT CHANNEL ---
Channel
    .fromPath(params.tsv)
    .ifEmpty { error "Input TSV file not found: ${params.tsv}" }
    .splitCsv(header: true, sep: '\t')
    .map { row -> tuple(row.sampleId, file(row.fasta)) }
    .set { fasta_ch }

// --- WORKFLOW ---
workflow {
    // 1. Setup eggNOG database if not present
    EGGNOG_DB()

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

    // 4. Align proteins to the genome with Miniprot
    MINIPROT(BAKTA_ANNOTATOR.out.for_miniprot)

    // 5. Classify CDS cases based on alignments
    TYPE_CASES(MINIPROT.out.paf_files.join(BAKTA_ANNOTATOR.out.for_cases))

    // 6. Correct CDS based on classified cases
    FIX_CDS(TYPE_CASES.out.for_fix.join(BAKTA_ANNOTATOR.out.for_fix))

    // 7. Annotate corrected proteins with eggNOG-mapper
    ANNOT_EGGNOGMAP(FIX_CDS.out.for_eggnog, EGGNOG_DB.out.db_path)

    // Combine results for downstream processing
    ch_for_downstream = FIX_CDS.out.for_downstream.join(ANNOT_EGGNOGMAP.out)

    // 8. Generate final output files (GFF, GBK, Master Table)
    GET_OUTPUTS(ch_for_downstream)

    // 9. Find and analyze gene duplications
    GET_DUPLICATIONS(ch_for_downstream)
    
    // 10. Runs pangenome for the genes of each genome
    PANGENE(ch_for_downstream)

    // 11. Count annotations per category
    GET_COUNT(GET_OUTPUTS.out.master_table)

    // 12. Generate final reports
    // Collect summaries from different processes
    summaries_cases_ch = TYPE_CASES.out.summaries.collect()
    summaries_annot_ch = GET_COUNT.out.collect()

    // 13. Create report inputs
    REPORT_INPUTS(summaries_annot_ch, summaries_cases_ch)

    // 14. Analyze duplications across all samples
    GET_ALLDUP(GET_DUPLICATIONS.out.dup_fastas.collect())

    // 15. Render the final HTML report
    FINAL_REPORT(REPORT_INPUTS.out)
}

// --- HELPER PROCESSES ---
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
