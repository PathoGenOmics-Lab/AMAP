#!/usr/bin/env nextflow
/*
========================================================================================
    A M A P   P I P E L I N E
========================================================================================
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

    Optional arguments:
    --outdir                    Path to the output directory. (Default: results)
    --skip_typing               Skip the PATHOTYPR lineage typing step. (Default: false)
    --skip_functional_annotation Skip all functional annotation steps. (Default: false)
    --genus                     Genus for Bakta annotation. (Default: "Mycobacterium")
    --species                   Species for Bakta annotation. (Default: "Mycobacterium tuberculosis")
    --help                      Display this help message.

    Reference File Arguments (override paths in nextflow.config):
    --eggnog_db                 Path to the eggNOG database. (Default: 'eggnog_db')
    --pathotyper_snps           Path to the SNPs TSV for Pathotyper.
    --pathotyper_ref            Path to the reference FASTA for Pathotyper.
    --miniprot_ref_aa           Path to the H37Rv amino acid reference for Miniprot.
    --miniprot_ref_bad_bakta    Path to the 'bad bakta' reference for Miniprot.
    --miniprot_ref_h37rv        Path to the H37Rv reference genome for Miniprot.
    --bakta_proteins            Path to the proteins GBK file for Bakta.
    --mycobrowser_ref           Path to the Mycobrowser reference TXT file.
    """
    exit 0
}

log.info logo
log.info """
    Input TSV                  : ${params.tsv}
    Output dir                 : ${params.outdir}
    Skip Typing                : ${params.skip_typing}
    Skip Functional Annotation : ${params.skip_functional_annotation}
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
    // These steps always run
    if (!params.skip_typing) {
        PATHOTYPR(fasta_ch)
        typed_ch = PATHOTYPR.out.typed_fasta
    } else {
        CREATE_EMPTY_FILE()
        typed_ch = fasta_ch.combine(CREATE_EMPTY_FILE.out).map { sampleId, fasta, empty_file ->
            tuple(sampleId, fasta, empty_file)
        }
    }
    BAKTA_ANNOTATOR(typed_ch.combine(Channel.value(params.genus)).combine(Channel.value(params.species)))
    MINIPROT(BAKTA_ANNOTATOR.out.for_miniprot)
    TYPE_CASES(MINIPROT.out.paf_files.join(BAKTA_ANNOTATOR.out.for_cases))
    FIX_CDS(TYPE_CASES.out.for_fix.join(BAKTA_ANNOTATOR.out.for_fix))


    // --- CONDITIONAL EXECUTION BLOCK ---
    if (!params.skip_functional_annotation) {
        // Run the full pipeline with functional annotation
        log.info "Running full pipeline with functional annotation."

        EGGNOG_DB()
        ANNOT_EGGNOGMAP(FIX_CDS.out.for_eggnog, EGGNOG_DB.out.db_path)
        ch_for_downstream = FIX_CDS.out.for_downstream.join(ANNOT_EGGNOGMAP.out)
        
        GET_OUTPUTS(ch_for_downstream)
        GET_DUPLICATIONS(ch_for_downstream)
        PANGENE(ch_for_downstream)
        GET_COUNT(GET_OUTPUTS.out.master_table)
        
        summaries_cases_ch = TYPE_CASES.out.summaries.collect()
        summaries_annot_ch = GET_COUNT.out.collect()
        
        REPORT_INPUTS(summaries_annot_ch, summaries_cases_ch)
        GET_ALLDUP(GET_DUPLICATIONS.out.dup_fastas.collect())
        FINAL_REPORT(REPORT_INPUTS.out)

    } else {
        // Skip functional annotation and all downstream steps
        log.info "Skipping functional annotation. The pipeline will stop after CDS correction."
        // The results from TYPE_CASES and FIX_CDS will be published automatically
        // as defined in their respective modules. No further steps are needed.
    }
}

// --- HELPER PROCESSES ---
process CREATE_EMPTY_FILE {
    output:
    path 'empty_summary.log'
    script:
    "touch empty_summary.log"
}
