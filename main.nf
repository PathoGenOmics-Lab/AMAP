#!/usr/bin/env nextflow
/*
========================================================================================
    A M A P   P I P E L I N E
========================================================================================
    Assisted MTBC Annotator Pipeline
----------------------------------------------------------------------------------------
    Publication:
    DOI:
----------------------------------------------------------------------------------------
    Usage:
    nextflow run PathoGenOmics-Lab/AMAP -profile singularity --tsv <your_samples.tsv>
----------------------------------------------------------------------------------------
*/

// Enable Nextflow DSL 2
nextflow.enable.dsl=2

//
// --- PARAMS ---
//
params.tsv = 'example_files/test.tsv'
params.outdir = 'results'
params.eggnog_db = 'eggnog_db'
params.help = false
params.skip_typing = false
params.skip_functional_annotation = false
params.genus = "Mycobacterium"
params.species = "Mycobacterium tuberculosis"

// Reference files
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
//
// --- MODULE INCLUDES ---
//
include { PATHOTYPR }           from './modules/local/pathotypr'
include { BAKTA_ANNOTATOR }     from './modules/local/bakta'
include { MINIPROT }            from './modules/local/miniprot'
include { TYPE_CASES }          from './modules/local/cases'
include { FIX_CDS }             from './modules/local/cases'
include { EGGNOG_DB }           from './modules/local/eggnog'
include { ANNOT_EGGNOGMAP }     from './modules/local/eggnog'
include { GET_COUNT }           from './modules/local/get'
include { GET_DUPLICATIONS }    from './modules/local/get'
include { GET_OUTPUTS }         from './modules/local/get'
include { GET_ALLDUP }          from './modules/local/get'
include { PANGENE }             from './modules/local/pangene'
include { REPORT_INPUTS }       from './modules/local/freport'
include { FINAL_REPORT }        from './modules/local/freport'
//
// --- INPUT CHANNEL ---
//
Channel
    .fromPath(params.tsv)
    .ifEmpty { error "Input TSV file not found: ${params.tsv}" }
    .splitCsv(header: true, sep: '\t')
    .map { row -> tuple(row.sampleId, file(row.fasta)) }
    .set { fasta_ch }

//
// --- WORKFLOW ---
//
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
        
        FIX_CDS.out.for_eggnog
            .map { sampleId, faa, fna -> tuple(sampleId, faa) }
            .set { for_eggnog_ch }

        ANNOT_EGGNOGMAP(for_eggnog_ch, EGGNOG_DB.out.db_path)
        
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
