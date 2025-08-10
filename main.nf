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
    nextflow run ${workflow.manifest.name} -profile singularity [options]

    Required arguments:
    --tsv                       Path to the input TSV file.
                                Format: sampleId<TAB>fasta_path

    Optional arguments:
    --outdir                    Path to the output directory. (Default: ${params.outdir})
    --eggnog_db                 Path to the eggNOG database directory. (Default: ${params.eggnog_db})
    --skip_typing               Skip the PATHOTYPR lineage typing step. (Default: ${params.skip_typing})
    --help                      Display this help message.

    Reference file arguments (advanced):
    --pathotyper_snps     Path to the SNPs TSV for MTB_TYPER.
    --pathotyper_ref      Path to the reference FASTA for MTB_TYPER.
    --miniprot_ref_aa     Path to the H37Rv amino acid reference for MINIPROT.
    ...and other reference files.
    """
    exit 0
}

log.info logo
log.info """
    fastas          : ${params.tsv}
    output          : ${params.outdir}
    eggnog_db       : ${params.eggnog_db}
    skip_typing     : ${params.skip_typing}
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

    // FIX 1: Transform the output of FIX_CDS to match what ANNOT_EGGNOGMAP expects.
    // The .map operator takes the 3-element tuple (id, faa, fna) and returns
    // the 2-element tuple (id, faa) that ANNOT_EGGNOGMAP needs.
    FIX_CDS.out.for_eggnog
        .map { sampleId, faa, fna -> tuple(sampleId, faa) }
        .set { for_eggnog_ch }

    // 7. Annotate corrected proteins with eggNOG-mapper
    ANNOT_EGGNOGMAP(for_eggnog_ch, EGGNOG_DB.out.db_path)

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

    // --- FIX 2: CHANNEL HANDLING FOR FINAL REPORT ---
    // 12. Generate final reports
    
    // First, collect the outputs from the upstream processes.
    // .collect() will gather all files emitted by each process into a single list.
    summaries_cases_list_ch = TYPE_CASES.out.summaries.collect()
    summaries_annot_list_ch = GET_COUNT.out.collect()

    // 13. Pass these two lists directly to REPORT_INPUTS.
    // This process is designed to take two inputs: the list of annotation files and the list of case files.
    REPORT_INPUTS(summaries_annot_list_ch, summaries_cases_list_ch)

    // 14. Analyze duplications across all samples
    GET_ALLDUP(GET_DUPLICATIONS.out.dup_fastas.collect())

    // 15. Render the final HTML report
    FINAL_REPORT(REPORT_INPUTS.out)

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
