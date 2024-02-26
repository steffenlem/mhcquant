/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowMhcquant.initialise(params, log)

// Input/output options
if (params.input)   { sample_sheet = file(params.input) }
if (params.fasta)   { params.fasta = params.fasta       }

// MHC affinity prediction
if (params.predict_class_1 || params.predict_class_2) {
    Channel.from(file(params.allele_sheet, checkIfExists: true))
        .splitCsv(header: true, sep:'\t')
        .multiMap { col ->
            classI: ["${col.Sample}", "${col.HLA_Alleles_Class_1}"]
            classII: ["${col.Sample}", "${col.HLA_Alleles_Class_2}"] }
        .set { ch_alleles_from_sheet }

        // Allele class 1
        if (params.predict_class_1) {
            ch_alleles_from_sheet.classI
                .ifEmpty { exit 1, "params.allele_sheet was empty - no allele input file supplied" }
                .flatMap { it -> [tuple(it[0].toString(), it[1])] }
                .set { peptides_class_1_alleles }
        }

        // Allele class 2
        if (params.predict_class_2) {
            ch_alleles_from_sheet.classII
                .ifEmpty { exit 1, "params.allele_sheet was empty - no allele input file supplied" }
                .flatMap { it -> [tuple(it[0].toString(), it[1])] }
                .set { peptides_class_2_alleles }
        }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { OPENMS_DECOYDATABASE }                                            from '../modules/local/openms_decoydatabase'
include { THERMORAWFILEPARSER }                                             from '../modules/local/thermorawfileparser'
include { TDF2MZML }                                                        from '../modules/local/tdf2mzml'
include { OPENMS_PEAKPICKERHIRES }                                          from '../modules/local/openms_peakpickerhires'
include { OPENMS_FILEFILTER }                                               from '../modules/local/openms_filefilter'
include { OPENMS_COMETADAPTER }                                             from '../modules/local/openms_cometadapter'
include { OPENMS_PEPTIDEINDEXER }                                           from '../modules/local/openms_peptideindexer'
include { MS2RESCORE }                                                      from '../modules/local/ms2rescore'
include { OPENMS_IDSCORESWITCHER }                                          from '../modules/local/openms_idscoreswitcher'

include { OPENMS_IDFILTER as OPENMS_IDFILTER_Q_VALUE }                      from '../modules/local/openms_idfilter'
include { OPENMS_IDMERGER }                                                 from '../modules/local/openms_idmerger'

include { OPENMS_PSMFEATUREEXTRACTOR }                                      from '../modules/local/openms_psmfeatureextractor'
include { OPENMS_PERCOLATORADAPTER }                                        from '../modules/local/openms_percolatoradapter'
include { PYOPENMS_IONANNOTATOR }                                           from '../modules/local/pyopenms_ionannotator'

include { OPENMS_TEXTEXPORTER }                                             from '../modules/local/openms_textexporter'
include { OPENMS_MZTABEXPORTER }                                            from '../modules/local/openms_mztabexporter'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []
// Sort closure for merging and splitting files
def sortById = { a, b -> a.id <=> b.id }

include { INCLUDE_PROTEINS }                                                from '../subworkflows/local/include_proteins'
include { REFINE_FDR }                                                      from '../subworkflows/local/refine_fdr'
include { QUANT }                                                           from '../subworkflows/local/quant'
include { PREDICT_CLASS1 }                                                  from '../subworkflows/local/predict_class1'
include { PREDICT_CLASS2 }                                                  from '../subworkflows/local/predict_class2'

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
workflow MHCQUANT {

    ch_versions = Channel.empty()

    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    INPUT_CHECK.out.ms_runs
        .branch {
            meta, filename ->
                raw : meta.ext == 'raw'
                    return [ meta.subMap('id', 'sample', 'condition'), filename ]
                mzml : meta.ext == 'mzml'
                    return [ meta.subMap('id', 'sample', 'condition'), filename ]
                tdf : meta.ext == 'd'
                    return [ meta.subMap('id', 'sample', 'condition'), filename ]
                idxml : meta.ext == 'idxml'
                    return [ meta.subMap('id', 'sample', 'condition'), filename ]
                other : true }
        .set { branched_ms_files }
    
    ch_ms_files = branched_ms_files.idxml.map{ meta, idxml -> [[id: meta.sample + '_' + meta.condition], idxml[0]]}.groupTuple()

    if (params.rescoring_engine == 'percolator') {
        // Run Percolator
        OPENMS_PERCOLATORADAPTER(ch_ms_files)
        ch_versions = ch_versions.mix(OPENMS_PERCOLATORADAPTER.out.versions.ifEmpty(null))
        ch_rescored_runs = OPENMS_PERCOLATORADAPTER.out.idxml
    } else {
        log.warn "The rescoring engine is set to mokapot. This rescoring engine currently only supports psm-level-fdr via ms2rescore."
        // TODO: remove whitelist argument from idscoreswitcher
        OPENMS_IDSCORESWITCHER(OPENMS_IDMERGER.out.idxml)
        ch_rescored_runs = OPENMS_IDSCORESWITCHER.out.switched_idxml.map { tuple -> tuple.findAll { it != [] }}
    }

    // Filter by percolator q-value
    // TODO: Use empty list instead of null
    OPENMS_IDFILTER_Q_VALUE(ch_rescored_runs.flatMap { it -> [tuple(it[0], it[1], null)] })
    ch_versions = ch_versions.mix(OPENMS_IDFILTER_Q_VALUE.out.versions.ifEmpty(null))

    //
    // SUBWORKFLOW: Refine the FDR values on the predicted subset
    //
    if (params.refine_fdr_on_predicted_subset && params.predict_class_1) {
        // Run the following subworkflow
        REFINE_FDR (
            OPENMS_IDFILTER_Q_VALUE.out.idxml,
            OPENMS_PSMFEATUREEXTRACTOR.out.idxml,
            peptides_class_1_alleles
        )
        ch_versions = ch_versions.mix(REFINE_FDR.out.versions.ifEmpty(null))
        // Define the outcome of the paramer to a fixed variable
        filter_q_value = REFINE_FDR.out.filter_refined_q_value
    } else {
        // Make sure that the columns that consists of the ID's, sample names and the idXML file names are returned
        filter_q_value = OPENMS_IDFILTER_Q_VALUE.out.idxml
    }

    //
    // SUBWORKFLOW: QUANT
    //
    if (!params.skip_quantification) {
        QUANT(merge_meta_map, ch_rescored_runs, filter_q_value, ch_clean_mzml_file)
        ch_versions = ch_versions.mix(QUANT.out.versions.ifEmpty(null))
        ch_output = QUANT.out.consensusxml
    } else {
        ch_output = filter_q_value
    }

    // Prepare for check if file is empty
    OPENMS_TEXTEXPORTER(ch_output)
    ch_versions = ch_versions.mix(OPENMS_TEXTEXPORTER.out.versions.ifEmpty(null))
    // Return an error message when there is only a header present in the document
    OPENMS_TEXTEXPORTER.out.tsv.map {
        meta, tsv -> if (tsv.size() < 130) {
        log.warn "It seems that there were no significant hits found for this sample: " + meta.sample + "\nPlease consider incrementing the '--fdr_threshold' after removing the work directory or to exclude this sample. "
        }
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
