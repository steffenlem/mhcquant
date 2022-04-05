/*
 * Perform the quantification of the samples when the parameter --skip_quantification is not provided
 */

include { OPENMS_MAPALIGNERIDENTIFICATION }                                 from '../../modules/local/openms_mapaligneridentification'
include {
    OPENMS_MAPRTTRANSFORMER as OPENMS_MAPRTTRANSFORMERMZML
    OPENMS_MAPRTTRANSFORMER as OPENMS_MAPRTTRANSFORMERIDXML }               from '../../modules/local/openms_maprttransformer'


workflow PRE_QUANTIFICATION {
    take:
        aligned_hits
        indexed_hits
        mzml_files

    main:
        ch_versions = Channel.empty()
        // Group samples together if they are replicates
        ch_grouped_fdr_filtered = aligned_hits
            .map {
                meta, raw ->
                    [[id:meta.sample + "_" + meta.condition, sample:meta.sample, condition:meta.condition, ext:meta.ext], raw]
                }
            .groupTuple(by: [0])
        // Compute alignment rt transformation
        OPENMS_MAPALIGNERIDENTIFICATION(ch_grouped_fdr_filtered)
        ch_versions = ch_versions.mix(OPENMS_MAPALIGNERIDENTIFICATION.out.versions.first().ifEmpty(null))
        // Intermediate step to join RT transformation files with mzml and idxml channels
        mzml_files
        .join(
            OPENMS_MAPALIGNERIDENTIFICATION.out.trafoxml
                .transpose()
                .flatMap {
                    meta, trafoxml ->
                        ident = trafoxml.baseName.split('_-_')[0]
                        [[[id:ident, sample:meta.sample, condition:meta.condition, ext:meta.ext], trafoxml]]
                }, by: [0] )
        .set { joined_trafos_mzmls }

        indexed_hits
        .join(
            OPENMS_MAPALIGNERIDENTIFICATION.out.trafoxml
                .transpose()
                .flatMap {
                    meta, trafoxml ->
                        ident = trafoxml.baseName.split('_-_')[0]
                        [[[id:ident, sample:meta.sample, condition:meta.condition, ext:meta.ext], trafoxml]]
                }, by: [0] )
        .set { joined_trafos_ids }
        // Align mzML files using trafoXMLs
        OPENMS_MAPRTTRANSFORMERMZML(joined_trafos_mzmls)
        ch_versions = ch_versions.mix(OPENMS_MAPRTTRANSFORMERMZML.out.versions.first().ifEmpty(null))
        // Align unfiltered idXMLfiles using trafoXMLs
        OPENMS_MAPRTTRANSFORMERIDXML(joined_trafos_ids)
        ch_versions = ch_versions.mix(OPENMS_MAPRTTRANSFORMERIDXML.out.versions.first().ifEmpty(null))
        ch_proceeding_idx = OPENMS_MAPRTTRANSFORMERIDXML.out.aligned
            .map {
                meta, raw ->
                [[id:meta.sample + "_" + meta.condition, sample:meta.sample, condition:meta.condition, ext:meta.ext], raw]
            }
            .groupTuple(by: [0])

    emit:
        // Define the information that is returned by this workflow
        versions = ch_versions
        ch_proceeding_idx
        aligned_mzml = OPENMS_MAPRTTRANSFORMERMZML.out.aligned
}