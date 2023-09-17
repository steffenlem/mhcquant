process OPENMS_FEATUREFINDERIDENTIFICATION  {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::openms=2.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:2.9.1--h135471a_1' :
        'biocontainers/openms:2.9.1--h135471a_1' }"

    input:
    tuple val(meta), path(mzml), path(id_int), path(id_ext)

    output:
        tuple val(meta), path("*_out.featureXML"), emit: featurexml
        path "versions.yml"                  , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix           = task.ext.prefix ?: "${meta.sample}_${meta.id}"
        def args             = task.ext.args  ?: ''
        def quant_fdr        = params.quantification_fdr ? "-id $id_int -id_ext $id_ext -svm:min_prob ${params.quantification_min_prob}" : "-id $id_ext"
        args = args + " $quant_fdr"

        """
        FeatureFinderIdentification -in $mzml \\
            -out ${prefix}_out.featureXML \\
            -threads $task.cpus \\
            -candidates_out ${prefix}_candidates.featureXML \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
