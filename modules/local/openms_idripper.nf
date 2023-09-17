process OPENMS_IDRIPPER {
    tag "$merge_id.id"
    label 'process_single'

    conda "bioconda::openms=2.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:2.9.1--h135471a_1' :
        'biocontainers/openms:2.9.1--h135471a_1' }"

    input:
        tuple val(merge_id), path(merged_idxml)

    output:
        tuple val(merge_id), path("*.idXML"), emit: ripped
        path "versions.yml"             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args             = task.ext.args  ?: ''

        """
        IDRipper -in $merged_idxml \\
            -out . \\
            -threads $task.cpus \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
