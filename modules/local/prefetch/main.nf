process PREFETCH {
    tag "${asc_id}"
    label "process_high"

    container 'https://depot.galaxyproject.org/singularity/sra-tools:3.1.0--h9f5acd7_0'

    input:
    tuple val(asc_id), val(experiment)

    output:
    tuple val(asc_id), path("${asc_id}/*.sra")           , emit: sra

    script:
    def args            = task.ext.args ?: ""
    """
    prefetch \\
        $asc_id \\
        $args
    """
}