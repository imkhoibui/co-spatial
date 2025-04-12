process PREFETCH {
    tag "${asc_id}"
    label "process_high"

    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/sra-tools:3.1.0--h9f5acd7_0' :
        'community.wave.seqera.io/library/sra-tools:3.2.1--2063130dadd340c5' }"

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