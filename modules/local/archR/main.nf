process ARCH_R {
    tag "${meta}"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' ?:
            'quay.io/biocontainers/r-archr' }"

    input:

    output:

    script:
    """

    """
}