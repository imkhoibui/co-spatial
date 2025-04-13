process CUTADAPT {
    tag "${meta}"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' ? 
            'oras://community.wave.seqera.io/library/cutadapt:5.0--57820bb065b39a99' :
            'community.wave.seqera.io/library/cutadapt:5.0--991bbd2e184b7014' }"

    input:

    output:

    script:
    """

    """
}