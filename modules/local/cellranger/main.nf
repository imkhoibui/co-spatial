process CELLRANGER {
    tag "${meta}"
    label "process_medium"

    container "nf-core/cellranger:8.0.0"

    input:
    tuple val(meta), path(R1), path(R2), path(R3)
    path ref_atac_genome

    output:
    tuple val(meta), path("outs/*")                 , emit: outputs

    script:
    """
    cellranger-atac count \\
        --id $meta \\
        --reference $ref_atac_genome \\
        --fastqs . \\
        --sample $meta \\
        --localcores $task.cpus \\
        --localmem $task.memory
    """
}