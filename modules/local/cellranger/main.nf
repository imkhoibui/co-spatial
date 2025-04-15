process CELLRANGER {
    tag "${meta}"
    label "process_high"

    container "quay.io/nf-core/cellranger-atac:2.1.0"

    input:
    tuple val(meta), path(R1), path(R2), path(R3)
    path ref_atac_genome

    output:
    tuple val(meta), path("outs/*")                 , emit: outputs

    script:
    def config_in               = task.ext.config_in ?: ""
    def version                 = task.ext.varsion ?: ""
    def cratac_txt              = task.ext.cratac_txt ?: ""
    def use_custom_cratac_file  = task.ext.use_cratac ?: ""
    def args                    = task.ext.args ?: ""
    """
    #!/bin/bash
    if [ $use_custom_cratac_file ] && cp $cratac_txt /opt/cellranger-atac-2.1.0/lib/python/atac/barcodes/
    mkdir fastq_files
    cp *fastq.gz fastq_files/

    cellranger-atac count \\
        --id $meta \\
        --reference $ref_atac_genome \\
        --fastqs fastq_files/ \\
        --sample $meta \\
        --localcores $task.cpus \\
        $args
    """
}