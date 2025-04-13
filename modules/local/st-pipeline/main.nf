process ST_PIPELINE {
    tag "${meta}"
    label "process_high_mem"

    container "${ workflow.containerEngine == 'singularity' ?: 'quay.io/imkhoibui/stpipeline:updated' }"

    input:
    tuple val(meta), path(fastq1), path(fastq2)
    path spatial_barcodes
    tuple val(meta2), path(ref_map)
    path ref_annotation
    
    output:
    tuple val(meta), path("st_pipeline")            , emit: st_pipeline

    script:
    def args            = task.ext.args ?: ""
    def outdir          = "st_pipeline/"
    """
    #!/bin/bash
    mkdir $outdir
    st_pipeline_run \\
        --output-folder $outdir \\
        --ids $spatial_barcodes \\
        --ref-map $ref_map \\
        --ref-annotation $ref_annotation \\
        --expName $meta \\
        --threads $task.cpus \\
        $args \\
        $fastq2 $fastq1
    """
}