process CUTADAPT {
    tag "${meta}"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' ? 
            'oras://community.wave.seqera.io/library/cutadapt:5.0--57820bb065b39a99' :
            'community.wave.seqera.io/library/cutadapt:5.0--991bbd2e184b7014' }"

    input:
    tuple val(meta), path(fastq1), path(fastq2)

    output:
    tuple val(meta), path("*_R1_fastq.gz"), path("*_R2.fastq.gz")       , emit: filter
    path "*cutadapt.log"                                                , emit: log

    script:
    def args        = task.ext.args ?: ""
    """
    cutadapt $args \\
        --discard-untrimmed \\
        -e 0.1 \\
        -j $task.cpus \\
        -o ${meta}_${task.name}_R1.fastq.gz \\
        -p ${meta}_${task.name}_R2.fastq.gz \\
        $fastq1 $fastq2 > ${meta}_${task.name}_cutadapt.log
    """
}