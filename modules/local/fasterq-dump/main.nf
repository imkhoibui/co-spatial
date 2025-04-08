process FASTERQ_DUMP {
    tag "${asc_id}"
    label 'process_high'

    container 'https://depot.galaxyproject.org/singularity/sra-tools:3.1.0--h9f5acd7_0'

    input:
    tuple val(asc_id), path(sra), path(fastq_output)

    output:
    tuple val(asc_id), path("*fastq.gz")                , emit: fastq

    script:
    """
    fasterq-dump $sra
    gzip *fastq
    """
}