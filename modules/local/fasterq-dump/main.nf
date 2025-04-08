process FASTERQ_DUMP {
    tag "${asc_id}"
    label 'process_high'

    container 'community.wave.seqera.io/library/sra-tools:3.2.0--7131354b4197d164'

    input:
    tuple val(asc_id), path(sra)

    output:
    tuple val(asc_id), path("*fastq.gz")                , emit: fastq

    script:
    """
    fasterq-dump $sra
    gzip *fastq
    """
}