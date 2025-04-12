process DBIT_PREPROCESS {
    tag "${meta}"
    label "processs_medium"

    container "${ workflow.containerEngine == 'singularity' ?: 
        'community.wave.seqera.io/library/biopython_zip:2b4b64a999cd5274'}"

    input:
    tuple val(meta), path(fastq1), path(fastq2)

    output:
    tuple val(meta), path(fastq1), path("*R2_processed.fastq.gz")      , emit: fastq_processed

    script:
    def raw_fastq1        = "${fastq1}".split("_")
    def processed_fastq2  = raw_fastq1[0] + "_R2_processed.fastq" 
    """
    python3 ${projectDir}/bin/fastq_process.py --input ${fastq2} --output ${processed_fastq2}
    gzip ${processed_fastq2}
    """
}