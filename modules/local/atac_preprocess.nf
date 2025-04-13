process ATAC_PREPROCESS {
    tag "${meta}"
    label "process_low"

    container "${ workflow.containerEngine == 'singularity' ?: 
        'community.wave.seqera.io/library/biopython_zip:2b4b64a999cd5274'}"

    input: 
    tuple val(meta), path(fastq1), path(fastq2)

    output:
    tuple val(meta), path("*_R1*.fastq.gz"), path("*_R2*.fastq.gz"), path("*_R3*.fastq.gz")       , emit: fastq

    script:
    def output_R1 = "${meta}_S1_L001_R1_001.fastq"
    def output_R2 = "${meta}_S1_L001_R2_001.fastq"
    def output_R3 = "${meta}_S1_L001_R3_001.fastq.gz"
    """
    #!/bin/bash 
    python3 ${projectDir}/bin/BC_process.py --input ${fastq2} 
        --output_R1 $output_R1 \\
        --output_R2 $output_R2
    
    gzip $output_R1
    gzip $output_R2
    mv $fastq3 $output_R3
    """
}