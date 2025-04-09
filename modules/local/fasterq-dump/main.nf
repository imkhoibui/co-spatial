process FASTERQ_DUMP {
    tag "${asc_id}"
    label 'process_high'

    container 'https://depot.galaxyproject.org/singularity/mulled-v2-5f89fe0cd045cb1d615630b9261a1d17943a9b6a:2f4a4c900edd6801ff0068c2b3048b4459d119eb-0'

    input:
    tuple val(asc_id), path(sra), path(fastq_output)

    output:
    tuple val(asc_id), path('*fastq.gz')                , emit: fastq

    script:
    """
    fasterq-dump --threads $task.cpus \\
        $sra

    pigz --no-name \\
        --processes $task.cpus \\
        *.fastq
    """
}