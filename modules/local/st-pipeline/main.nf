process ST_PIPELINE {
    tag "${meta}"
    label "process_medium"

    container ""

    input:

    output:

    script:
    """
    #!/usr/env/bin python3
    st_pipeline_run.py --ids \\
        --ref-map \\
        --ref-annotation \\
        --expName \\
        $fastq1 $fastq2
    """
}