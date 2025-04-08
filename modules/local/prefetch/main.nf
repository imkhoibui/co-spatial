process PREFETCH {
    tag "${asc_id}"
    label "process_medium"

    container "community.wave.seqera.io/library/sra-tools:3.2.0--7131354b4197d164"

    input:
    val asc_id

    output:
    tuple val(asc_id), path("${asc_id}/.sra")           , emit: sra

    script:
    """
    prefetch $asc_id 
    """
}