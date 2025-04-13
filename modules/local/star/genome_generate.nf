process STAR_GENOMEGENERATE {
    tag "${meta}"
    label "process_high_mem"

    container "${ workflow.containerEngine == 'singularity' ?
        'oras://community.wave.seqera.io/library/star:2.7.11b--84fcc19fdfab53a4' :
        'community.wave.seqera.io/library/star:2.7.11b--822039d47adf19a7' }"

    input:
    tuple val(meta), path(fasta)
    path gtf

    output:
    tuple val(meta), path("star")           , emit: index

    script:
    """
    mkdir star
    STAR --runThreadN $task.cpus \\
        --runMode genomeGenerate \\
        --genomeDir star/ \\
        --genomeFastaFiles $fasta \\
        --sjdbGTFfile $gtf
    """
}