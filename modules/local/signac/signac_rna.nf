process SIGNAC_RNA {
    tag "${meta}"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' ?
        'oras://community.wave.seqera.io/library/r-signac_r-seurat:b07da4f886126de7' :
        'quay.io/imkhoibui/r-container-sp' }"

    input:
    tuple val(meta), path(rna_outputs)
    tuple val(meta), path(atac_outputs)
    tuple val(meta2), path(tissue_dir)
    path spatial_barcodes

    output:
    tuple val(meta), path(rds)             , emit: rds

    script:
    def seed                = task.ext.seed ?: 1234
    def genome_ref          = task.ext.genome_ref ?: "EnsDb.Mmusculus.v79"
    def args                = task.ext.args ?: ""
    """
    #!/usr/bin/env Rscript
    Rscript "${projectDir}/bin/spatial_preprocess.R" \\
        --meta ${meta} \\
        --rnaseq ${rna_outputs} \\
        --atacseq ${atac_outputs}
        --tissue_dir ${tissue_dir} \\
        --spatial_barcodes ${spatial_barcodes}
    """
    
}