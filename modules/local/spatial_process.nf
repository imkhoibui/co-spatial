process SPATIAL_PROCESS {
    tag "${meta}"
    label "process_high"

    container "${ workflow.containerEngine == 'singularity' ?
        'oras://community.wave.seqera.io/library/r-signac_r-seurat:b07da4f886126de7' :
        'quay.io/imkhoibui/r-container-sp:updated' }"

    input:
    tuple val(meta), path(rna_outputs)
    tuple val(meta2), path(atac_outputs)
    tuple val(meta3), path(tissue_dir)
    path spatial_barcodes

    output:
    tuple val(meta), path("*spatial.rds")             , emit: rds
    path ("*RNA_UMAP.png")                            , optional: true, emit: rna_umap

    script:
    def seed                = task.ext.seed ?: 1234
    def genome_ref          = task.ext.genome_ref ?: "EnsDb.Mmusculus.v79"
    def args                = task.ext.args ?: ""
    """
    #!/bin/bash
    mkdir callpeaks
    chmod +rw callpeaks
    Rscript ${projectDir}/bin/spatial_preprocess.R --meta_rna ${meta} --meta_atac ${meta2} --meta_genome ${meta3} --rnaseq ${rna_outputs} --atacseq ${atac_outputs} --tissue_dir ${tissue_dir} --spatial_barcodes ${spatial_barcodes}
    """
    
}