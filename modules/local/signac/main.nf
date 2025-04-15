process SIGNAC {
    tag "${meta}"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' ?
        'oras://community.wave.seqera.io/library/r-signac_r-seurat:b07da4f886126de7' :
        'community.wave.seqera.io/library/r-signac_r-seurat:29ac9966e457cf5a' }"

    input:
    tuple val(meta), path(counts)
    tuple val(meta), path(fragments)
    tuple val(meta2), path(tissue_dir)
    path spatial_barcodes

    output:
    tuple val(meta), path(peaks)             , emit: peaks

    script:
    """
    #!/usr/bin/env Rscript

    library(Signac)
    library(Seurat)

    # Create Seurat Object
    ${meta} <- CreateSeuratObject(counts = ${counts})

    # Load spatial barcodes

    spatial_barcodes <- read.csv(${spatial_barcodes})
    
    """
    
}