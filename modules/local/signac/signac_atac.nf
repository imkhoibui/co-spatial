process SIGNAC_ATAC {
    tag "${meta}"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' ?
        'oras://community.wave.seqera.io/library/r-signac_r-seurat:b07da4f886126de7' :
        'community.wave.seqera.io/library/r-signac_r-seurat:29ac9966e457cf5a' }"

    input:
    tuple val(meta), path(rds)
    tuple val(meta), path(fragment_outputs)
    tuple val(meta), path(tissue_dir)
    path spatial_barcodes

    output:
    tuple val(meta), path(peaks)             , emit: peaks

    script:
    """
    #!/usr/bin/env Rscript

    library(Signac)
    library(Seurat)

    # Create Seurat Object
    data <- readRDS("${rds}"}
    DefaultAssay(data) <- "ATAC"
    
    atac_counts <- Read10X_h5("${fragment_outputs}/raw_peak_bc_matrix.h5")
    atac_counts <- atac_counts[, bc_in_tissue]

    NucleosomeSignal(data)
    TSSEnrichment(data)

    
    """
    
}