process SIGNAC_RNA {
    tag "${meta}"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' ?
        'oras://community.wave.seqera.io/library/r-signac_r-seurat:b07da4f886126de7' :
        'community.wave.seqera.io/library/r-signac_r-seurat:29ac9966e457cf5a' }"

    input:
    tuple val(meta), path(counts)
    tuple val(meta), path(tissue_dir)
    path spatial_barcodes

    output:
    tuple val(meta), path(rds)             , emit: rds

    script:
    def genome_ref          = task.ext.genome_ref ?: ""
    """
    #!/usr/bin/env Rscript

    library(Signac)
    library(Seurat)
    library(${genome_ref})

    # Creating Seurat Object
    data <- CreateSeuratObject(counts = ${counts}, header = TRUE, row.names = 1, as.is = TRUE)
    data <- as.data.frame(t(data))

    data\$V1 <- paste0(row.names(data))
    data\$V1 <- paste0(data\$V1, "-1")
    row.names(data) <- data\$V1
    data <- data[,!(colnames(data) %in% c("V1"))]

    RNA_counts <- as.data.frame(t(data))
    combined_object <- CreateSeuratObject(counts = RNA_counts, assay = "RNA")

    # Loading spatial barcodes
    spatial_barcodes <- read.table("${spatial_barcodes}", sep =",", header = FALSE, dec =".", stringsAsFactors = F)
    x <- as.character(spatial_barcodes[1,])
    x = x[-1]
    spatial_barcodes = spatial_barcodes[-1]
    spatial_barcodes = as.data.frame(t(spatial_barcodes))

    positions <- read.table(file = '${tissue_dir}/tissue_positions_list.csv', 
                        sep = '\t', stringsAsFactors=FALSE)
    positions\$V4 <- paste0(positions\$V2,"x",positions\$V3)
    row.names(positions) <- positions\$V4
    positions_in_tissue <- positions[spatial_barcodes\$V1, ]
    positions_in_tissue\$V1 <- paste0(positions_in_tissue\$V1, "-1")
    bc_in_tissue <- positions_in_tissue\$V1 
    
    """
    
}