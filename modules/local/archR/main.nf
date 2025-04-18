process ARCH_R {
    tag "${meta}"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' ?
       'oras://community.wave.seqera.io/library/r-signac_r-seurat:b07da4f886126de7' :
        'quay.io/imkhoibui/r-container-sp:updated' }"

    input:
    tuple val(meta), path(rds)
    path spatial_barcodes

    output:
    tuple val(meta), path(peaks)             , emit: peaks

    script:
    """
    #!/usr/bin/env Rscript
    library(paletteer)#
    library(Signac)
    library(Seurat)
    library(GenomeInfoDb)
    library(EnsDb.Mmusculus.v79)
    library(ggplot2)
    library(patchwork)
    library(SeuratObject)
    library(dplyr)
    library(ggpubr)
    library(Matrix)
    set.seed(1234)

    col <- paletteer_d("ggthemes::Classic_20")[c(20,6,9,19,3,18,8,7,5,1,13,2,17,15,10,11,12,16)]

    sampleNames <- '${meta}'
    data <- readRDS("${rds}")

    DefaultAssay(data) <- 'ATAC'
    data <- FindTopFeatures(data, min.cutoff = 10)
    data <- RunTFIDF(data)
    data <- RunSVD(data)
    data <- FindNeighbors(
        object = data,
        reduction = 'lsi',
        dims = 2:20
    )
    data <- FindClusters(
        object = data,
        algorithm = 3,
        resolution = 0.6,
        verbose = FALSE
    )
    data <- RunUMAP(data, reduction = 'lsi', dims = 2:20, reduction.name = 'umap.cut')

    p2 <- DimPlot(data, reduction = 'umap.cut', label = TRUE) + NoLegend() + ggtitle("ATAC UMAP")

    write.csv(data@reductions\$umap.cut@cell.embeddings, "/data_all_ATAC_UMAP_coordinates.csv", row.names = T, col.names = T, quote=F)

    DefaultAssay(data) <- "peaks"
    data <- FindTopFeatures(data, min.cutoff = 5)
    data <- RunTFIDF(data)
    data <- RunSVD(data)

    # build a joint neighbor graph using both assays
    data <- FindMultiModalNeighbors(
        object = data,
        reduction.list = list("pca", "lsi"), 
        dims.list = list(1:10, 2:20),
        modality.weight.name = "RNA.weight",
        verbose = TRUE
    )

    data <- RunUMAP(
        object = data,
        nn.name = "weighted.nn",
        assay = "RNA",
        verbose = TRUE,
        reduction.name = "wnn.umap", 
        reduction.key = "wnnUMAP_"
    )

    data <- FindClusters(data, graph.name = "wsnn", algorithm = 3, resolution = 0.8, verbose = FALSE)
    meta.data <- data@meta.data
    new_row_names <- row.names(meta.data)
    new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
    row.names(meta.data) <- new_row_names

    # calculating gene score matrix
    addGeneScoreMatrix(input=,
    
    )

    spatial_archr_data <- read.csv(file = "${meta}_gene_score_matrix.csv")
    row.names(spatial_archr_data) <- spatial_archr_data\$X
    spatial_archr_data <- spatial_archr_data[,-1]

    assay = "Spatial"
    filter.matrix = TRUE
    slice = "slice1"

    object <- CreateSeuratObject(counts = spatial_archr_data, assay = assay, meta.data = meta.data)
    image <- Read10X_Image(image.dir = file.path(${tissue_dir}, "spatial"), filter.matrix = filter.matrix)
    image <- image[Cells(x = object)]
    DefaultAssay(object = image) <- assay
    object[[slice]] <- image

    spatial.obj <- object
    table(data\$wsnn_res.0.8)
    table(data\$SCT_snn_res.0.8)
    table(data\$ATAC_snn_res.0.6)

    p6 <- SpatialPlot(spatial.obj, label = FALSE, label.size = 1.8, group.by = 'wsnn_res.0.8', pt.size.factor = 1, cols = col, image.alpha = 1, stroke = 0)
    p6\$layers[[1]]\$aes_params <- c(p6\$layers[[1]]\$aes_params, shape=22)

    ggsave("Spatial-wsnn_res.0.8-UMP.png", plot = p6, width = 9, height = 9)
    p7 <- SpatialPlot(spatial.obj, label = FALSE, label.size = 1.8, group.by = 'SCT_snn_res.0.8', pt.size.factor = 1, cols = col, image.alpha = 1, stroke = 0)
    p7\$layers[[1]]\$aes_params <- c(p7\$layers[[1]]\$aes_params, shape=22)

    ggsave("Spatial-RNA_res.0.8-UMP.png", plot = p7, width = 9, height = 9)

    p8 <- SpatialPlot(spatial.obj, label = FALSE, label.size = 1.8, group.by = 'ATAC_snn_res.0.6', pt.size.factor = 1, cols = col, image.alpha = 1, stroke = 0)
    p8\$layers[[1]]\$aes_params <- c(p8\$layers[[1]]\$aes_params, shape=22)

    ggsave("ATAC_spatial-UMP.png", plot = p8, width = 9, height = 9)
    p2 <- p7+p8+p6

    ggsave("wnn-3-spatial-UMP.png", plot = p2, width = 15, height = 9)

    p1 <- DimPlot(data, label = TRUE, cols = col, group.by = "SCT_snn_res.0.8", pt.size = 1.8) + NoLegend() + ggtitle("RNA UMAP")
    p2 <- DimPlot(data, label = TRUE, cols = col, group.by = "ATAC_snn_res.0.6", pt.size = 1.8) + NoLegend() + ggtitle("ATAC UMAP")
    p3 <- DimPlot(data, label = TRUE, cols = col, group.by = "wsnn_res.0.8", pt.size = 1.8) + NoLegend() + ggtitle("WNN UMAP")
    p <- p1 + p2 + p3
    ggsave("3-UMP-wnn.png", plot = p, width = 15, height = 9)

    write.csv(data@reductions\$umap.rna@cell.embeddings, "${meta}_RNA_UMAP_coordinates-3.csv", row.names = T, col.names = T, quote=F)
    write.csv(data@reductions\$umap.cut@cell.embeddings, "${meta}_ATAC_UMAP_coordinates-3.csv", row.names = T, col.names = T, quote=F)
    write.csv(data@reductions\$wnn.umap@cell.embeddings, "${meta}_wnn_UMAP_coordinates-3.csv", row.names = T, col.names = T, quote=F)

    """
    
}