# install.packages("optparse", repos = "https://cloud.r-project.org")
remove.packages("Matrix")
remove.packages("irlba")
remotes::install_github("cran/Matrix@1.6-5", force = TRUE)
install.packages("irlba", type = "source", force = TRUE, repos = "https://cloud.r-project.org")
suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratObject))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(EnsDb.Mmusculus.v79))
suppressMessages(library(patchwork))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(ggpubr))
suppressMessages(library(Matrix))
suppressMessages(library(paletteer))
suppressMessages(library(optparse))
suppressMessages(library(biovizBase))

## Extracting arguments
set.seed(1234)
option_list <- list(
    make_option("--meta_rna", type = "character"),
    make_option("--meta_atac", type = "character"),
    make_option("--meta_genome", type = "character"),
    make_option("--rnaseq", type = "character"),
    make_option("--atacseq", type = "character"),
    make_option("--tissue_dir", type = "character"),
    make_option("--spatial_barcodes", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
meta_rna         <- opt$meta_rna
meta_atac        <- opt$meta_atac
meta_genome      <- opt$meta_genome
rna_counts       <- opt$rnaseq
atac_fragments   <- opt$atacseq
tissue_dir       <- opt$tissue_dir
spatial_barcodes <- opt$spatial_barcodes

## Loading files
rna_data <- read.table(file = paste0("st_pipeline/", meta_rna, "_stdata.tsv"),
                   header = TRUE, row.names = 1, as.is = TRUE)

spatial_barcodes <- read.table(file="spatial_barcodes.txt", header=FALSE, row.names=1, as.is=TRUE)
spatial_barcodes$V4 <- paste0(spatial_barcodes$V2,"x",spatial_barcodes$V3)
matching <- match(rownames(rna_data), spatial_barcodes$V4)
rownames(rna_data)[!is.na(matching)] <- rownames(spatial_barcodes)[matching[!is.na(matching)]]
rna_data$V1 <- paste0(row.names(rna_data))
rna_data$V1 <- paste0(rna_data$V1, "-1")
row.names(rna_data) <- rna_data$V1
rna_data <- rna_data[,!(colnames(rna_data) %in% c("V1"))]

rna_counts <- as.data.frame(t(rna_data))
obj <- CreateSeuratObject(counts = rna_counts, assay = "RNA")
meta.data <- obj@meta.data
positions <- read.table("spatial/tissue_positions_list.csv", sep =",", header = FALSE, row.names = 1, dec =".", stringsAsFactors = F)

rownames(spatial_barcodes) <- as.character(rownames(spatial_barcodes))
rownames(positions) <- as.character(rownames(positions))
my_data_in_tissue <- spatial_barcodes[rownames(positions), ]
my_data_in_tissue$V1 <- paste0(rownames(my_data_in_tissue), "-1")
bc_in_tissue <- my_data_in_tissue$V1 

## Gene annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
atac_data <- Read10X_h5("outs/raw_peak_bc_matrix.h5")
atac_data <- atac_data[, bc_in_tissue]
obj[["ATAC"]] <- CreateChromatinAssay(
    counts = atac_data,
    sep = c(":", "-"),
    genome = "mm10",
    fragments = "outs/fragments.tsv.gz",
    min.cells = 0
)
Annotation(obj[["ATAC"]]) <- annotations

## ATAC analysis
DefaultAssay(obj) <- "ATAC"
obj <- NucleosomeSignal(obj)
obj <- TSSEnrichment(obj)
obj$blacklist_fraction <- FractionCountsInRegion(
    object = obj,
    assay = 'ATAC',
    regions = blacklist_mm10
)

## Call Peaks
peaks <- CallPeaks(
    obj, 
    macs2.path = "/opt/conda/bin/macs2",
    outdir = "callpeaks"
)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

## Quantify Counts in each Peak
macs2_counts <- FeatureMatrix(
    fragments = Fragments(obj),
    features = peaks,
    cells = colnames(obj)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
obj[["peaks"]] <- CreateChromatinAssay(
    counts = macs2_counts,
    fragments = "outs/fragments.tsv.gz",
    genome = "mm10"
)
Annotation(obj[["peaks"]]) <- annotations

# Perform Seurat Clustering for RNAseq data
DefaultAssay(obj) <- "RNA"
obj <- FindVariableFeatures(obj, nfeatures = 3000)
obj <- SCTransform(obj)
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:10, reduction.name = "umap.rna")
obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj, resolution = 0.8, algorithm = 3)
p1  <- DimPlot(obj, label = TRUE) + NoLegend() + ggtitle("RNA UMAP")
ggsave("RNA_UMAP.png", plot = p1, width = 6, height = 5, dpi = 300)
saveRDS(obj, file = "spatial.rds")

## Plotting cell images 
meta.data <- obj@meta.data
new_row_names <- row.names(meta.data)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data) <- new_row_names

print(head(meta.data))

object <- CreateSeuratObject(counts=rna_data, assay="Spatial", meta.data=meta.data)
image <- Read10X_Image(image.dir=tissue_dir, image.name="tissue_lowres_image.png", assay="Spatial", slice="slice1", filter.matrix = TRUE)
image <- image[Cells(x = object)]
DefaultAssay(image) <- "Spatial"
object[["slice1"]] <- image

spatial.obj <- object
col <- paletteer_d("ggthemes::Classic_20")[c(20,6,9,19,3,18,8,7,5,1,13,2,17,15,10,11,12,16)]
p2 <- DimPlot(data, label = TRUE,cols = col,group.by = "SCT_snn_res.0.8",pt.size = 1.8) + NoLegend() + ggtitle("RNA UMAP")
ggsave("RNA_UMAP_scn_res_0.8.png", plot = p1, width = 6, height = 5, dpi = 300)
p3 <- SpatialPlot(spatial.obj, label = FALSE, label.size = 1.8, group.by = 'SCT_snn_res.0.8', pt.size.factor = 1, cols = col, image.alpha = 1, stroke = 0)
p3$layers[[1]]$aes_params <- c(p3$layers[[1]]$aes_params, shape=22)
ggsave("Spatial_UMAP_scn_res_0.8.png", plot = p1, width = 6, height = 5, dpi = 300)