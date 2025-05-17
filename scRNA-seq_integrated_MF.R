#### Public Domain data exploration of midface tissues

# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)


################################################################################
################### Pre-analysis of anterior_palate_E13.5_rep1 dataset  ######################

# Set out folder (to store results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_midface_mandible_palate/pre-analysis_anterior_palate_E13.5_rep1/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
anterior_palate_E13.5_rep1 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/anterior_palate_E13.5_FB00001364_rep1')

# Create Seurat object
so_anterior_palate_E13.5_rep1 <- CreateSeuratObject(counts = anterior_palate_E13.5_rep1, 
                                       project = "anterior_palate_E13.5_rep1", 
                                       min.cells = 3, min.features = 1000)      # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_anterior_palate_E13.5_rep1 <- PercentageFeatureSet(so_anterior_palate_E13.5_rep1, 
                                         pattern = "^mt-", 
                                         col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_anterior_palate_E13.5_rep1, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "anterior_palate_E13.5_rep1.qc.unfiltered.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

### Selecting cells for further analysis
percent.mt_max <- 5          # maximum percentage of mitochondrial genes (adjust as needed)
nFeature_RNA_min <- 100      # minimum number of features per cell
nFeature_RNA_max <- 5000     # minimum number of features per cell
nCount_RNA_min <- 100        # minimum number of RNA counts per cell
nCount_RNA_max <- 50000      # maximum number of RNA counts per cell

# Subset the Seurat object (filter based on thresholds above)
so_anterior_palate_E13.5_rep1 <- subset(so_anterior_palate_E13.5_rep1,
                           subset = percent.mt <= percent.mt_max &
                             nFeature_RNA >= nFeature_RNA_min &
                             nFeature_RNA <= nFeature_RNA_max &
                             nCount_RNA >= nCount_RNA_min &
                             nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_anterior_palate_E13.5_rep1,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "anterior_palate_E13.5_rep1.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_anterior_palate_E13.5_rep1[["RNA"]]), invert = TRUE)
so_anterior_palate_E13.5_rep1 <- CreateSeuratObject(counts = GetAssayData(so_anterior_palate_E13.5_rep1[["RNA"]], layer = "counts")[keep, ], 
                                       project = "anterior_palate_E13.5_rep1", 
                                       meta.data = so_anterior_palate_E13.5_rep1@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_anterior_palate_E13.5_rep1 <- SCTransform(so_anterior_palate_E13.5_rep1,
                                verbose = TRUE,
                                variable.features.n = n_features)

############## Cell Cycle regression (Script can be run from here)
# 1. Use Seurat's predefined cell cycle genes
cc.genes <- Seurat::cc.genes

# 2. Perform cell cycle scoring
# Function to capitalize only the first letter and make the rest lowercase
capitalize_first_letter <- function(x) {
  # Make the entire string lowercase, then capitalize the first letter
  return(paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2))))
}

# Apply this function to the Seurat predefined cell cycle genes
cc.genes_corrected <- list(
  s.genes = sapply(cc.genes$s.genes, capitalize_first_letter),   # S phase genes
  g2m.genes = sapply(cc.genes$g2m.genes, capitalize_first_letter)  # G2M phase genes
)

# Perform cell cycle scoring using the corrected gene sets
so_anterior_palate_E13.5_rep1 <- CellCycleScoring(so_anterior_palate_E13.5_rep1,
                                     s.features = cc.genes_corrected$s.genes,   # S phase genes
                                     g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                     set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_anterior_palate_E13.5_rep1 <- ScaleData(so_anterior_palate_E13.5_rep1,
                              features = rownames(so_anterior_palate_E13.5_rep1),
                              vars.to.regress = c("S.Score", "G2M.Score"),
                              verbose = TRUE
)

# PCA
so_anterior_palate_E13.5_rep1 <- RunPCA(so_anterior_palate_E13.5_rep1,
                           verbose = FALSE, 
                           npcs = 100,                # Number of principal components to compute
                           ndims.print = 1:5,         # Print details for the first 5 PCs
                           nfeatures.print = 30       # Print details for the top 30 features
)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_anterior_palate_E13.5_rep1, ndims = 20)
out_file <- paste(output_folder, "anterior_palate_E13.5_rep1.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 17   

#UMAP
so_anterior_palate_E13.5_rep1 <- RunUMAP(so_anterior_palate_E13.5_rep1, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_anterior_palate_E13.5_rep1 <- FindNeighbors(so_anterior_palate_E13.5_rep1, dims = 1:pca_dim_sel)         
so_anterior_palate_E13.5_rep1 <- FindClusters(so_anterior_palate_E13.5_rep1, 
                                 resolution = 0.5,
                                 algorithm = 4)            

# Visualize clusters as Dimplot
p <- DimPlot(object = so_anterior_palate_E13.5_rep1,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/anterior_palate_E13.5_rep1.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_anterior_palate_E13.5_rep1, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/anterior_palate_E13.5_rep1.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_anterior_palate_E13.5_rep1, file = "~/Documents/postdoc/bioinformatics/results/integrated_midface_mandible_palate/pre-analysis_anterior_palate_E13.5_rep1/so_anterior_palate_E13.5_rep1.rds")

# Load data
so_anterior_palate_E13.5_rep1 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_midface_mandible_palate/pre-analysis_anterior_palate_E13.5_rep1/so_anterior_palate_E13.5_rep1.rds")

FeaturePlot(so_anterior_palate_E13.5_rep1, 
            features = "Col2a1",
            pt.size = 0.2)


################################################################################
###################### INTEGRATION OF DATASETS ON WYNTON #######################

# Define output folder (for results)
output_folder <- "/wynton/home/selleri/veritasnondatur/scRNA-seq/midface_integration/analysis/"

# Load raw data (pre-processed Seurat objects produced in previous analysis)
so_midface_E9.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre_analysis_midfaceE9.5/so_midface_E9.5.rds")
so_midface_E10.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre_analysis_midfaceE10.5/so_midface_E10.5.rds")
so_midface_E11.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre_analysis_midfaceE11.5/so_midface_E11.5.rds")

so_mandible_E9.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE9.5/so_mandible_E9.5.rds")
so_mandible_E10.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE10.5/so_mandible_E10.5.rds")
so_mandible_E11.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE11.5/so_mandible_E11.5.rds")
