#### Public Domain data exploration of midface tissues

# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)


################################################################################
############# Pre-analysis of anterior_palate_E13.5_rep1 dataset  ##############

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
            features = "Vps25",
            pt.size = 0.2)

################################################################################
############ Pre-analysis of anterior_palate_E13.5_rep2 dataset  ###############

# Set out folder (to store results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_midface_mandible_palate/pre-analysis_anterior_palate_E13.5_rep2/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
anterior_palate_E13.5_rep2 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/anterior_palate_E13.5_FB00001364_rep2')

# Create Seurat object
so_anterior_palate_E13.5_rep2 <- CreateSeuratObject(counts = anterior_palate_E13.5_rep2, 
                                                    project = "anterior_palate_E13.5_rep2", 
                                                    min.cells = 3, min.features = 1000)      # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_anterior_palate_E13.5_rep2 <- PercentageFeatureSet(so_anterior_palate_E13.5_rep2, 
                                                      pattern = "^mt-", 
                                                      col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_anterior_palate_E13.5_rep2, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "anterior_palate_E13.5_rep2.qc.unfiltered.pdf", sep = "")
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
so_anterior_palate_E13.5_rep2 <- subset(so_anterior_palate_E13.5_rep2,
                                        subset = percent.mt <= percent.mt_max &
                                          nFeature_RNA >= nFeature_RNA_min &
                                          nFeature_RNA <= nFeature_RNA_max &
                                          nCount_RNA >= nCount_RNA_min &
                                          nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_anterior_palate_E13.5_rep2,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "anterior_palate_E13.5_rep2.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_anterior_palate_E13.5_rep2[["RNA"]]), invert = TRUE)
so_anterior_palate_E13.5_rep2 <- CreateSeuratObject(counts = GetAssayData(so_anterior_palate_E13.5_rep2[["RNA"]], layer = "counts")[keep, ], 
                                                    project = "anterior_palate_E13.5_rep2", 
                                                    meta.data = so_anterior_palate_E13.5_rep2@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_anterior_palate_E13.5_rep2 <- SCTransform(so_anterior_palate_E13.5_rep2,
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
so_anterior_palate_E13.5_rep2 <- CellCycleScoring(so_anterior_palate_E13.5_rep2,
                                                  s.features = cc.genes_corrected$s.genes,   # S phase genes
                                                  g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                                  set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_anterior_palate_E13.5_rep2 <- ScaleData(so_anterior_palate_E13.5_rep2,
                                           features = rownames(so_anterior_palate_E13.5_rep2),
                                           vars.to.regress = c("S.Score", "G2M.Score"),
                                           verbose = TRUE
)

# PCA
so_anterior_palate_E13.5_rep2 <- RunPCA(so_anterior_palate_E13.5_rep2,
                                        verbose = FALSE, 
                                        npcs = 100,                # Number of principal components to compute
                                        ndims.print = 1:5,         # Print details for the first 5 PCs
                                        nfeatures.print = 30       # Print details for the top 30 features
)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_anterior_palate_E13.5_rep2, ndims = 20)
out_file <- paste(output_folder, "anterior_palate_E13.5_rep2.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 15  

#UMAP
so_anterior_palate_E13.5_rep2 <- RunUMAP(so_anterior_palate_E13.5_rep2, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_anterior_palate_E13.5_rep2 <- FindNeighbors(so_anterior_palate_E13.5_rep2, dims = 1:pca_dim_sel)         
so_anterior_palate_E13.5_rep2 <- FindClusters(so_anterior_palate_E13.5_rep2, 
                                              resolution = 0.5,
                                              algorithm = 4)            

# Visualize clusters as Dimplot
p <- DimPlot(object = so_anterior_palate_E13.5_rep2,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/anterior_palate_E13.5_rep2.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_anterior_palate_E13.5_rep2, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/anterior_palate_E13.5_rep2.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_anterior_palate_E13.5_rep2, file = "~/Documents/postdoc/bioinformatics/results/integrated_midface_mandible_palate/pre-analysis_anterior_palate_E13.5_rep2/so_anterior_palate_E13.5_rep2.rds")

# Load data
so_anterior_palate_E13.5_rep2 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_midface_mandible_palate/pre-analysis_anterior_palate_E13.5_rep2/so_anterior_palate_E13.5_rep2.rds")

FeaturePlot(so_anterior_palate_E13.5_rep2, 
            features = "Vps25",
            pt.size = 0.2)

################################################################################
############ Pre-analysis of anterior_palate_E13.5_rep3 dataset  ###############

# Set out folder (to store results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_midface_mandible_palate/pre-analysis_anterior_palate_E13.5_rep3/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
anterior_palate_E13.5_rep3 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/anterior_palate_E13.5_FB00001364_rep3')

# Create Seurat object
so_anterior_palate_E13.5_rep3 <- CreateSeuratObject(counts = anterior_palate_E13.5_rep3, 
                                                    project = "anterior_palate_E13.5_rep3", 
                                                    min.cells = 3, min.features = 1000)      # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_anterior_palate_E13.5_rep3 <- PercentageFeatureSet(so_anterior_palate_E13.5_rep3, 
                                                      pattern = "^mt-", 
                                                      col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_anterior_palate_E13.5_rep3, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "anterior_palate_E13.5_rep3.qc.unfiltered.pdf", sep = "")
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
so_anterior_palate_E13.5_rep3 <- subset(so_anterior_palate_E13.5_rep3,
                                        subset = percent.mt <= percent.mt_max &
                                          nFeature_RNA >= nFeature_RNA_min &
                                          nFeature_RNA <= nFeature_RNA_max &
                                          nCount_RNA >= nCount_RNA_min &
                                          nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_anterior_palate_E13.5_rep3,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "anterior_palate_E13.5_rep3.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_anterior_palate_E13.5_rep3[["RNA"]]), invert = TRUE)
so_anterior_palate_E13.5_rep3 <- CreateSeuratObject(counts = GetAssayData(so_anterior_palate_E13.5_rep3[["RNA"]], layer = "counts")[keep, ], 
                                                    project = "anterior_palate_E13.5_rep3", 
                                                    meta.data = so_anterior_palate_E13.5_rep3@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_anterior_palate_E13.5_rep3 <- SCTransform(so_anterior_palate_E13.5_rep3,
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
so_anterior_palate_E13.5_rep3 <- CellCycleScoring(so_anterior_palate_E13.5_rep3,
                                                  s.features = cc.genes_corrected$s.genes,   # S phase genes
                                                  g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                                  set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_anterior_palate_E13.5_rep3 <- ScaleData(so_anterior_palate_E13.5_rep3,
                                           features = rownames(so_anterior_palate_E13.5_rep3),
                                           vars.to.regress = c("S.Score", "G2M.Score"),
                                           verbose = TRUE
)

# PCA
so_anterior_palate_E13.5_rep3 <- RunPCA(so_anterior_palate_E13.5_rep3,
                                        verbose = FALSE, 
                                        npcs = 100,                # Number of principal components to compute
                                        ndims.print = 1:5,         # Print details for the first 5 PCs
                                        nfeatures.print = 30       # Print details for the top 30 features
)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_anterior_palate_E13.5_rep3, ndims = 20)
out_file <- paste(output_folder, "anterior_palate_E13.5_rep3.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 13   

#UMAP
so_anterior_palate_E13.5_rep3 <- RunUMAP(so_anterior_palate_E13.5_rep3, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_anterior_palate_E13.5_rep3 <- FindNeighbors(so_anterior_palate_E13.5_rep3, dims = 1:pca_dim_sel)         
so_anterior_palate_E13.5_rep3 <- FindClusters(so_anterior_palate_E13.5_rep3, 
                                              resolution = 0.5,
                                              algorithm = 4)            

# Visualize clusters as Dimplot
p <- DimPlot(object = so_anterior_palate_E13.5_rep3,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/anterior_palate_E13.5_rep3.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_anterior_palate_E13.5_rep3, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/anterior_palate_E13.5_rep3.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_anterior_palate_E13.5_rep3, file = "~/Documents/postdoc/bioinformatics/results/integrated_midface_mandible_palate/pre-analysis_anterior_palate_E13.5_rep3/so_anterior_palate_E13.5_rep3.rds")

# Load data
so_anterior_palate_E13.5_rep3 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_midface_mandible_palate/pre-analysis_anterior_palate_E13.5_rep3/so_anterior_palate_E13.5_rep3.rds")

FeaturePlot(so_anterior_palate_E13.5_rep3, 
            features = "Vps25",
            pt.size = 0.2)

################################################################################
############## Pre-analysis of secondary_palate_E14.5 dataset  #################

# Set out folder (to store results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_midface_mandible_palate/pre-analysis_secondary_palate_E14.5/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
secondary_palate_E14.5 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/secondary_palate_E14.5_1-YKAG_FB00001124')

# Create Seurat object
so_secondary_palate_E14.5 <- CreateSeuratObject(counts = secondary_palate_E14.5, 
                                                project = "secondary_palate_E14.5", 
                                                min.cells = 3, min.features = 1000)      # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_secondary_palate_E14.5 <- PercentageFeatureSet(so_secondary_palate_E14.5, 
                                                  pattern = "^mt-", 
                                                  col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_secondary_palate_E14.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "secondary_palate_E14.5.qc.unfiltered.pdf", sep = "")
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
so_secondary_palate_E14.5 <- subset(so_secondary_palate_E14.5,
                                    subset = percent.mt <= percent.mt_max &
                                      nFeature_RNA >= nFeature_RNA_min &
                                      nFeature_RNA <= nFeature_RNA_max &
                                      nCount_RNA >= nCount_RNA_min &
                                      nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_secondary_palate_E14.5,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "secondary_palate_E14.5.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_secondary_palate_E14.5[["RNA"]]), invert = TRUE)
so_secondary_palate_E14.5 <- CreateSeuratObject(counts = GetAssayData(so_secondary_palate_E14.5[["RNA"]], layer = "counts")[keep, ], 
                                                project = "secondary_palate_E14.5", 
                                                meta.data = so_secondary_palate_E14.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_secondary_palate_E14.5 <- SCTransform(so_secondary_palate_E14.5,
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
so_secondary_palate_E14.5 <- CellCycleScoring(so_secondary_palate_E14.5,
                                              s.features = cc.genes_corrected$s.genes,   # S phase genes
                                              g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                              set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_secondary_palate_E14.5 <- ScaleData(so_secondary_palate_E14.5,
                                       features = rownames(so_secondary_palate_E14.5),
                                       vars.to.regress = c("S.Score", "G2M.Score"),
                                       verbose = TRUE
)

# PCA
so_secondary_palate_E14.5 <- RunPCA(so_secondary_palate_E14.5,
                                    verbose = FALSE, 
                                    npcs = 100,                # Number of principal components to compute
                                    ndims.print = 1:5,         # Print details for the first 5 PCs
                                    nfeatures.print = 30       # Print details for the top 30 features
)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_secondary_palate_E14.5, ndims = 20)
out_file <- paste(output_folder, "secondary_palate_E14.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 14  

#UMAP
so_secondary_palate_E14.5 <- RunUMAP(so_secondary_palate_E14.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_secondary_palate_E14.5 <- FindNeighbors(so_secondary_palate_E14.5, dims = 1:pca_dim_sel)         
so_secondary_palate_E14.5 <- FindClusters(so_secondary_palate_E14.5, 
                                          resolution = 0.5,
                                          algorithm = 4)            

# Visualize clusters as Dimplot
p <- DimPlot(object = so_secondary_palate_E14.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/secondary_palate_E14.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_secondary_palate_E14.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/secondary_palate_E14.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_secondary_palate_E14.5, file = "~/Documents/postdoc/bioinformatics/results/integrated_midface_mandible_palate/pre-analysis_secondary_palate_E14.5/so_secondary_palate_E14.5.rds")

# Load data
so_secondary_palate_E14.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_midface_mandible_palate/pre-analysis_secondary_palate_E14.5/so_secondary_palate_E14.5.rds")

FeaturePlot(so_secondary_palate_E14.5, 
            features = "Vps25",
            pt.size = 0.2)


################################################################################
############## Mouse facial tissue E10.5-E15.5 dataset  #################

# Load data
so_facial_tissue_E10.5_E15.5 <- readRDS("/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/mouse_facial_tissue_scRNA-seq-E11.5_E12.5_E13.5_snATAC-seq-E10.5-E15.5_FB00001359/2021-05-23_seurat.face.rds")

DefaultAssay(so_facial_tissue_E10.5_E15.5) <- "SCT"

FeaturePlot(so_facial_tissue_E10.5_E15.5, 
            features = "Hand2",
            split.by = 'orig.ident',
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

### Gene list to explore
goi <- c("Runx2", "Alpl","Bglap", "Col1a1", "Sp7", "Ibsp", "Spp1", "Mmp13",     # osteoblast markers (according to ChatGPT, so careful!)
         "Sox9", "Col2a1", "Acan", "Comp", "Col11a1", "Cdh2", "Ihh", "Mgp"      # chondroblast markers (according to ChatGPT, so careful!)
        )