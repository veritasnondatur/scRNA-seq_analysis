#### Integrated analysis of hindlimb E10.5 and midface E9.5, E10.5 and E11.5 datasets

# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)

################################################################################
################## Pre-analysis of hindlimb_E10.5 dataset  #####################

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE10.5/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
hindlimb_E10.5 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_E10.5_hindlimb/data_zip')

# Create Seurat object
so_hindlimb_E10.5 <- CreateSeuratObject(counts = hindlimb_E10.5, 
                                        project = "hindlimb_E10.5", 
                                        min.cells = 3, min.features = 1000)  # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_hindlimb_E10.5 <- PercentageFeatureSet(so_hindlimb_E10.5, 
                                          pattern = "^mt-", 
                                          col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_hindlimb_E10.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "hindlimb_E10.5.qc.unfiltered.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

### Selecting cells for further analysis
percent.mt_max <- 5  # maximum percentage of mitochondrial genes (adjust as needed)
nFeature_RNA_min <- 100  # minimum number of features per cell
nFeature_RNA_max <- 5000 # minimum number of features per cell
nCount_RNA_min <- 100  # minimum number of RNA counts per cell
nCount_RNA_max <- 50000  # maximum number of RNA counts per cell

# Subset the Seurat object (filter based on thresholds above)
so_hindlimb_E10.5 <- subset(so_hindlimb_E10.5,
                            subset = percent.mt <= percent.mt_max &
                              nFeature_RNA >= nFeature_RNA_min &
                              nFeature_RNA <= nFeature_RNA_max &
                              nCount_RNA >= nCount_RNA_min &
                              nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_hindlimb_E10.5,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "hindlimb_E10.5.data.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_hindlimb_E10.5[["RNA"]]), invert = TRUE)
so_hindlimb_E10.5 <- CreateSeuratObject(counts = GetAssayData(so_hindlimb_E10.5[["RNA"]], layer = "counts")[keep, ], 
                                        project = "hindlimb_E10.5", 
                                        meta.data = so_hindlimb_E10.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_hindlimb_E10.5 <- SCTransform(so_hindlimb_E10.5,
                                 verbose = TRUE,
                                 variable.features.n = n_features)

############## Cell Cycle regression
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
so_hindlimb_E10.5 <- CellCycleScoring(so_hindlimb_E10.5,
                                      s.features = cc.genes_corrected$s.genes,   # S phase genes
                                      g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                      set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_hindlimb_E10.5 <- ScaleData(so_hindlimb_E10.5,
                               features = rownames(so_hindlimb_E10.5),
                               vars.to.regress = c("S.Score", "G2M.Score"),
                               verbose = TRUE
)

# PCA
so_hindlimb_E10.5 <- RunPCA(so_hindlimb_E10.5,
                            verbose = FALSE, 
                            npcs = 100,                # Number of principal components to compute
                            ndims.print = 1:5,         # Print details for the first 5 PCs
                            nfeatures.print = 30       # Print details for the top 30 features
)


# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_hindlimb_E10.5, ndims = 20)
out_file <- paste(output_folder, "hindlimb_E10.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 8                                      

#UMAP
so_hindlimb_E10.5 <- RunUMAP(so_hindlimb_E10.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_hindlimb_E10.5 <- FindNeighbors(so_hindlimb_E10.5, 
                                   dims = 1:pca_dim_sel)     # "Cell clusters were identified by constructing a shared nearest neighbor graph followed by a modularity optimization-based clustering algorithm (Leiden algorithm) using the top 60 principal components as determined by PCA." (Losa et al., 2023)
so_hindlimb_E10.5 <- FindClusters(so_hindlimb_E10.5, 
                                  resolution = 0.8,          # "Clustering was performed at multiple resolutions between 0.2 and 2, and optimal resolution was determined empirically based on the expression of known population markers (resolution = 0.8)." (Losa et al., 2023)
                                  algorithm = 4)             # Leiden clustering

# Visualize clusters as Dimplot
p <- DimPlot(object = so_hindlimb_E10.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/hindlimb_E10.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_hindlimb_E10.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/hindlimb_E10.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_hindlimb_E10.5, file = "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5_midfaceE9.5,E10.5,E11.5/pre-analysis_hindlimbE10.5/so_hindlimb_E10.5.rds")

# Load data
so_hindlimb_E10.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5_midfaceE9.5,E10.5,E11.5/pre-analysis_hindlimbE10.5/so_hindlimb_E10.5.rds")

FeaturePlot(so_hindlimb_E10.5, 
            features = "Runx2",
            #reduction = "umap.integrated",
            #split.by = "orig.ident",
            pt.size = 0.2)

################################################################################
################## Pre-analysis of hindlimb_E11.5 dataset  #####################

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE11.5/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
# Downloaded from GEO: GSM4227224
hindlimb_E11.5 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_mouse_hindlimb_publicDomain/GSM4227224_hindlimb_E11.5')

# Create Seurat object
so_hindlimb_E11.5 <- CreateSeuratObject(counts = hindlimb_E11.5, 
                                        project = "hindlimb_E11.5", 
                                        min.cells = 3, min.features = 1000)  # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_hindlimb_E11.5 <- PercentageFeatureSet(so_hindlimb_E11.5, 
                                          pattern = "^mt-", 
                                          col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_hindlimb_E11.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "hindlimb_E11.5.qc.unfiltered.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

### Selecting cells for further analysis
percent.mt_max <- 5  # maximum percentage of mitochondrial genes (adjust as needed)
nFeature_RNA_min <- 100  # minimum number of features per cell
nFeature_RNA_max <- 5000 # minimum number of features per cell
nCount_RNA_min <- 100  # minimum number of RNA counts per cell
nCount_RNA_max <- 30000  # maximum number of RNA counts per cell

# Subset the Seurat object (filter based on thresholds above)
so_hindlimb_E11.5 <- subset(so_hindlimb_E11.5,
                            subset = percent.mt <= percent.mt_max &
                              nFeature_RNA >= nFeature_RNA_min &
                              nFeature_RNA <= nFeature_RNA_max &
                              nCount_RNA >= nCount_RNA_min &
                              nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_hindlimb_E11.5,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "hindlimb_E11.5.data.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_hindlimb_E11.5[["RNA"]]), invert = TRUE)
so_hindlimb_E11.5 <- CreateSeuratObject(counts = GetAssayData(so_hindlimb_E11.5[["RNA"]], layer = "counts")[keep, ], 
                                        project = "hindlimb_E11.5", 
                                        meta.data = so_hindlimb_E11.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_hindlimb_E11.5 <- SCTransform(so_hindlimb_E11.5,
                                 verbose = TRUE,
                                 variable.features.n = n_features)

############## Cell Cycle regression
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
so_hindlimb_E11.5 <- CellCycleScoring(so_hindlimb_E11.5,
                                      s.features = cc.genes_corrected$s.genes,   # S phase genes
                                      g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                      set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_hindlimb_E11.5 <- ScaleData(so_hindlimb_E11.5,
                               features = rownames(so_hindlimb_E11.5),
                               vars.to.regress = c("S.Score", "G2M.Score"),
                               verbose = TRUE
)

# PCA
so_hindlimb_E11.5 <- RunPCA(so_hindlimb_E11.5,
                            verbose = FALSE, 
                            npcs = 100,                # Number of principal components to compute
                            ndims.print = 1:5,         # Print details for the first 5 PCs
                            nfeatures.print = 30       # Print details for the top 30 features
)


# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_hindlimb_E11.5, ndims = 20)
out_file <- paste(output_folder, "hindlimb_E11.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 10                                      

#UMAP
so_hindlimb_E11.5 <- RunUMAP(so_hindlimb_E11.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_hindlimb_E11.5 <- FindNeighbors(so_hindlimb_E11.5, 
                                   dims = 1:pca_dim_sel)     # "Cell clusters were identified by constructing a shared nearest neighbor graph followed by a modularity optimization-based clustering algorithm (Leiden algorithm) using the top 60 principal components as determined by PCA." (Losa et al., 2023)
so_hindlimb_E11.5 <- FindClusters(so_hindlimb_E11.5, 
                                  resolution = 0.8,          # "Clustering was performed at multiple resolutions between 0.2 and 2, and optimal resolution was determined empirically based on the expression of known population markers (resolution = 0.8)." (Losa et al., 2023)
                                  algorithm = 4)             # Leiden clustering

# Visualize clusters as Dimplot
p <- DimPlot(object = so_hindlimb_E11.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/hindlimb_E11.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_hindlimb_E11.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/hindlimb_E11.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_hindlimb_E11.5, file = "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE11.5/so_hindlimb_E11.5.rds")

# Load data
so_hindlimb_E11.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE11.5/so_hindlimb_E11.5.rds")

FeaturePlot(so_hindlimb_E11.5, 
            features = "Hand2",
            #reduction = "umap.integrated",
            #split.by = "orig.ident",
            pt.size = 0.2)

################################################################################
################## Pre-analysis of hindlimb_E13.5 dataset  #####################

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE13.5/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
hindlimb_E13.5 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_mouse_hindlimb_publicDomain/GSM4227225_hindlimb_E13.5')

# Create Seurat object
so_hindlimb_E13.5 <- CreateSeuratObject(counts = hindlimb_E13.5, 
                                        project = "hindlimb_E13.5", 
                                        min.cells = 3, min.features = 1000)  # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_hindlimb_E13.5 <- PercentageFeatureSet(so_hindlimb_E13.5, 
                                          pattern = "^mt-", 
                                          col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_hindlimb_E13.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "hindlimb_E13.5.qc.unfiltered.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

### Selecting cells for further analysis
percent.mt_max <- 5  # maximum percentage of mitochondrial genes (adjust as needed)
nFeature_RNA_min <- 1000  # minimum number of features per cell
nFeature_RNA_max <- 5000 # minimum number of features per cell
nCount_RNA_min <- 100  # minimum number of RNA counts per cell
nCount_RNA_max <- 50000  # maximum number of RNA counts per cell

# Subset the Seurat object (filter based on thresholds above)
so_hindlimb_E13.5 <- subset(so_hindlimb_E13.5,
                            subset = percent.mt <= percent.mt_max &
                              nFeature_RNA >= nFeature_RNA_min &
                              nFeature_RNA <= nFeature_RNA_max &
                              nCount_RNA >= nCount_RNA_min &
                              nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_hindlimb_E13.5,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "hindlimb_E13.5.data.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_hindlimb_E13.5[["RNA"]]), invert = TRUE)
so_hindlimb_E13.5 <- CreateSeuratObject(counts = GetAssayData(so_hindlimb_E13.5[["RNA"]], layer = "counts")[keep, ], 
                                        project = "hindlimb_E13.5", 
                                        meta.data = so_hindlimb_E13.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_hindlimb_E13.5 <- SCTransform(so_hindlimb_E13.5,
                                 verbose = TRUE,
                                 variable.features.n = n_features)

############## Cell Cycle regression
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
so_hindlimb_E13.5 <- CellCycleScoring(so_hindlimb_E13.5,
                                      s.features = cc.genes_corrected$s.genes,   # S phase genes
                                      g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                      set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_hindlimb_E13.5 <- ScaleData(so_hindlimb_E13.5,
                               features = rownames(so_hindlimb_E13.5),
                               vars.to.regress = c("S.Score", "G2M.Score"),
                               verbose = TRUE
)

# PCA
so_hindlimb_E13.5 <- RunPCA(so_hindlimb_E13.5,
                            verbose = FALSE, 
                            npcs = 100,                # Number of principal components to compute
                            ndims.print = 1:5,         # Print details for the first 5 PCs
                            nfeatures.print = 30       # Print details for the top 30 features
)


# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_hindlimb_E13.5, ndims = 20)
out_file <- paste(output_folder, "hindlimb_E13.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 9                                    

#UMAP
so_hindlimb_E13.5 <- RunUMAP(so_hindlimb_E13.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_hindlimb_E13.5 <- FindNeighbors(so_hindlimb_E13.5, 
                                   dims = 1:pca_dim_sel)     # "Cell clusters were identified by constructing a shared nearest neighbor graph followed by a modularity optimization-based clustering algorithm (Leiden algorithm) using the top 60 principal components as determined by PCA." (Losa et al., 2023)
so_hindlimb_E13.5 <- FindClusters(so_hindlimb_E13.5, 
                                  resolution = 0.8,          # "Clustering was performed at multiple resolutions between 0.2 and 2, and optimal resolution was determined empirically based on the expression of known population markers (resolution = 0.8)." (Losa et al., 2023)
                                  algorithm = 4)             # Leiden clustering

# Visualize clusters as Dimplot
p <- DimPlot(object = so_hindlimb_E13.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/hindlimb_E13.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_hindlimb_E13.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/hindlimb_E13.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_hindlimb_E13.5, file = "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE13.5/so_hindlimb_E13.5.rds")

# Load data
so_hindlimb_E13.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE13.5/so_hindlimb_E13.5.rds")

FeaturePlot(so_hindlimb_E13.5, 
            features = "Pbx1",
            #reduction = "umap.integrated",
            #split.by = "orig.ident",
            pt.size = 0.2)

################################################################################
################## Pre-analysis of hindlimb_E15.5 dataset  #####################

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE15.5/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
hindlimb_E15.5 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_mouse_hindlimb_publicDomain/GSM4227226_hindlimb_E15.5')

# Create Seurat object
so_hindlimb_E15.5 <- CreateSeuratObject(counts = hindlimb_E15.5, 
                                        project = "hindlimb_E15.5", 
                                        min.cells = 3, min.features = 1000)  # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_hindlimb_E15.5 <- PercentageFeatureSet(so_hindlimb_E15.5, 
                                          pattern = "^mt-", 
                                          col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_hindlimb_E15.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "hindlimb_E15.5.qc.unfiltered.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

### Selecting cells for further analysis
percent.mt_max <- 5  # maximum percentage of mitochondrial genes (adjust as needed)
nFeature_RNA_min <- 100  # minimum number of features per cell
nFeature_RNA_max <- 5000 # minimum number of features per cell
nCount_RNA_min <- 100  # minimum number of RNA counts per cell
nCount_RNA_max <- 50000  # maximum number of RNA counts per cell

# Subset the Seurat object (filter based on thresholds above)
so_hindlimb_E15.5 <- subset(so_hindlimb_E15.5,
                            subset = percent.mt <= percent.mt_max &
                              nFeature_RNA >= nFeature_RNA_min &
                              nFeature_RNA <= nFeature_RNA_max &
                              nCount_RNA >= nCount_RNA_min &
                              nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_hindlimb_E15.5,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "hindlimb_E15.5.data.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_hindlimb_E15.5[["RNA"]]), invert = TRUE)
so_hindlimb_E15.5 <- CreateSeuratObject(counts = GetAssayData(so_hindlimb_E15.5[["RNA"]], layer = "counts")[keep, ], 
                                        project = "hindlimb_E15.5", 
                                        meta.data = so_hindlimb_E15.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_hindlimb_E15.5 <- SCTransform(so_hindlimb_E15.5,
                                 verbose = TRUE,
                                 variable.features.n = n_features)

############## Cell Cycle regression
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
so_hindlimb_E15.5 <- CellCycleScoring(so_hindlimb_E15.5,
                                      s.features = cc.genes_corrected$s.genes,   # S phase genes
                                      g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                      set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_hindlimb_E15.5 <- ScaleData(so_hindlimb_E15.5,
                               features = rownames(so_hindlimb_E15.5),
                               vars.to.regress = c("S.Score", "G2M.Score"),
                               verbose = TRUE
)

# PCA
so_hindlimb_E15.5 <- RunPCA(so_hindlimb_E15.5,
                            verbose = FALSE, 
                            npcs = 100,                # Number of principal components to compute
                            ndims.print = 1:5,         # Print details for the first 5 PCs
                            nfeatures.print = 30       # Print details for the top 30 features
)


# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_hindlimb_E15.5, ndims = 20)
out_file <- paste(output_folder, "hindlimb_E15.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 11                                      

#UMAP
so_hindlimb_E15.5 <- RunUMAP(so_hindlimb_E15.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_hindlimb_E15.5 <- FindNeighbors(so_hindlimb_E15.5, 
                                   dims = 1:pca_dim_sel)     # "Cell clusters were identified by constructing a shared nearest neighbor graph followed by a modularity optimization-based clustering algorithm (Leiden algorithm) using the top 60 principal components as determined by PCA." (Losa et al., 2023)
so_hindlimb_E15.5 <- FindClusters(so_hindlimb_E15.5, 
                                  resolution = 0.8,          # "Clustering was performed at multiple resolutions between 0.2 and 2, and optimal resolution was determined empirically based on the expression of known population markers (resolution = 0.8)." (Losa et al., 2023)
                                  algorithm = 4)             # Leiden clustering

# Visualize clusters as Dimplot
p <- DimPlot(object = so_hindlimb_E15.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/hindlimb_E15.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_hindlimb_E15.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/hindlimb_E15.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_hindlimb_E15.5, file = "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE15.5/so_hindlimb_E15.5.rds")

# Load data
so_hindlimb_E15.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE15.5/so_hindlimb_E15.5.rds")

FeaturePlot(so_hindlimb_E15.5, 
            features = "Pbx1",
            #reduction = "umap.integrated",
            #split.by = "orig.ident",
            pt.size = 0.2)

################################################################################
################## Pre-analysis of hindlimb_E18.5 dataset  #####################

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE18.5/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
hindlimb_E18.5 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_mouse_hindlimb_publicDomain/GSM4227227_hindlimb_E18.5')

# Create Seurat object
so_hindlimb_E18.5 <- CreateSeuratObject(counts = hindlimb_E18.5, 
                                        project = "hindlimb_E18.5", 
                                        min.cells = 3, min.features = 1000)  # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_hindlimb_E18.5 <- PercentageFeatureSet(so_hindlimb_E18.5, 
                                          pattern = "^mt-", 
                                          col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_hindlimb_E18.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "hindlimb_E18.5.qc.unfiltered.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

### Selecting cells for further analysis
percent.mt_max <- 5  # maximum percentage of mitochondrial genes (adjust as needed)
nFeature_RNA_min <- 100  # minimum number of features per cell
nFeature_RNA_max <- 5000 # minimum number of features per cell
nCount_RNA_min <- 100  # minimum number of RNA counts per cell
nCount_RNA_max <- 50000  # maximum number of RNA counts per cell

# Subset the Seurat object (filter based on thresholds above)
so_hindlimb_E18.5 <- subset(so_hindlimb_E18.5,
                            subset = percent.mt <= percent.mt_max &
                              nFeature_RNA >= nFeature_RNA_min &
                              nFeature_RNA <= nFeature_RNA_max &
                              nCount_RNA >= nCount_RNA_min &
                              nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_hindlimb_E18.5,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "hindlimb_E18.5.data.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_hindlimb_E18.5[["RNA"]]), invert = TRUE)
so_hindlimb_E18.5 <- CreateSeuratObject(counts = GetAssayData(so_hindlimb_E18.5[["RNA"]], layer = "counts")[keep, ], 
                                        project = "hindlimb_E18.5", 
                                        meta.data = so_hindlimb_E18.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_hindlimb_E18.5 <- SCTransform(so_hindlimb_E18.5,
                                 verbose = TRUE,
                                 variable.features.n = n_features)

############## Cell Cycle regression
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
so_hindlimb_E18.5 <- CellCycleScoring(so_hindlimb_E18.5,
                                      s.features = cc.genes_corrected$s.genes,   # S phase genes
                                      g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                      set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_hindlimb_E18.5 <- ScaleData(so_hindlimb_E18.5,
                               features = rownames(so_hindlimb_E18.5),
                               vars.to.regress = c("S.Score", "G2M.Score"),
                               verbose = TRUE
)

# PCA
so_hindlimb_E18.5 <- RunPCA(so_hindlimb_E18.5,
                            verbose = FALSE, 
                            npcs = 100,                # Number of principal components to compute
                            ndims.print = 1:5,         # Print details for the first 5 PCs
                            nfeatures.print = 30       # Print details for the top 30 features
)


# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_hindlimb_E18.5, ndims = 20)
out_file <- paste(output_folder, "hindlimb_E18.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 11                                      

#UMAP
so_hindlimb_E18.5 <- RunUMAP(so_hindlimb_E18.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_hindlimb_E18.5 <- FindNeighbors(so_hindlimb_E18.5, 
                                   dims = 1:pca_dim_sel)     # "Cell clusters were identified by constructing a shared nearest neighbor graph followed by a modularity optimization-based clustering algorithm (Leiden algorithm) using the top 60 principal components as determined by PCA." (Losa et al., 2023)
so_hindlimb_E18.5 <- FindClusters(so_hindlimb_E18.5, 
                                  resolution = 0.8,          # "Clustering was performed at multiple resolutions between 0.2 and 2, and optimal resolution was determined empirically based on the expression of known population markers (resolution = 0.8)." (Losa et al., 2023)
                                  algorithm = 4)             # Leiden clustering

# Visualize clusters as Dimplot
p <- DimPlot(object = so_hindlimb_E18.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/hindlimb_E18.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_hindlimb_E18.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/hindlimb_E18.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_hindlimb_E18.5, file = "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE18.5/so_hindlimb_E18.5.rds")

# Load data
so_hindlimb_E18.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE18.5/so_hindlimb_E18.5.rds")

FeaturePlot(so_hindlimb_E18.5, 
            features = "Pbx1",
            #reduction = "umap.integrated",
            #split.by = "orig.ident",
            pt.size = 0.2)


################################################################################
################## Pre-analysis of midface_E9.5 dataset  #####################

# Set out folder (to store results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre_analysis_midfaceE9.5/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
midface_E9.5 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/MF-E9.5_ML11_071719/filtered_feature_bc_matrix')

# Create Seurat object
so_midface_E9.5 <- CreateSeuratObject(counts = midface_E9.5, 
                                      project = "midface_E9.5", 
                                      min.cells = 3, min.features = 1000)      # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_midface_E9.5 <- PercentageFeatureSet(so_midface_E9.5, 
                                        pattern = "^mt-", 
                                        col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_midface_E9.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "midface_E9.5.qc.unfiltered.pdf", sep = "")
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
so_midface_E9.5 <- subset(so_midface_E9.5,
                          subset = percent.mt <= percent.mt_max &
                            nFeature_RNA >= nFeature_RNA_min &
                            nFeature_RNA <= nFeature_RNA_max &
                            nCount_RNA >= nCount_RNA_min &
                            nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_midface_E9.5,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "midface_E9.5.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_midface_E9.5[["RNA"]]), invert = TRUE)
so_midface_E9.5 <- CreateSeuratObject(counts = GetAssayData(so_midface_E9.5[["RNA"]], layer = "counts")[keep, ], 
                                      project = "midface_E9.5", 
                                      meta.data = so_midface_E9.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_midface_E9.5 <- SCTransform(so_midface_E9.5,
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
so_midface_E9.5 <- CellCycleScoring(so_midface_E9.5,
                                    s.features = cc.genes_corrected$s.genes,   # S phase genes
                                    g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                    set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_midface_E9.5 <- ScaleData(so_midface_E9.5,
                             features = rownames(so_midface_E9.5),
                             vars.to.regress = c("S.Score", "G2M.Score"),
                             verbose = TRUE
)

# PCA
so_midface_E9.5 <- RunPCA(so_midface_E9.5,
                          verbose = FALSE, 
                          npcs = 100,                # Number of principal components to compute
                          ndims.print = 1:5,         # Print details for the first 5 PCs
                          nfeatures.print = 30       # Print details for the top 30 features
)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_midface_E9.5, ndims = 100)
out_file <- paste(output_folder, "midface_E9.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <-8   

#UMAP
so_midface_E9.5 <- RunUMAP(so_midface_E9.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_midface_E9.5 <- FindNeighbors(so_midface_E9.5, dims = 1:pca_dim_sel)         
so_midface_E9.5 <- FindClusters(so_midface_E9.5, 
                                resolution = 0.8,
                                algorithm = 4)            

# Visualize clusters as Dimplot
p <- DimPlot(object = so_midface_E9.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/midface_E9.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_midface_E9.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/midface_E9.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_midface_E9.5, file = "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5_midfaceE9.5,E10.5,E11.5/pre_analysis_midfaceE9.5/so_midface_E9.5.rds")
