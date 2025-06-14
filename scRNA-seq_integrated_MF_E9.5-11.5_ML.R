#### Integrated analysis of midface E9.5, E10.5 and E11.5 datasets

# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)


################################################################################
################### Pre-analysis of midface_E9.5 dataset  ######################

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
p <- ElbowPlot(so_midface_E9.5, ndims = 20)
out_file <- paste(output_folder, "midface_E9.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 15   

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
saveRDS(so_midface_E9.5, file = "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre_analysis_midfaceE9.5/so_midface_E9.5.rds")

# Load data
so_midface_E9.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre_analysis_midfaceE9.5/so_midface_E9.5.rds")

FeaturePlot(so_midface_E9.5, 
            features = "Pbx1",
            #reduction = "umap.integrated",
            #split.by = "orig.ident",
            pt.size = 0.2)


################################################################################
################## Pre-analysis of midface_E10.5 dataset  #####################

# Set out folder (to store results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre_analysis_midfaceE10.5/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
midface_E10.5 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/MF-E10.5_ML10_071719/filtered_feature_bc_matrix')

# Create Seurat object
so_midface_E10.5 <- CreateSeuratObject(counts = midface_E10.5, 
                                       project = "midface_E10.5", 
                                       min.cells = 3, min.features = 1000)      # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_midface_E10.5 <- PercentageFeatureSet(so_midface_E10.5, 
                                         pattern = "^mt-", 
                                         col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_midface_E10.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "midface_E10.5.qc.unfiltered.pdf", sep = "")
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
so_midface_E10.5 <- subset(so_midface_E10.5,
                           subset = percent.mt <= percent.mt_max &
                             nFeature_RNA >= nFeature_RNA_min &
                             nFeature_RNA <= nFeature_RNA_max &
                             nCount_RNA >= nCount_RNA_min &
                             nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_midface_E10.5,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "midface_E10.5.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_midface_E10.5[["RNA"]]), invert = TRUE)
so_midface_E10.5 <- CreateSeuratObject(counts = GetAssayData(so_midface_E10.5[["RNA"]], layer = "counts")[keep, ], 
                                       project = "midface_E10.5", 
                                       meta.data = so_midface_E10.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_midface_E10.5 <- SCTransform(so_midface_E10.5,
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
so_midface_E10.5 <- CellCycleScoring(so_midface_E10.5,
                                     s.features = cc.genes_corrected$s.genes,   # S phase genes
                                     g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                     set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_midface_E10.5 <- ScaleData(so_midface_E10.5,
                              features = rownames(so_midface_E10.5),
                              vars.to.regress = c("S.Score", "G2M.Score"),
                              verbose = TRUE
)

# PCA
so_midface_E10.5 <- RunPCA(so_midface_E10.5,
                           verbose = FALSE, 
                           npcs = 100,                # Number of principal components to compute
                           ndims.print = 1:5,         # Print details for the first 5 PCs
                           nfeatures.print = 30       # Print details for the top 30 features
)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_midface_E10.5, ndims = 20)
out_file <- paste(output_folder, "midface_E10.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 10   

#UMAP
so_midface_E10.5 <- RunUMAP(so_midface_E10.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_midface_E10.5 <- FindNeighbors(so_midface_E10.5, dims = 1:pca_dim_sel)         
so_midface_E10.5 <- FindClusters(so_midface_E10.5, 
                                 resolution = 0.8,
                                 algorithm = 4)            

# Visualize clusters as Dimplot
p <- DimPlot(object = so_midface_E10.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/midface_E10.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_midface_E10.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/midface_E10.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_midface_E10.5, file = "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre_analysis_midfaceE10.5/so_midface_E10.5.rds")

# Load data
so_midface_E10.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre_analysis_midfaceE10.5/so_midface_E10.5.rds")

FeaturePlot(so_midface_E10.5, 
            features = "Pbx1",
            #reduction = "umap.integrated",
            #split.by = "orig.ident",
            pt.size = 0.2)


################################################################################
################## Pre-analysis of midface_E11.5 dataset  #####################

# Set out folder (to store results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre_analysis_midfaceE11.5/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
midface_E11.5 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/MF-E11.5_rep1_ML8_071719/filtered_feature_bc_matrix')

# Create Seurat object
so_midface_E11.5 <- CreateSeuratObject(counts = midface_E11.5, 
                                       project = "midface_E11.5", 
                                       min.cells = 3, min.features = 1000)      # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_midface_E11.5 <- PercentageFeatureSet(so_midface_E11.5, 
                                         pattern = "^mt-", 
                                         col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_midface_E11.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "midface_E11.5.qc.unfiltered.pdf", sep = "")
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
so_midface_E11.5 <- subset(so_midface_E11.5,
                           subset = percent.mt <= percent.mt_max &
                             nFeature_RNA >= nFeature_RNA_min &
                             nFeature_RNA <= nFeature_RNA_max &
                             nCount_RNA >= nCount_RNA_min &
                             nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_midface_E11.5,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "midface_E11.5.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_midface_E11.5[["RNA"]]), invert = TRUE)
so_midface_E11.5 <- CreateSeuratObject(counts = GetAssayData(so_midface_E11.5[["RNA"]], layer = "counts")[keep, ], 
                                       project = "midface_E11.5", 
                                       meta.data = so_midface_E11.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_midface_E11.5 <- SCTransform(so_midface_E11.5,
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
so_midface_E11.5 <- CellCycleScoring(so_midface_E11.5,
                                     s.features = cc.genes_corrected$s.genes,   # S phase genes
                                     g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                     set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_midface_E11.5 <- ScaleData(so_midface_E11.5,
                              features = rownames(so_midface_E11.5),
                              vars.to.regress = c("S.Score", "G2M.Score"),
                              verbose = TRUE
)

# PCA
so_midface_E11.5 <- RunPCA(so_midface_E11.5,
                           verbose = FALSE, 
                           npcs = 100,                # Number of principal components to compute
                           ndims.print = 1:5,         # Print details for the first 5 PCs
                           nfeatures.print = 30       # Print details for the top 30 features
)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_midface_E11.5, ndims = 20)
out_file <- paste(output_folder, "midface_E11.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 8   

#UMAP
so_midface_E11.5 <- RunUMAP(so_midface_E11.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_midface_E11.5 <- FindNeighbors(so_midface_E11.5, dims = 1:pca_dim_sel)         
so_midface_E11.5 <- FindClusters(so_midface_E11.5, 
                                 resolution = 0.8,
                                 algorithm = 4)            

# Visualize clusters as Dimplot
p <- DimPlot(object = so_midface_E11.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/midface_E11.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_midface_E11.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/midface_E11.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_midface_E11.5, file = "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre_analysis_midfaceE11.5/so_midface_E11.5.rds")

# Load data
so_midface_E11.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre_analysis_midfaceE11.5/so_midface_E11.5.rds")

FeaturePlot(so_midface_E11.5, 
            features = "Pbx1",
            #reduction = "umap.integrated",
            #split.by = "orig.ident",
            pt.size = 0.2)


######################## Plots for pre-analysis datasets #######################

##### midface E9.5
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre_analysis_midfaceE9.5/"

# UMAP of goi
out_file <- paste(output_folder, "/midface_E9.5.UMAP.goi.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi) {
  if (gene %in% rownames(so_midface_E9.5[["SCT"]])) {
    p <- FeaturePlot(so_midface_E9.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi
p <- FeaturePlot(so_midface_E9.5, features = goi,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/midface_E9.5.UMAP.goi.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()

# UMAP of goi_celltype
out_file <- paste(output_folder, "/midface_E9.5.UMAP.goi_celltype.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi_celltype) {
  if (gene %in% rownames(so_midface_E9.5[["SCT"]])) {
    p <- FeaturePlot(so_midface_E9.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi_celltype
p <- FeaturePlot(so_midface_E9.5, features = goi_celltype,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/midface_E9.5.UMAP.goi_celltype.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()

# UMAP of goi_spatial
out_file <- paste(output_folder, "/midface_E9.5.UMAP.goi_spatial.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi_spatial) {
  if (gene %in% rownames(so_midface_E9.5[["SCT"]])) {
    p <- FeaturePlot(so_midface_E9.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi_spatial
p <- FeaturePlot(so_midface_E9.5, features = goi_spatial,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/midface_E9.5.UMAP.goi_spatial.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()


##### midface E10.5
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre_analysis_midfaceE10.5/"

# UMAP of goi
out_file <- paste(output_folder, "/midface_E10.5.UMAP.goi.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi) {
  if (gene %in% rownames(so_midface_E10.5[["SCT"]])) {
    p <- FeaturePlot(so_midface_E10.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi
p <- FeaturePlot(so_midface_E10.5, features = goi,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/midface_E10.5.UMAP.goi.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()

# UMAP of goi_celltype
out_file <- paste(output_folder, "/midface_E10.5.UMAP.goi_celltype.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi_celltype) {
  if (gene %in% rownames(so_midface_E10.5[["SCT"]])) {
    p <- FeaturePlot(so_midface_E10.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi_celltype
p <- FeaturePlot(so_midface_E10.5, features = goi_celltype,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/midface_E10.5.UMAP.goi_celltype.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()

# UMAP of goi_spatial
out_file <- paste(output_folder, "/midface_E10.5.UMAP.goi_spatial.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi_spatial) {
  if (gene %in% rownames(so_midface_E10.5[["SCT"]])) {
    p <- FeaturePlot(so_midface_E10.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi_spatial
p <- FeaturePlot(so_midface_E10.5, features = goi_spatial,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/midface_E10.5.UMAP.goi_spatial.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()


##### midface E11.5
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre_analysis_midfaceE11.5/"

# UMAP of goi
out_file <- paste(output_folder, "/midface_E11.5.UMAP.goi.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi) {
  if (gene %in% rownames(so_midface_E11.5[["SCT"]])) {
    p <- FeaturePlot(so_midface_E11.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi
p <- FeaturePlot(so_midface_E11.5, features = goi,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/midface_E11.5.UMAP.goi.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()

# UMAP of goi_celltype
out_file <- paste(output_folder, "/midface_E11.5.UMAP.goi_celltype.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi_celltype) {
  if (gene %in% rownames(so_midface_E11.5[["SCT"]])) {
    p <- FeaturePlot(so_midface_E11.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi_celltype
p <- FeaturePlot(so_midface_E11.5, features = goi_celltype,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/midface_E11.5.UMAP.goi_celltype.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()

# UMAP of goi_spatial
out_file <- paste(output_folder, "/midface_E11.5.UMAP.goi_spatial.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi_spatial) {
  if (gene %in% rownames(so_midface_E11.5[["SCT"]])) {
    p <- FeaturePlot(so_midface_E11.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi_spatial
p <- FeaturePlot(so_midface_E11.5, features = goi_spatial,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/midface_E11.5.UMAP.goi_spatial.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()


################################################################################
####################### INTEGRATION OF MIDFACE DATASETS ########################

# Define output folder (for results)
output_folder <- "/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_midfaceE9.5-11.5/analysis/"

# Load raw data (pre-processed Seurat objects produced above)
so_midface_E9.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_midfaceE9.5-11.5/data/pre_analysis_midfaceE9.5/so_midface_E9.5.rds")
so_midface_E10.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_midfaceE9.5-11.5/data/pre_analysis_midfaceE10.5/so_midface_E10.5.rds")
so_midface_E11.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_midfaceE9.5-11.5/data/pre_analysis_midfaceE11.5/so_midface_E11.5.rds")

# Change the default assay to "SCT"
DefaultAssay(so_midface_E9.5) <- "SCT"
DefaultAssay(so_midface_E10.5) <- "SCT"
DefaultAssay(so_midface_E11.5) <- "SCT"

# Check memory usage before and after merge (due to Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Identify variable features for both datasets
so_midface_E9.5 <- FindVariableFeatures(so_midface_E9.5, 
                                        selection.method = "vst", 
                                        nfeatures = 2000)
so_midface_E10.5 <- FindVariableFeatures(so_midface_E10.5, 
                                         selection.method = "vst", 
                                         nfeatures = 2000)
so_midface_E11.5 <- FindVariableFeatures(so_midface_E11.5, 
                                         selection.method = "vst", 
                                         nfeatures = 2000)


# Select integration features
SelectIntegrationFeatures(object.list = list(so_midface_E9.5, 
                                             so_midface_E10.5, 
                                             so_midface_E11.5),
                          nfeatures = 2000,
                          verbose = TRUE)

# Step 1: Get the variable features for both datasets
var_features_E9.5 <- so_midface_E9.5@assays[["SCT"]]@var.features
var_features_E10.5 <- so_midface_E10.5@assays[["SCT"]]@var.features
var_features_E11.5 <- so_midface_E11.5@assays[["SCT"]]@var.features

# Step 2: Find the common variable features between the two datasets
common_var_features <- Reduce(intersect, list(var_features_E9.5, 
                                              var_features_E10.5, 
                                              var_features_E11.5))

# Step 3: Prepare the objects for integration using the common features
objects <- list(so_midface_E9.5, 
                so_midface_E10.5, 
                so_midface_E11.5)

# Prepare objects for integration
objects <- PrepSCTIntegration(object.list = objects,
                              anchor.features = common_var_features,
                              verbose = TRUE)

# Step 4: Find integration anchors - make sure to specify the common features
anchors <- FindIntegrationAnchors(object.list = objects, 
                                  normalization.method = "SCT", 
                                  dims = 1:10, 
                                  anchor.features = common_var_features,  # explicitly specify the features
                                  k.anchor = 3,
                                  verbose = TRUE)

# Clear objects from workspace or clear all/close R and retrieve saved IntegrationAnchorSet
# [necessary due to Memory capacity error]
rm(so_midface_E9.5)
rm(so_midface_E10.5)
rm(so_midface_E11.5)
rm(objects)

# Step 5: Integrate the datasets using the found anchors
so_midface_E9.5_E10.5_E11.5_integrated <- IntegrateData(anchorset = anchors, 
                                                        normalization.method = "SCT",
                                                        dims = 1:10)

# Step 6: Perform scaling and PCA on the integrated data
so_midface_E9.5_E10.5_E11.5_integrated <- ScaleData(so_midface_E9.5_E10.5_E11.5_integrated)
so_midface_E9.5_E10.5_E11.5_integrated <- RunPCA(so_midface_E9.5_E10.5_E11.5_integrated, verbose = FALSE)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_midface_E9.5_E10.5_E11.5_integrated, ndims = 20)
out_path <- paste(output_folder, "integrated_midface_E9.5_E10.5_E11.5.data.qc.ellbowplot.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

# Store number of principle components in new variable (to be used later)
pca_dim_sel <- 15

# Perform UMAP on integrated data
so_midface_E9.5_E10.5_E11.5_integrated <- RunUMAP(so_midface_E9.5_E10.5_E11.5_integrated, 
                                                  dims = 1:pca_dim_sel, 
                                                  reduction = "pca", 
                                                  reduction.name = "umap.integrated")

# Change the default assay to "SCT" (normalized dataset)
DefaultAssay(so_midface_E9.5_E10.5_E11.5_integrated) <- "SCT"

# Clustering (Leiden) - Seurat v5 should work similarly
so_midface_E9.5_E10.5_E11.5_integrated <- FindNeighbors(so_midface_E9.5_E10.5_E11.5_integrated,
                                                        dims = 1:pca_dim_sel)
so_midface_E9.5_E10.5_E11.5_integrated <- FindClusters(so_midface_E9.5_E10.5_E11.5_integrated,
                                                       resolution = 0.5,
                                                       algorithm = 4,
                                                       graph.name = "integrated_snn")

# Set the desired order of orig.ident
so_midface_E9.5_E10.5_E11.5_integrated$orig.ident <- factor(
  so_midface_E9.5_E10.5_E11.5_integrated$orig.ident,
  levels = c("midface_E9.5", "midface_E10.5", "midface_E11.5")
)

# Visualize datasets as UMAP after Integration
p <- DimPlot(so_midface_E9.5_E10.5_E11.5_integrated, 
             reduction = "umap.integrated", 
             group.by = c("orig.ident", "seurat_clusters"))
outFile <- paste(output_folder, "/midface_E9.5_E10.5_E11.5_integrated.UMAP.orig.ident.clusters.pdf", sep = "")
pdf(outFile, width = 12, height = 5)
plot(p)
dev.off()

p <- DimPlot(object = so_midface_E9.5_E10.5_E11.5_integrated,
             reduction = "umap.integrated",
             group.by = 'seurat_clusters',
             split.by = 'orig.ident',
             label = TRUE)
outFile <- paste(output_folder, "/midface_E9.5_E10.5_E11.5_integrated.UMAP.orig.ident.clusters_split.pdf", sep = "")
pdf(outFile, width = 20, height = 5)
plot(p)
dev.off()


# Save the Seurat object
outFile <- paste(output_folder, "so_midface_E9.5_E10.5_E11.5_integrated.rds", sep = "")
saveRDS(so_midface_E9.5_E10.5_E11.5_integrated, file = outFile)


######## Visualize GOI as FeaturePlot

# Load pre-processed, integrated Seurat objects produced above
so_midface_E9.5_E10.5_E11.5_integrated <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_midfaceE9.5-11.5/analysis/so_midface_E9.5_E10.5_E11.5_integrated.rds")

# Set the desired order of orig.ident
so_midface_E9.5_E10.5_E11.5_integrated$orig.ident <- factor(
  so_midface_E9.5_E10.5_E11.5_integrated$orig.ident,
  levels = c("midface_E9.5", "midface_E10.5", "midface_E11.5")
)

# Definition of GOI (TFs and Senesence pathway genes)
goi <- c("Pbx1", "Pbx2", "Pbx3", "Pbx4", "Zfhx3", "Zfhx4",           
         "Grhl3", "Snai2", "Tfap2a", "Tfap2b", "Twist1", "Barx1", "Cxxc4", 
         "Aldh1a2", "Dlx1", "Dlx2", "Sox9",                                     # Marker Genes for Epithelial Seam Cells at the Lambdoidal Junction
         "Tgfb1", "Smad2", "Smad4", "Smad5",                                    # TGF-beta pathway
         "Mtor", "Akt1",                                                        # mTOR pathway
         "Foxo1", "Foxo3", "Foxo4", "Foxo6",                                    # FOXO TFs
         "Cdkn2a")                                                              # p19ARF

# Change the default assay to "SCT"
DefaultAssay(so_midface_E9.5_E10.5_E11.5_integrated) <- "SCT"

# Plots for goi
outFile <- paste(output_folder,
                 "/midface_E9.5_E10.5_E11.5_integrated.UMAP.goi.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 20, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi) {
  if (gene %in% rownames(so_midface_E9.5_E10.5_E11.5_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_midface_E9.5_E10.5_E11.5_integrated, 
                     features = gene,
                     reduction = "umap.integrated",
                     split.by = "orig.ident")
    plot(p)
  } else {
    # Print a message for missing genes (optional)
    message(paste("Gene not found in data: ", gene))
  }
}

dev.off()

# Definition of GOI (TFs and Senesence pathway genes)
goi_nicotine <- c("Pbx1", "Pbx2", "Pbx3", "Pbx4", "Zfhx3", "Zfhx4",           
                  "Chrna3", "Chrna5", "Chrnb4",                                 # nAchRs
                  "Cyp2a5", "Ankk1",                                            # nicotine pathway and metabolism
                  "Slc6a3", "Slc18a2", "Th", "Ddc", "Drd2")                     # dopamine release

# Change the default assay to "SCT"
DefaultAssay(so_midface_E9.5_E10.5_E11.5_integrated) <- "SCT"

# Plots for goi
outFile <- paste(output_folder,
                 "/midface_E9.5_E10.5_E11.5_integrated.UMAP.goi_nicotine.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 20, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi_nicotine) {
  if (gene %in% rownames(so_midface_E9.5_E10.5_E11.5_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_midface_E9.5_E10.5_E11.5_integrated, 
                     features = gene,
                     reduction = "umap.integrated",
                     split.by = "orig.ident")
    plot(p)
  } else {
    # Print a message for missing genes (optional)
    message(paste("Gene not found in data: ", gene))
  }
}

dev.off()


##### UMAP Coexpression of Pbx1 and Zfhx3
# Define a list of gene sets (co-expressed genes)
gene_set_1 <- list(Coexpression = c("Pbx1",
                                    "Zfhx3"))

# Add module scores to the Seurat object
data <- AddModuleScore(so_midface_E9.5_E10.5_E11.5_integrated, 
                       features = gene_set_1, name = "CoexpressionScore")

# Visualize the module score in UMAP
p <- FeaturePlot(data, features = "CoexpressionScore1",
                 pt.size = 0.5) +
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = median(data$CoexpressionScore1)) +
  labs(title = "Co-expression of Pbx1 and Zfhx3 in E9.5-11.5 midface scRNA-seq", 
       x = "UMAP 1", 
       y = "UMAP 2", 
       color = "Co-expression\n(Pbx1, Zfhx3)") +  # Custom legend title
  theme_minimal()
outFile <- paste(output_folder, "/midface_E9.5_E10.5_E11.5_integrated.UMAP.coexpressionPbx1+Zfhx3.pdf", sep = "")
pdf(outFile, width = 10, height = 5)
plot(p)
dev.off()


### Determine marker genes for each cluster
so_midface_E9.5_E10.5_E11.5_integrated <- PrepSCTFindMarkers(so_midface_E9.5_E10.5_E11.5_integrated)

markers <- FindAllMarkers(
  object = so_midface_E9.5_E10.5_E11.5_integrated,
  only.pos = TRUE,               # Return only positive markers (upregulated in the cluster)
  min.pct = 0.25,                # Gene expressed in at least 25% of cells in either group
  logfc.threshold = 0.25         # Minimum log fold change
)

# Save the marker list to a CSV file
write.csv(markers, file = "/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_midfaceE9.5-11.5/analysis/midface_E9.5_E10.5_E11.5_integrated.markers.csv")

