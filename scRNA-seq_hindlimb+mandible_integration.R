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
saveRDS(so_hindlimb_E10.5, file = "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE10.5/so_hindlimb_E10.5.rds")

# Load data
so_hindlimb_E10.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE10.5/so_hindlimb_E10.5.rds")

FeaturePlot(so_hindlimb_E10.5, 
            features = "Zfhx3",
            #reduction = "umap.integrated",
            #split.by = "orig.ident",
            pt.size = 0.2)

# UMAP overlay of goi
goi <- c("Pbx1", "Pbx2", "Pbx3", "Hand2", "Irx3","Irx5", "Meis1", "Meis2",         # TALE-HD and Hand2
         "Hoxa1", "Hoxa2", "Hoxa3", "Hoxa4", "Hoxa5", "Hoxa6", "Hoxa7", "Hoxa9",   # Hoxa cluster
         "Hoxa10", "Hoxa11", "Hoxa13", "Hoxd13")
out_file <- paste(out_folder, "/hindlimb_E10.5.UMAP.goi.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi) {
  p <- FeaturePlot(so_hindlimb_E10.5, gene)
  plot(p)
}
dev.off()

# Multipannel UMAP of goi
p <- FeaturePlot(so_hindlimb_E10.5, features = goi,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(out_folder, "/hindlimb_E10.5.UMAP.goi.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()


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

FeaturePlot(so_hindlimb_E11.5, features = goi,
            cols = c('lightgray', 'blue'),
            pt.size = 0.01)   # Adjust pt.size to your desired value


################################################################################
################## Pre-analysis of autopod_E12.5 dataset  #####################

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE12.5_autopod/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
# Downloaded from ArrayExpress: E-MTAB-10514, https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10514#processed-data
autopod_E12.5 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_mouse_hindlimb_publicDomain/WSSS_THYst9807813_hindlimb-autopod_E12.5/filtered_feature_bc_matrix')

# Create Seurat object
so_autopod_E12.5 <- CreateSeuratObject(counts = autopod_E12.5, 
                                       project = "autopod_E12.5", 
                                       min.cells = 3, min.features = 1000)  # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_autopod_E12.5 <- PercentageFeatureSet(so_autopod_E12.5, 
                                         pattern = "^mt-", 
                                         col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_autopod_E12.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "autopod_E12.5.qc.unfiltered.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

### Selecting cells for further analysis
percent.mt_max <- 5  # maximum percentage of mitochondrial genes (adjust as needed)
nFeature_RNA_min <- 100  # minimum number of features per cell
nFeature_RNA_max <- 10000 # minimum number of features per cell
nCount_RNA_min <- 100  # minimum number of RNA counts per cell
nCount_RNA_max <- 30000  # maximum number of RNA counts per cell

# Subset the Seurat object (filter based on thresholds above)
so_autopod_E12.5 <- subset(so_autopod_E12.5,
                           subset = percent.mt <= percent.mt_max &
                             nFeature_RNA >= nFeature_RNA_min &
                             nFeature_RNA <= nFeature_RNA_max &
                             nCount_RNA >= nCount_RNA_min &
                             nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_autopod_E12.5,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "autopod_E12.5.data.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_autopod_E12.5[["RNA"]]), invert = TRUE)
so_autopod_E12.5 <- CreateSeuratObject(counts = GetAssayData(so_autopod_E12.5[["RNA"]], layer = "counts")[keep, ], 
                                       project = "autopod_E12.5", 
                                       meta.data = so_autopod_E12.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_autopod_E12.5 <- SCTransform(so_autopod_E12.5,
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
so_autopod_E12.5 <- CellCycleScoring(so_autopod_E12.5,
                                     s.features = cc.genes_corrected$s.genes,   # S phase genes
                                     g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                     set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_autopod_E12.5 <- ScaleData(so_autopod_E12.5,
                              features = rownames(so_autopod_E12.5),
                              vars.to.regress = c("S.Score", "G2M.Score"),
                              verbose = TRUE
)

# PCA
so_autopod_E12.5 <- RunPCA(so_autopod_E12.5,
                           verbose = FALSE, 
                           npcs = 100,                # Number of principal components to compute
                           ndims.print = 1:5,         # Print details for the first 5 PCs
                           nfeatures.print = 30       # Print details for the top 30 features
)


# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_autopod_E12.5, ndims = 20)
out_file <- paste(output_folder, "autopod_E12.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 10                                      

#UMAP
so_autopod_E12.5 <- RunUMAP(so_autopod_E12.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_autopod_E12.5 <- FindNeighbors(so_autopod_E12.5, 
                                  dims = 1:pca_dim_sel)     # "Cell clusters were identified by constructing a shared nearest neighbor graph followed by a modularity optimization-based clustering algorithm (Leiden algorithm) using the top 60 principal components as determined by PCA." (Losa et al., 2023)
so_autopod_E12.5 <- FindClusters(so_autopod_E12.5, 
                                 resolution = 0.8,          # "Clustering was performed at multiple resolutions between 0.2 and 2, and optimal resolution was determined empirically based on the expression of known population markers (resolution = 0.8)." (Losa et al., 2023)
                                 algorithm = 4)             # Leiden clustering

# Visualize clusters as Dimplot
p <- DimPlot(object = so_autopod_E12.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/autopod_E12.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_autopod_E12.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/autopod_E12.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_autopod_E12.5, file = "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE12.5_autopod/so_autopod_E12.5.rds")

# Load data
so_autopod_E12.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE12.5_autopod/so_autopod_E12.5.rds")

FeaturePlot(so_autopod_E12.5, 
            features = "Pbx1",
            #reduction = "umap.integrated",
            #split.by = "orig.ident",
            pt.size = 0.2)

################################################################################
################## Pre-analysis of stylopod_E12.5 dataset  #####################

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE12.5_stylopod/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
# Downloaded from ArrayExpress: E-MTAB-10514, https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10514#processed-data
stylopod_E12.5 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_mouse_hindlimb_publicDomain/WSSS_THYst9807811_hindlimb_stylopod_E12.5/filtered_feature_bc_matrix')

# Create Seurat object
so_stylopod_E12.5 <- CreateSeuratObject(counts = stylopod_E12.5, 
                                        project = "stylopod_E12.5", 
                                        min.cells = 3, min.features = 1000)  # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_stylopod_E12.5 <- PercentageFeatureSet(so_stylopod_E12.5, 
                                          pattern = "^mt-", 
                                          col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_stylopod_E12.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "stylopod_E12.5.qc.unfiltered.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

### Selecting cells for further analysis
percent.mt_max <- 5  # maximum percentage of mitochondrial genes (adjust as needed)
nFeature_RNA_min <- 100  # minimum number of features per cell
nFeature_RNA_max <- 10000 # minimum number of features per cell
nCount_RNA_min <- 100  # minimum number of RNA counts per cell
nCount_RNA_max <- 30000  # maximum number of RNA counts per cell

# Subset the Seurat object (filter based on thresholds above)
so_stylopod_E12.5 <- subset(so_stylopod_E12.5,
                            subset = percent.mt <= percent.mt_max &
                              nFeature_RNA >= nFeature_RNA_min &
                              nFeature_RNA <= nFeature_RNA_max &
                              nCount_RNA >= nCount_RNA_min &
                              nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_stylopod_E12.5,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "stylopod_E12.5.data.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_stylopod_E12.5[["RNA"]]), invert = TRUE)
so_stylopod_E12.5 <- CreateSeuratObject(counts = GetAssayData(so_stylopod_E12.5[["RNA"]], layer = "counts")[keep, ], 
                                        project = "stylopod_E12.5", 
                                        meta.data = so_stylopod_E12.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_stylopod_E12.5 <- SCTransform(so_stylopod_E12.5,
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
so_stylopod_E12.5 <- CellCycleScoring(so_stylopod_E12.5,
                                      s.features = cc.genes_corrected$s.genes,   # S phase genes
                                      g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                      set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_stylopod_E12.5 <- ScaleData(so_stylopod_E12.5,
                               features = rownames(so_stylopod_E12.5),
                               vars.to.regress = c("S.Score", "G2M.Score"),
                               verbose = TRUE
)

# PCA
so_stylopod_E12.5 <- RunPCA(so_stylopod_E12.5,
                            verbose = FALSE, 
                            npcs = 100,                # Number of principal components to compute
                            ndims.print = 1:5,         # Print details for the first 5 PCs
                            nfeatures.print = 30       # Print details for the top 30 features
)


# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_stylopod_E12.5, ndims = 20)
out_file <- paste(output_folder, "stylopod_E12.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 12                                      

#UMAP
so_stylopod_E12.5 <- RunUMAP(so_stylopod_E12.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_stylopod_E12.5 <- FindNeighbors(so_stylopod_E12.5, 
                                   dims = 1:pca_dim_sel)     # "Cell clusters were identified by constructing a shared nearest neighbor graph followed by a modularity optimization-based clustering algorithm (Leiden algorithm) using the top 60 principal components as determined by PCA." (Losa et al., 2023)
so_stylopod_E12.5 <- FindClusters(so_stylopod_E12.5, 
                                  resolution = 0.8,          # "Clustering was performed at multiple resolutions between 0.2 and 2, and optimal resolution was determined empirically based on the expression of known population markers (resolution = 0.8)." (Losa et al., 2023)
                                  algorithm = 4)             # Leiden clustering

# Visualize clusters as Dimplot
p <- DimPlot(object = so_stylopod_E12.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/stylopod_E12.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_stylopod_E12.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/stylopod_E12.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_stylopod_E12.5, file = "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE12.5_stylopod/so_stylopod_E12.5.rds")

# Load data
so_stylopod_E12.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE12.5_stylopod/so_stylopod_E12.5.rds")

FeaturePlot(so_stylopod_E12.5, 
            features = "Pbx1",
            #reduction = "umap.integrated",
            #split.by = "orig.ident",
            pt.size = 0.2)

################################################################################
################## Pre-analysis of stylopod_zeugopod_E12.5 dataset  #####################

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE12.5_stylopod_zeugopod/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
# Downloaded from ArrayExpress: E-MTAB-10514, https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10514#processed-data
stylopod_zeugopod_E12.5 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_mouse_hindlimb_publicDomain/WSSS_THYst9807812_hindlimb_stylopod_zeugopod_E12.5/filtered_feature_bc_matrix')

# Create Seurat object
so_stylopod_zeugopod_E12.5 <- CreateSeuratObject(counts = stylopod_zeugopod_E12.5, 
                                                 project = "stylopod_zeugopod_E12.5", 
                                                 min.cells = 3, min.features = 1000)  # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_stylopod_zeugopod_E12.5 <- PercentageFeatureSet(so_stylopod_zeugopod_E12.5, 
                                                   pattern = "^mt-", 
                                                   col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_stylopod_zeugopod_E12.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "stylopod_zeugopod_E12.5.qc.unfiltered.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

### Selecting cells for further analysis
percent.mt_max <- 5  # maximum percentage of mitochondrial genes (adjust as needed)
nFeature_RNA_min <- 100  # minimum number of features per cell
nFeature_RNA_max <- 10000 # minimum number of features per cell
nCount_RNA_min <- 100  # minimum number of RNA counts per cell
nCount_RNA_max <- 30000  # maximum number of RNA counts per cell

# Subset the Seurat object (filter based on thresholds above)
so_stylopod_zeugopod_E12.5 <- subset(so_stylopod_zeugopod_E12.5,
                                     subset = percent.mt <= percent.mt_max &
                                       nFeature_RNA >= nFeature_RNA_min &
                                       nFeature_RNA <= nFeature_RNA_max &
                                       nCount_RNA >= nCount_RNA_min &
                                       nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_stylopod_zeugopod_E12.5,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "stylopod_zeugopod_E12.5.data.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_stylopod_zeugopod_E12.5[["RNA"]]), invert = TRUE)
so_stylopod_zeugopod_E12.5 <- CreateSeuratObject(counts = GetAssayData(so_stylopod_zeugopod_E12.5[["RNA"]], layer = "counts")[keep, ], 
                                                 project = "stylopod_zeugopod_E12.5", 
                                                 meta.data = so_stylopod_zeugopod_E12.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_stylopod_zeugopod_E12.5 <- SCTransform(so_stylopod_zeugopod_E12.5,
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
so_stylopod_zeugopod_E12.5 <- CellCycleScoring(so_stylopod_zeugopod_E12.5,
                                               s.features = cc.genes_corrected$s.genes,   # S phase genes
                                               g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                               set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_stylopod_zeugopod_E12.5 <- ScaleData(so_stylopod_zeugopod_E12.5,
                                        features = rownames(so_stylopod_zeugopod_E12.5),
                                        vars.to.regress = c("S.Score", "G2M.Score"),
                                        verbose = TRUE
)

# PCA
so_stylopod_zeugopod_E12.5 <- RunPCA(so_stylopod_zeugopod_E12.5,
                                     verbose = FALSE, 
                                     npcs = 100,                # Number of principal components to compute
                                     ndims.print = 1:5,         # Print details for the first 5 PCs
                                     nfeatures.print = 30       # Print details for the top 30 features
)


# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_stylopod_zeugopod_E12.5, ndims = 20)
out_file <- paste(output_folder, "stylopod_zeugopod_E12.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 9                                      

#UMAP
so_stylopod_zeugopod_E12.5 <- RunUMAP(so_stylopod_zeugopod_E12.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_stylopod_zeugopod_E12.5 <- FindNeighbors(so_stylopod_zeugopod_E12.5, 
                                            dims = 1:pca_dim_sel)     # "Cell clusters were identified by constructing a shared nearest neighbor graph followed by a modularity optimization-based clustering algorithm (Leiden algorithm) using the top 60 principal components as determined by PCA." (Losa et al., 2023)
so_stylopod_zeugopod_E12.5 <- FindClusters(so_stylopod_zeugopod_E12.5, 
                                           resolution = 0.8,          # "Clustering was performed at multiple resolutions between 0.2 and 2, and optimal resolution was determined empirically based on the expression of known population markers (resolution = 0.8)." (Losa et al., 2023)
                                           algorithm = 4)             # Leiden clustering

# Visualize clusters as Dimplot
p <- DimPlot(object = so_stylopod_zeugopod_E12.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/stylopod_zeugopod_E12.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_stylopod_zeugopod_E12.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/stylopod_zeugopod_E12.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_stylopod_zeugopod_E12.5, file = "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE12.5_stylopod_zeugopod/so_stylopod_zeugopod_E12.5.rds")

# Load data
so_stylopod_zeugopod_E12.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE12.5_stylopod_zeugopod/so_stylopod_zeugopod_E12.5.rds")

FeaturePlot(so_stylopod_zeugopod_E12.5, 
            features = "Pbx1",
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
# Downloaded from GEO: GSM4227225
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

FeaturePlot(so_hindlimb_E13.5, features = goi,
            cols = c('lightgray', 'blue'),
            pt.size = 0.01)   # Adjust pt.size to your desired value



################################################################################
################## Pre-analysis of autopod_E13.5 dataset  #####################

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE13.5_autopod/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
# Downloaded from ArrayExpress: E-MTAB-10514, https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10514#processed-data
autopod_E13.5 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_mouse_hindlimb_publicDomain/WSSS_THYst9807819_hindlimb-autopod_E13.5/filtered_feature_bc_matrix')

# Create Seurat object
so_autopod_E13.5 <- CreateSeuratObject(counts = autopod_E13.5, 
                                       project = "autopod_E13.5", 
                                       min.cells = 3, min.features = 1000)  # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_autopod_E13.5 <- PercentageFeatureSet(so_autopod_E13.5, 
                                         pattern = "^mt-", 
                                         col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_autopod_E13.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "autopod_E13.5.qc.unfiltered.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

### Selecting cells for further analysis
percent.mt_max <- 5  # maximum percentage of mitochondrial genes (adjust as needed)
nFeature_RNA_min <- 100  # minimum number of features per cell
nFeature_RNA_max <- 10000 # minimum number of features per cell
nCount_RNA_min <- 100  # minimum number of RNA counts per cell
nCount_RNA_max <- 30000  # maximum number of RNA counts per cell

# Subset the Seurat object (filter based on thresholds above)
so_autopod_E13.5 <- subset(so_autopod_E13.5,
                           subset = percent.mt <= percent.mt_max &
                             nFeature_RNA >= nFeature_RNA_min &
                             nFeature_RNA <= nFeature_RNA_max &
                             nCount_RNA >= nCount_RNA_min &
                             nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_autopod_E13.5,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "autopod_E13.5.data.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_autopod_E13.5[["RNA"]]), invert = TRUE)
so_autopod_E13.5 <- CreateSeuratObject(counts = GetAssayData(so_autopod_E13.5[["RNA"]], layer = "counts")[keep, ], 
                                       project = "autopod_E13.5", 
                                       meta.data = so_autopod_E13.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_autopod_E13.5 <- SCTransform(so_autopod_E13.5,
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
so_autopod_E13.5 <- CellCycleScoring(so_autopod_E13.5,
                                     s.features = cc.genes_corrected$s.genes,   # S phase genes
                                     g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                     set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_autopod_E13.5 <- ScaleData(so_autopod_E13.5,
                              features = rownames(so_autopod_E13.5),
                              vars.to.regress = c("S.Score", "G2M.Score"),
                              verbose = TRUE
)

# PCA
so_autopod_E13.5 <- RunPCA(so_autopod_E13.5,
                           verbose = FALSE, 
                           npcs = 100,                # Number of principal components to compute
                           ndims.print = 1:5,         # Print details for the first 5 PCs
                           nfeatures.print = 30       # Print details for the top 30 features
)


# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_autopod_E13.5, ndims = 20)
out_file <- paste(output_folder, "autopod_E13.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 13                                      

#UMAP
so_autopod_E13.5 <- RunUMAP(so_autopod_E13.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_autopod_E13.5 <- FindNeighbors(so_autopod_E13.5, 
                                  dims = 1:pca_dim_sel)     # "Cell clusters were identified by constructing a shared nearest neighbor graph followed by a modularity optimization-based clustering algorithm (Leiden algorithm) using the top 60 principal components as determined by PCA." (Losa et al., 2023)
so_autopod_E13.5 <- FindClusters(so_autopod_E13.5, 
                                 resolution = 0.8,          # "Clustering was performed at multiple resolutions between 0.2 and 2, and optimal resolution was determined empirically based on the expression of known population markers (resolution = 0.8)." (Losa et al., 2023)
                                 algorithm = 4)             # Leiden clustering

# Visualize clusters as Dimplot
p <- DimPlot(object = so_autopod_E13.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/autopod_E13.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_autopod_E13.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/autopod_E13.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_autopod_E13.5, file = "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE13.5_autopod/so_autopod_E13.5.rds")

# Load data
so_autopod_E13.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE13.5_autopod/so_autopod_E13.5.rds")

FeaturePlot(so_autopod_E13.5, 
            features = "Hoxa13",
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
# Downloaded from GEO: GSM4227226
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
            features = "Fgf10",
            #reduction = "umap.integrated",
            #split.by = "orig.ident",
            pt.size = 0.2)

FeaturePlot(so_hindlimb_E15.5, features = goi,
            cols = c('lightgray', 'blue'),
            pt.size = 0.01)   # Adjust pt.size to your desired value

################################################################################
################## Pre-analysis of hindlimb_E18.5 dataset  #####################

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE18.5/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
# Downloaded from GEO: GSM4227227
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


################################################################################
################### Pre-analysis of mandible_E9.5 dataset  ######################

# Set out folder (to store results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE9.5/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
mandible_E9.5 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/mandible_E9.5_1-YSGM_FB00001195')

# Create Seurat object
so_mandible_E9.5 <- CreateSeuratObject(counts = mandible_E9.5, 
                                       project = "mandible_E9.5", 
                                       min.cells = 3, min.features = 1000)      # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_mandible_E9.5 <- PercentageFeatureSet(so_mandible_E9.5, 
                                         pattern = "^mt-", 
                                         col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_mandible_E9.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "mandible_E9.5.qc.unfiltered.pdf", sep = "")
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
so_mandible_E9.5 <- subset(so_mandible_E9.5,
                           subset = percent.mt <= percent.mt_max &
                             nFeature_RNA >= nFeature_RNA_min &
                             nFeature_RNA <= nFeature_RNA_max &
                             nCount_RNA >= nCount_RNA_min &
                             nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_mandible_E9.5,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "mandible_E9.5.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_mandible_E9.5[["RNA"]]), invert = TRUE)
so_mandible_E9.5 <- CreateSeuratObject(counts = GetAssayData(so_mandible_E9.5[["RNA"]], layer = "counts")[keep, ], 
                                       project = "mandible_E9.5", 
                                       meta.data = so_mandible_E9.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_mandible_E9.5 <- SCTransform(so_mandible_E9.5,
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
so_mandible_E9.5 <- CellCycleScoring(so_mandible_E9.5,
                                     s.features = cc.genes_corrected$s.genes,   # S phase genes
                                     g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                     set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_mandible_E9.5 <- ScaleData(so_mandible_E9.5,
                              features = rownames(so_mandible_E9.5),
                              vars.to.regress = c("S.Score", "G2M.Score"),
                              verbose = TRUE
)

# PCA
so_mandible_E9.5 <- RunPCA(so_mandible_E9.5,
                           verbose = FALSE, 
                           npcs = 100,                # Number of principal components to compute
                           ndims.print = 1:5,         # Print details for the first 5 PCs
                           nfeatures.print = 30       # Print details for the top 30 features
)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_mandible_E9.5, ndims = 20)
out_file <- paste(output_folder, "mandible_E9.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 11   

#UMAP
so_mandible_E9.5 <- RunUMAP(so_mandible_E9.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_mandible_E9.5 <- FindNeighbors(so_mandible_E9.5, dims = 1:pca_dim_sel)         
so_mandible_E9.5 <- FindClusters(so_mandible_E9.5, 
                                 resolution = 0.5,
                                 algorithm = 4)            

# Visualize clusters as Dimplot
p <- DimPlot(object = so_mandible_E9.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/mandible_E9.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_mandible_E9.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/mandible_E9.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_mandible_E9.5, file = "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE9.5/so_mandible_E9.5.rds")

# Load data
so_mandible_E9.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE9.5/so_mandible_E9.5.rds")

FeaturePlot(so_mandible_E9.5, 
            features = "Hand2",
            #reduction = "umap.integrated",
            #split.by = "orig.ident",
            pt.size = 0.2)

################################################################################
################## Pre-analysis of mandible_E10.5 dataset  #####################

# Set out folder (to store results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE10.5/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
mandible_E10.5 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/mandible_E10.5_P-R4YP_FB00001195')

# Create Seurat object
so_mandible_E10.5 <- CreateSeuratObject(counts = mandible_E10.5, 
                                        project = "mandible_E10.5", 
                                        min.cells = 3, min.features = 1000)      # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_mandible_E10.5 <- PercentageFeatureSet(so_mandible_E10.5, 
                                          pattern = "^mt-", 
                                          col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_mandible_E10.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "mandible_E10.5.qc.unfiltered.pdf", sep = "")
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
so_mandible_E10.5 <- subset(so_mandible_E10.5,
                            subset = percent.mt <= percent.mt_max &
                              nFeature_RNA >= nFeature_RNA_min &
                              nFeature_RNA <= nFeature_RNA_max &
                              nCount_RNA >= nCount_RNA_min &
                              nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_mandible_E10.5,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "mandible_E10.5.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_mandible_E10.5[["RNA"]]), invert = TRUE)
so_mandible_E10.5 <- CreateSeuratObject(counts = GetAssayData(so_mandible_E10.5[["RNA"]], layer = "counts")[keep, ], 
                                        project = "mandible_E10.5", 
                                        meta.data = so_mandible_E10.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_mandible_E10.5 <- SCTransform(so_mandible_E10.5,
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
so_mandible_E10.5 <- CellCycleScoring(so_mandible_E10.5,
                                      s.features = cc.genes_corrected$s.genes,   # S phase genes
                                      g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                      set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_mandible_E10.5 <- ScaleData(so_mandible_E10.5,
                               features = rownames(so_mandible_E10.5),
                               vars.to.regress = c("S.Score", "G2M.Score"),
                               verbose = TRUE
)

# PCA
so_mandible_E10.5 <- RunPCA(so_mandible_E10.5,
                            verbose = FALSE, 
                            npcs = 100,                # Number of principal components to compute
                            ndims.print = 1:5,         # Print details for the first 5 PCs
                            nfeatures.print = 30       # Print details for the top 30 features
)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_mandible_E10.5, ndims = 20)
out_file <- paste(output_folder, "mandible_E10.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 9   

#UMAP
so_mandible_E10.5 <- RunUMAP(so_mandible_E10.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_mandible_E10.5 <- FindNeighbors(so_mandible_E10.5, dims = 1:pca_dim_sel)         
so_mandible_E10.5 <- FindClusters(so_mandible_E10.5, 
                                  resolution = 0.5,
                                  algorithm = 4)            

# Visualize clusters as Dimplot
p <- DimPlot(object = so_mandible_E10.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/mandible_E10.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_mandible_E10.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/mandible_E10.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_mandible_E10.5, file = "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE10.5/so_mandible_E10.5.rds")

# Load data
so_mandible_E10.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE10.5/so_mandible_E10.5.rds")

FeaturePlot(so_mandible_E10.5, 
            features = "Hand2",
            #reduction = "umap.integrated",
            #split.by = "orig.ident",
            pt.size = 0.2)

################################################################################
################### Pre-analysis of mandible_E11.5 dataset  ######################

# Set out folder (to store results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE11.5/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
mandible_E11.5 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/mandible_E11.5_P-R4YR_FB00001195')

# Create Seurat object
so_mandible_E11.5 <- CreateSeuratObject(counts = mandible_E11.5, 
                                        project = "mandible_E11.5", 
                                        min.cells = 3, min.features = 1000)      # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_mandible_E11.5 <- PercentageFeatureSet(so_mandible_E11.5, 
                                          pattern = "^mt-", 
                                          col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_mandible_E11.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "mandible_E11.5.qc.unfiltered.pdf", sep = "")
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
so_mandible_E11.5 <- subset(so_mandible_E11.5,
                            subset = percent.mt <= percent.mt_max &
                              nFeature_RNA >= nFeature_RNA_min &
                              nFeature_RNA <= nFeature_RNA_max &
                              nCount_RNA >= nCount_RNA_min &
                              nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_mandible_E11.5,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "mandible_E11.5.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_mandible_E11.5[["RNA"]]), invert = TRUE)
so_mandible_E11.5 <- CreateSeuratObject(counts = GetAssayData(so_mandible_E11.5[["RNA"]], layer = "counts")[keep, ], 
                                        project = "mandible_E11.5", 
                                        meta.data = so_mandible_E11.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 2000
so_mandible_E11.5 <- SCTransform(so_mandible_E11.5,
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
so_mandible_E11.5 <- CellCycleScoring(so_mandible_E11.5,
                                      s.features = cc.genes_corrected$s.genes,   # S phase genes
                                      g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                      set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_mandible_E11.5 <- ScaleData(so_mandible_E11.5,
                               features = rownames(so_mandible_E11.5),
                               vars.to.regress = c("S.Score", "G2M.Score"),
                               verbose = TRUE
)

# PCA
so_mandible_E11.5 <- RunPCA(so_mandible_E11.5,
                            verbose = FALSE, 
                            npcs = 100,                # Number of principal components to compute
                            ndims.print = 1:5,         # Print details for the first 5 PCs
                            nfeatures.print = 30       # Print details for the top 30 features
)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_mandible_E11.5, ndims = 20)
out_file <- paste(output_folder, "mandible_E11.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 8   

#UMAP
so_mandible_E11.5 <- RunUMAP(so_mandible_E11.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_mandible_E11.5 <- FindNeighbors(so_mandible_E11.5, dims = 1:pca_dim_sel)         
so_mandible_E11.5 <- FindClusters(so_mandible_E11.5, 
                                  resolution = 0.5,
                                  algorithm = 4)            

# Visualize clusters as Dimplot
p <- DimPlot(object = so_mandible_E11.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/mandible_E11.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_mandible_E11.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/mandible_E11.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_mandible_E11.5, file = "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE11.5/so_mandible_E11.5.rds")

# Load data
so_mandible_E11.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE11.5/so_mandible_E11.5.rds")

FeaturePlot(so_mandible_E11.5, 
            features = "Hand2",
            #reduction = "umap.integrated",
            #split.by = "orig.ident",
            pt.size = 0.2)


################################################################################
############################ DEFINITION OF CLUSTERS ############################

### GOI for TALE-HD/Hand2/Hox cluster, cell type and spatial distribution
# TALE-HD, Hand2 and Hox expression
goi <- c("Pbx1", "Pbx2", "Pbx3", "Hand2", "Irx3","Irx5", "Meis1", "Meis2",      # TALE-HD and Hand2
         "Hoxa1", "Hoxa2", "Hoxa3", "Hoxa4", "Hoxa5", "Hoxa6", "Hoxa7",         # Hox cluster genes
         "Hoxa9", "Hoxa10", "Hoxa11", "Hoxa13", "Hoxd13")

# Celltype-clusters via marker genes
goi_celltype <- c("Lin28a", "Sall4", "Fzd7",                                    # undifferentiated cells; Fernandez-Guerrero et al., 2021
                  "Bmp2", "Col2a1", "Sox9",                                     # chondrogenic differentiation and limb skeletal, digit and joint morphogenesis; Fernandez-Guerrero et al., 2021
                  "Runx2", "Dlx5", "Sp7",                                       # osteoblast progenitor and differentiation marker 
                  "Alx3", "Alx4",                                               # proximal anterior mesoderm; Fernandez-Guerrero et al., 2021                                             # autopod-associated genes; Fernandez-Guerrero et al., 2021
                  "Asb4",                                                       # distal anterior mesoderm; Fernandez-Guerrero et al., 2021 
                  "Krt8", "Krt14", "Krt15",                                     # epidermis/ epidermal keratins; Kelly et al., 2020; Fernandez-Guerrero et al., 2021 
                  "Cdh5",                                                       # vasculature @E11.5; Kelly et al., 2020
                  "Lyz2",                                                       # blood @E11.5; Kelly et al., 2020
                  "Bcan", "Dsg2", "Esrp1"                                       # epithelial markers; Fernandez-Guerrero et al., 2021 
)

# Spatial distribution via marker genes
goi_spatial <- c("Hoxa9", "Hoxd9", "Shox2",                                     # stylopod-associated genes; Fernandez-Guerrero et al., 2021
                 "Hoxa11",                                                      # zeugopod-associated genes; Fernandez-Guerrero et al., 2021
                 "Sall1", "Sall3", "Hoxa13",                                    # autopod-associated genes; Fernandez-Guerrero et al., 2021
                 "Hoxd11", "Hoxd12", "Hoxd13",                                  # digit patterning; Fernandez-Guerrero et al., 2021
                 "Fgf8", "Wwc1", "Pdgfa"                                        # AER; Fernandez-Guerrero et al., 2021
)


######################## Plots for pre-analysis datasets #######################

##### hindlimb E10.5
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE10.5/"

# UMAP of goi
out_file <- paste(output_folder, "/hindlimb_E10.5.UMAP.goi.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi) {
  p <- FeaturePlot(so_hindlimb_E10.5, gene)
  plot(p)
}
dev.off()

# Multipannel UMAP of goi
p <- FeaturePlot(so_hindlimb_E10.5, features = goi,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/hindlimb_E10.5.UMAP.goi.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()

# UMAP of goi_celltype
out_file <- paste(output_folder, "/hindlimb_E10.5.UMAP.goi_celltype.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi_celltype) {
  p <- FeaturePlot(so_hindlimb_E10.5, gene)
  plot(p)
}
dev.off()

# Multipannel UMAP of goi_celltype
p <- FeaturePlot(so_hindlimb_E10.5, features = goi_celltype,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/hindlimb_E10.5.UMAP.goi_celltype.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()

# UMAP of goi_spatial
out_file <- paste(output_folder, "/hindlimb_E10.5.UMAP.goi_spatial.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi_spatial) {
  p <- FeaturePlot(so_hindlimb_E10.5, gene)
  plot(p)
}
dev.off()

# Multipannel UMAP of goi_spatial
p <- FeaturePlot(so_hindlimb_E10.5, features = goi_spatial,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/hindlimb_E10.5.UMAP.goi_spatial.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()


##### E11.5
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE11.5/"

# UMAP of goi
out_file <- paste(output_folder, "/hindlimb_E11.5.UMAP.goi.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi) {
  if (gene %in% rownames(so_hindlimb_E11.5[["SCT"]])) {
    p <- FeaturePlot(so_hindlimb_E11.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi
p <- FeaturePlot(so_hindlimb_E11.5, features = goi,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/hindlimb_E11.5.UMAP.goi.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()

# UMAP of goi_celltype
out_file <- paste(output_folder, "/hindlimb_E11.5.UMAP.goi_celltype.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi_celltype) {
  if (gene %in% rownames(so_hindlimb_E11.5[["SCT"]])) {
    p <- FeaturePlot(so_hindlimb_E11.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi_celltype
p <- FeaturePlot(so_hindlimb_E11.5, features = goi_celltype,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/hindlimb_E11.5.UMAP.goi_celltype.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()

# UMAP of goi_spatial
out_file <- paste(output_folder, "/hindlimb_E11.5.UMAP.goi_spatial.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi_spatial) {
  if (gene %in% rownames(so_hindlimb_E11.5[["SCT"]])) {
    p <- FeaturePlot(so_hindlimb_E11.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi_spatial
p <- FeaturePlot(so_hindlimb_E11.5, features = goi_spatial,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/hindlimb_E11.5.UMAP.goi_spatial.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()


##### hindlimb E13.5
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE13.5/"

# UMAP of goi
out_file <- paste(output_folder, "/hindlimb_E13.5.UMAP.goi.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi) {
  if (gene %in% rownames(so_hindlimb_E13.5[["SCT"]])) {
    p <- FeaturePlot(so_hindlimb_E13.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi
p <- FeaturePlot(so_hindlimb_E13.5, features = goi,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/hindlimb_E13.5.UMAP.goi.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()

# UMAP of goi_celltype
out_file <- paste(output_folder, "/hindlimb_E13.5.UMAP.goi_celltype.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi_celltype) {
  if (gene %in% rownames(so_hindlimb_E13.5[["SCT"]])) {
    p <- FeaturePlot(so_hindlimb_E13.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi_celltype
p <- FeaturePlot(so_hindlimb_E13.5, features = goi_celltype,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/hindlimb_E13.5.UMAP.goi_celltype.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()

# UMAP of goi_spatial
out_file <- paste(output_folder, "/hindlimb_E13.5.UMAP.goi_spatial.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi_spatial) {
  if (gene %in% rownames(so_hindlimb_E13.5[["SCT"]])) {
    p <- FeaturePlot(so_hindlimb_E13.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi_spatial
p <- FeaturePlot(so_hindlimb_E13.5, features = goi_spatial,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/hindlimb_E13.5.UMAP.goi_spatial.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()


##### hindlimb E15.5
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE15.5/"

# UMAP of goi
out_file <- paste(output_folder, "/hindlimb_E15.5.UMAP.goi.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi) {
  if (gene %in% rownames(so_hindlimb_E15.5[["SCT"]])) {
    p <- FeaturePlot(so_hindlimb_E15.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi
p <- FeaturePlot(so_hindlimb_E15.5, features = goi,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/hindlimb_E15.5.UMAP.goi.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()

# UMAP of goi_celltype
out_file <- paste(output_folder, "/hindlimb_E15.5.UMAP.goi_celltype.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi_celltype) {
  if (gene %in% rownames(so_hindlimb_E15.5[["SCT"]])) {
    p <- FeaturePlot(so_hindlimb_E15.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi_celltype
p <- FeaturePlot(so_hindlimb_E15.5, features = goi_celltype,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/hindlimb_E15.5.UMAP.goi_celltype.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()

# UMAP of goi_spatial
out_file <- paste(output_folder, "/hindlimb_E15.5.UMAP.goi_spatial.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi_spatial) {
  if (gene %in% rownames(so_hindlimb_E15.5[["SCT"]])) {
    p <- FeaturePlot(so_hindlimb_E15.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi_spatial
p <- FeaturePlot(so_hindlimb_E15.5, features = goi_spatial,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/hindlimb_E15.5.UMAP.goi_spatial.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()


##### hindlimb E18.5
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE18.5/"

# UMAP of goi
out_file <- paste(output_folder, "/hindlimb_E18.5.UMAP.goi.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi) {
  if (gene %in% rownames(so_hindlimb_E18.5[["SCT"]])) {
    p <- FeaturePlot(so_hindlimb_E18.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi
p <- FeaturePlot(so_hindlimb_E18.5, features = goi,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/hindlimb_E18.5.UMAP.goi.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()

# UMAP of goi_celltype
out_file <- paste(output_folder, "/hindlimb_E18.5.UMAP.goi_celltype.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi_celltype) {
  if (gene %in% rownames(so_hindlimb_E18.5[["SCT"]])) {
    p <- FeaturePlot(so_hindlimb_E18.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi_celltype
p <- FeaturePlot(so_hindlimb_E18.5, features = goi_celltype,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/hindlimb_E18.5.UMAP.goi_celltype.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()

# UMAP of goi_spatial
out_file <- paste(output_folder, "/hindlimb_E18.5.UMAP.goi_spatial.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi_spatial) {
  if (gene %in% rownames(so_hindlimb_E18.5[["SCT"]])) {
    p <- FeaturePlot(so_hindlimb_E18.5, features = gene)
    print(p)
  } else {
    message(paste("Gene", gene, "not found in dataset. Skipping."))
  }
}
dev.off()

# Multipannel UMAP of goi_spatial
p <- FeaturePlot(so_hindlimb_E18.5, features = goi_spatial,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(output_folder, "/hindlimb_E18.5.UMAP.goi_spatial.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()


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
####################### INTEGRATION OF HINDLIMB DATASETS #######################

# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/hindlimb_integration/"

# Load raw data (pre-processed Seurat objects produced above)
so_hindlimb_E10.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE10.5/so_hindlimb_E10.5.rds")
so_hindlimb_E11.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE11.5/so_hindlimb_E11.5.rds")
so_hindlimb_E13.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE13.5/so_hindlimb_E13.5.rds")
so_hindlimb_E15.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE15.5/so_hindlimb_E15.5.rds")
so_hindlimb_E18.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE18.5/so_hindlimb_E18.5.rds")


# Change the default assay to "SCT"
DefaultAssay(so_hindlimb_E10.5) <- "SCT"
DefaultAssay(so_hindlimb_E11.5) <- "SCT"
DefaultAssay(so_hindlimb_E13.5) <- "SCT"
DefaultAssay(so_hindlimb_E15.5) <- "SCT"
DefaultAssay(so_hindlimb_E18.5) <- "SCT"

# Check memory usage before and after merge (due to Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Identify variable features for both datasets
so_hindlimb_E10.5 <- FindVariableFeatures(so_hindlimb_E10.5, 
                                          selection.method = "vst", 
                                          nfeatures = 2000)
so_hindlimb_E11.5 <- FindVariableFeatures(so_hindlimb_E11.5, 
                                          selection.method = "vst", 
                                          nfeatures = 2000)
so_hindlimb_E13.5 <- FindVariableFeatures(so_hindlimb_E13.5, 
                                          selection.method = "vst", 
                                          nfeatures = 2000)
so_hindlimb_E15.5 <- FindVariableFeatures(so_hindlimb_E15.5, 
                                          selection.method = "vst", 
                                          nfeatures = 2000)
so_hindlimb_E18.5 <- FindVariableFeatures(so_hindlimb_E18.5, 
                                          selection.method = "vst", 
                                          nfeatures = 2000)

# Select integration features
SelectIntegrationFeatures(object.list = list(so_hindlimb_E10.5, 
                                             so_hindlimb_E11.5, 
                                             so_hindlimb_E13.5, 
                                             so_hindlimb_E15.5,
                                             so_hindlimb_E18.5),
                          nfeatures = 2000,
                          verbose = TRUE)

# Step 1: Get the variable features for both datasets
var_features_E10.5 <- so_hindlimb_E10.5@assays[["SCT"]]@var.features
var_features_E11.5 <- so_hindlimb_E11.5@assays[["SCT"]]@var.features
var_features_E13.5 <- so_hindlimb_E13.5@assays[["SCT"]]@var.features
var_features_E15.5 <- so_hindlimb_E15.5@assays[["SCT"]]@var.features
var_features_E18.5 <- so_hindlimb_E18.5@assays[["SCT"]]@var.features

# Step 2: Find the common variable features between the two datasets
common_var_features <- Reduce(intersect, list(var_features_E10.5, 
                                              var_features_E11.5, 
                                              var_features_E13.5,
                                              var_features_E15.5,
                                              var_features_E18.5))

# Step 3: Prepare the objects for integration using the common features
objects <- list(so_hindlimb_E10.5, 
                so_hindlimb_E11.5, 
                so_hindlimb_E13.5,
                so_hindlimb_E15.5,
                so_hindlimb_E18.5)

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
rm(so_hindlimb_E10.5)
rm(so_hindlimb_E11.5)
rm(so_hindlimb_E13.5)
rm(so_hindlimb_E15.5)
rm(so_hindlimb_E18.5)
rm(objects)

# Step 5: Integrate the datasets using the found anchors
so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated <- IntegrateData(anchorset = anchors, 
                                                                      normalization.method = "SCT",
                                                                      dims = 1:10)

# Step 6: Perform scaling and PCA on the integrated data
so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated <- ScaleData(so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated)
so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated <- RunPCA(so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated, verbose = FALSE)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated, ndims = 20)
out_path <- paste(output_folder, "integrated_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5.data.qc.ellbowplot.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

# Store number of principle components in new variable (to be used later)
pca_dim_sel <- 11

# Perform UMAP on integrated data
so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated <- RunUMAP(so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated, 
                                                                dims = 1:pca_dim_sel, 
                                                                reduction = "pca", 
                                                                reduction.name = "umap.integrated")

# Change the default assay to "SCT" (normalized dataset)
DefaultAssay(so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated) <- "SCT"

# Clustering (Leiden) - Seurat v5 should work similarly
so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated <- FindNeighbors(so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated,
                                                                      dims = 1:pca_dim_sel)
so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated <- FindClusters(so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated,
                                                                     resolution = 0.3,
                                                                     algorithm = 4,
                                                                     graph.name = "integrated_snn")

# Visualize datasets as UMAP after Integration
p <- DimPlot(so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated, 
             reduction = "umap.integrated", 
             group.by = c("orig.ident", "seurat_clusters"))
outFile <- paste(output_folder, "/hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated.UMAP.orig.ident.clusters.pdf", sep = "")
pdf(outFile, width = 12, height = 5)
plot(p)
dev.off()

p <- DimPlot(object = so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated,
             reduction = "umap.integrated",
             group.by = 'seurat_clusters',
             split.by = 'orig.ident',
             label = TRUE)
outFile <- paste(output_folder, "/hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated.UMAP.orig.ident.clusters_split.pdf", sep = "")
pdf(outFile, width = 20, height = 5)
plot(p)
dev.off()


# Save the Seurat object
outFile <- paste(output_folder, "so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated.rds", sep = "")
saveRDS(so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated, file = outFile)


######## Visualize as FeaturePlot

# Change the default assay to "SCT"
DefaultAssay(so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated) <- "SCT"

# Plots for goi
outFile <- paste(output_folder,
                 "/hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated.UMAP.goi.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 20, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi) {
  if (gene %in% rownames(so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated, 
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

# Plots for goi_celltype
outFile <- paste(output_folder,
                 "/hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated.UMAP.goi_celltype.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 20, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi_celltype) {
  if (gene %in% rownames(so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated, 
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

# Plots for goi_spatial
outFile <- paste(output_folder,
                 "/hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated.UMAP.goi_spatial.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 20, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi_spatial) {
  if (gene %in% rownames(so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_hindlimb_E10.5_E11.5_E13.5_E15.5_E18.5_integrated, 
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


################################################################################
####################### INTEGRATION OF MANDIBLE DATASETS #######################

# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/hindlimb_integration/"

# Load raw data (pre-processed Seurat objects produced above)
so_mandible_E9.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE9.5/so_mandible_E9.5.rds")
so_mandible_E10.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE10.5/so_mandible_E10.5.rds")
so_mandible_E11.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE11.5/so_mandible_E11.5.rds")


# Change the default assay to "SCT"
DefaultAssay(so_mandible_E9.5) <- "SCT"
DefaultAssay(so_mandible_E10.5) <- "SCT"
DefaultAssay(so_mandible_E11.5) <- "SCT"

# Check memory usage before and after merge (due to Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Identify variable features for both datasets
so_mandible_E9.5 <- FindVariableFeatures(so_mandible_E9.5, 
                                         selection.method = "vst", 
                                         nfeatures = 2000)
so_mandible_E10.5 <- FindVariableFeatures(so_mandible_E10.5, 
                                          selection.method = "vst", 
                                          nfeatures = 2000)
so_mandible_E11.5 <- FindVariableFeatures(so_mandible_E11.5, 
                                          selection.method = "vst", 
                                          nfeatures = 2000)


# Select integration features
SelectIntegrationFeatures(object.list = list(so_mandible_E9.5, 
                                             so_mandible_E10.5, 
                                             so_mandible_E11.5),
                          nfeatures = 2000,
                          verbose = TRUE)

# Step 1: Get the variable features for both datasets
var_features_E9.5 <- so_mandible_E9.5@assays[["SCT"]]@var.features
var_features_E10.5 <- so_mandible_E10.5@assays[["SCT"]]@var.features
var_features_E11.5 <- so_mandible_E11.5@assays[["SCT"]]@var.features

# Step 2: Find the common variable features between the two datasets
common_var_features <- Reduce(intersect, list(var_features_E9.5, 
                                              var_features_E10.5, 
                                              var_features_E11.5))

# Step 3: Prepare the objects for integration using the common features
objects <- list(so_mandible_E9.5, 
                so_mandible_E10.5, 
                so_mandible_E11.5)

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
rm(so_mandible_E9.5)
rm(so_mandible_E10.5)
rm(so_mandible_E11.5)
rm(objects)

# Step 5: Integrate the datasets using the found anchors
so_mandible_E9.5_E10.5_E11.5_integrated <- IntegrateData(anchorset = anchors, 
                                                         normalization.method = "SCT",
                                                         dims = 1:10)

# Step 6: Perform scaling and PCA on the integrated data
so_mandible_E9.5_E10.5_E11.5_integrated <- ScaleData(so_mandible_E9.5_E10.5_E11.5_integrated)
so_mandible_E9.5_E10.5_E11.5_integrated <- RunPCA(so_mandible_E9.5_E10.5_E11.5_integrated, verbose = FALSE)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_mandible_E9.5_E10.5_E11.5_integrated, ndims = 20)
out_path <- paste(output_folder, "integrated_mandible_E9.5_E10.5_E11.5.data.qc.ellbowplot.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

# Store number of principle components in new variable (to be used later)
pca_dim_sel <- 8

# Perform UMAP on integrated data
so_mandible_E9.5_E10.5_E11.5_integrated <- RunUMAP(so_mandible_E9.5_E10.5_E11.5_integrated, 
                                                   dims = 1:pca_dim_sel, 
                                                   reduction = "pca", 
                                                   reduction.name = "umap.integrated")

# Change the default assay to "SCT" (normalized dataset)
DefaultAssay(so_mandible_E9.5_E10.5_E11.5_integrated) <- "SCT"

# Clustering (Leiden) - Seurat v5 should work similarly
so_mandible_E9.5_E10.5_E11.5_integrated <- FindNeighbors(so_mandible_E9.5_E10.5_E11.5_integrated,
                                                         dims = 1:pca_dim_sel)
so_mandible_E9.5_E10.5_E11.5_integrated <- FindClusters(so_mandible_E9.5_E10.5_E11.5_integrated,
                                                        resolution = 0.5,
                                                        algorithm = 4,
                                                        graph.name = "integrated_snn")

# Visualize datasets as UMAP after Integration
p <- DimPlot(so_mandible_E9.5_E10.5_E11.5_integrated, 
             reduction = "umap.integrated", 
             group.by = c("orig.ident", "seurat_clusters"))
outFile <- paste(output_folder, "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.orig.ident.clusters.pdf", sep = "")
pdf(outFile, width = 12, height = 5)
plot(p)
dev.off()

p <- DimPlot(object = so_mandible_E9.5_E10.5_E11.5_integrated,
             reduction = "umap.integrated",
             group.by = 'seurat_clusters',
             split.by = 'orig.ident',
             label = TRUE)
outFile <- paste(output_folder, "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.orig.ident.clusters_split.pdf", sep = "")
pdf(outFile, width = 20, height = 5)
plot(p)
dev.off()


# Save the Seurat object
outFile <- paste(output_folder, "so_mandible_E9.5_E10.5_E11.5_integrated.rds", sep = "")
saveRDS(so_mandible_E9.5_E10.5_E11.5_integrated, file = outFile)


######## Visualize as FeaturePlot

# Change the default assay to "SCT"
DefaultAssay(so_mandible_E9.5_E10.5_E11.5_integrated) <- "SCT"

# Plots for goi
outFile <- paste(output_folder,
                 "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.goi.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 20, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi) {
  if (gene %in% rownames(so_mandible_E9.5_E10.5_E11.5_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_mandible_E9.5_E10.5_E11.5_integrated, 
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

# Plots for goi_celltype
outFile <- paste(output_folder,
                 "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.goi_celltype.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 20, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi_celltype) {
  if (gene %in% rownames(so_mandible_E9.5_E10.5_E11.5_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_mandible_E9.5_E10.5_E11.5_integrated, 
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

# Plots for goi_spatial
outFile <- paste(output_folder,
                 "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.goi_spatial.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 20, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi_spatial) {
  if (gene %in% rownames(so_mandible_E9.5_E10.5_E11.5_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_mandible_E9.5_E10.5_E11.5_integrated, 
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


################################################################################
##### INTEGRATION OF MANDIBLE E9.5-E11.5 WITH HINDLIMB E10.5-E13.5 DATASETS ####

# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/hindlimb_integration/"

# Load raw data (pre-processed Seurat objects produced above)
so_mandible_E9.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE9.5/so_mandible_E9.5.rds")
so_mandible_E10.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE10.5/so_mandible_E10.5.rds")
so_mandible_E11.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE11.5/so_mandible_E11.5.rds")
so_hindlimb_E10.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE10.5/so_hindlimb_E10.5.rds")
so_hindlimb_E11.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE11.5/so_hindlimb_E11.5.rds")
so_hindlimb_E13.5 <- readRDS("~/Documents/postdoc/bioinformatics/results/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE13.5/so_hindlimb_E13.5.rds")

# Change the default assay to "SCT"
DefaultAssay(so_mandible_E9.5) <- "SCT"
DefaultAssay(so_mandible_E10.5) <- "SCT"
DefaultAssay(so_mandible_E11.5) <- "SCT"
DefaultAssay(so_hindlimb_E10.5) <- "SCT"
DefaultAssay(so_hindlimb_E11.5) <- "SCT"
DefaultAssay(so_hindlimb_E13.5) <- "SCT"

# Check memory usage before and after merge (due to Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Identify variable features for both datasets
so_mandible_E9.5 <- FindVariableFeatures(so_mandible_E9.5, 
                                         selection.method = "vst", 
                                         nfeatures = 2000)
so_mandible_E10.5 <- FindVariableFeatures(so_mandible_E10.5, 
                                          selection.method = "vst", 
                                          nfeatures = 2000)
so_mandible_E11.5 <- FindVariableFeatures(so_mandible_E11.5, 
                                          selection.method = "vst", 
                                          nfeatures = 2000)
so_hindlimb_E10.5 <- FindVariableFeatures(so_hindlimb_E10.5, 
                                          selection.method = "vst", 
                                          nfeatures = 2000)
so_hindlimb_E11.5 <- FindVariableFeatures(so_hindlimb_E11.5, 
                                          selection.method = "vst", 
                                          nfeatures = 2000)
so_hindlimb_E13.5 <- FindVariableFeatures(so_hindlimb_E13.5, 
                                          selection.method = "vst", 
                                          nfeatures = 2000)


# Select integration features
SelectIntegrationFeatures(object.list = list(so_mandible_E9.5, 
                                             so_mandible_E10.5, 
                                             so_mandible_E11.5,
                                             so_hindlimb_E10.5, 
                                             so_hindlimb_E11.5, 
                                             so_hindlimb_E13.5),
                          nfeatures = 2000,
                          verbose = TRUE)

# Step 1: Get the variable features for both datasets
var_features_E9.5 <- so_mandible_E9.5@assays[["SCT"]]@var.features
var_features_E10.5 <- so_mandible_E10.5@assays[["SCT"]]@var.features
var_features_E11.5 <- so_mandible_E11.5@assays[["SCT"]]@var.features
var_features_E10.5 <- so_hindlimb_E10.5@assays[["SCT"]]@var.features
var_features_E11.5 <- so_hindlimb_E11.5@assays[["SCT"]]@var.features
var_features_E13.5 <- so_hindlimb_E13.5@assays[["SCT"]]@var.features

# Step 2: Find the common variable features between the two datasets
common_var_features <- Reduce(intersect, list(var_features_E9.5, 
                                              var_features_E10.5, 
                                              var_features_E11.5,
                                              var_features_E10.5, 
                                              var_features_E11.5, 
                                              var_features_E13.5))

# Step 3: Prepare the objects for integration using the common features
objects <- list(so_mandible_E9.5, 
                so_mandible_E10.5, 
                so_mandible_E11.5,
                so_hindlimb_E10.5, 
                so_hindlimb_E11.5, 
                so_hindlimb_E13.5)

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
rm(so_mandible_E9.5)
rm(so_mandible_E10.5)
rm(so_mandible_E11.5)
rm(so_hindlimb_E10.5)
rm(so_hindlimb_E11.5)
rm(so_hindlimb_E13.5)
rm(objects)

# Step 5: Integrate the datasets using the found anchors
so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated <- IntegrateData(anchorset = anchors, 
                                                                                    normalization.method = "SCT",
                                                                                    dims = 1:10)

# Step 6: Perform scaling and PCA on the integrated data
so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated <- ScaleData(so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated)
so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated <- RunPCA(so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated, verbose = FALSE)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated, ndims = 20)
out_path <- paste(output_folder, "integrated_mandible_E9.5-E11.5_hindlimb_E10.5-E13.5.data.qc.ellbowplot.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

# Store number of principle components in new variable (to be used later)
pca_dim_sel <- 8

# Perform UMAP on integrated data
so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated <- RunUMAP(so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated, 
                                                                              dims = 1:pca_dim_sel, 
                                                                              reduction = "pca", 
                                                                              reduction.name = "umap.integrated")

# Change the default assay to "SCT" (normalized dataset)
DefaultAssay(so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated) <- "SCT"

# Clustering (Leiden) - Seurat v5 should work similarly
so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated <- FindNeighbors(so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated,
                                                                                    dims = 1:pca_dim_sel)
so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated <- FindClusters(so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated,
                                                                                   resolution = 0.4,
                                                                                   algorithm = 4,
                                                                                   graph.name = "integrated_snn")

# Visualize datasets as UMAP after Integration
p <- DimPlot(so_so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated, 
             reduction = "umap.integrated", 
             group.by = c("orig.ident", "seurat_clusters"))
outFile <- paste(output_folder, "/integrated_mandible_E9.5-E11.5_hindlimb_E10.5-E13.5.UMAP.orig.ident.clusters.pdf", sep = "")
pdf(outFile, width = 12, height = 5)
plot(p)
dev.off()

p <- DimPlot(object = so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated,
             reduction = "umap.integrated",
             group.by = 'seurat_clusters',
             split.by = 'orig.ident',
             label = TRUE)
outFile <- paste(output_folder, "/integrated_mandible_E9.5-E11.5_hindlimb_E10.5-E13.5.UMAP.orig.ident.clusters_split.pdf", sep = "")
pdf(outFile, width = 20, height = 5)
plot(p)
dev.off()


# Save the Seurat object
outFile <- paste(output_folder, "so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated.rds", sep = "")
saveRDS(so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated, file = outFile)


######## Visualize as FeaturePlot
# Change the default assay to "SCT"
DefaultAssay(so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated) <- "SCT"

# Set the desired order of orig.ident
so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated$orig.ident <- factor(
  so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated$orig.ident,
  levels = c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5", 
             "hindlimb_E10.5", "hindlimb_E11.5", "hindlimb_E13.5")
)

### UMAP Plots for goi
outFile <- paste(output_folder,
                 "/integrated_mandible_E9.5-E11.5_hindlimb_E10.5-E13.5.UMAP.goi.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 25, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi) {
  if (gene %in% rownames(so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated, 
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

# Violin Plot for goi
so <- so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated

# Ensure clusters are identities
Idents(so) <- "seurat_clusters"

# Create composite cluster + sample identity
so$cluster_identity <- paste0("Cluster_", Idents(so), "_", so$orig.ident)

# Get all unique combinations, sorted by desired sample order and cluster
sample_order <- c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5", 
                  "hindlimb_E10.5", "hindlimb_E11.5", "hindlimb_E13.5")

# Build desired factor levels
all_combos <- expand.grid(
  cluster = sort(unique(Idents(so))),
  sample = sample_order
)
cluster_levels <- paste0("Cluster_", all_combos$cluster, "_", all_combos$sample)

# Set the cluster_identity as a factor with desired order
so$cluster_identity <- factor(so$cluster_identity, levels = cluster_levels)

outFile <- paste(output_folder,
                 "/integrated_mandible_E9.5-E11.5_hindlimb_E10.5-E13.5.Violin.goi.by.cluster.and.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 25, height = 6)

for (gene in goi) {
  if (gene %in% rownames(so)) {
    p <- VlnPlot(
      so,
      features = gene,
      group.by = "cluster_identity",
      pt.size = 0.1
    ) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(p)
  } else {
    message(paste("Gene not found in data:", gene))
  }
}

dev.off()


### UMAP Plots for goi_celltype
outFile <- paste(output_folder,
                 "/integrated_mandible_E9.5-E11.5_hindlimb_E10.5-E13.5.UMAP.goi_celltype.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 25, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi_celltype) {
  if (gene %in% rownames(so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated, 
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

### UMAP Plots for goi_spatial
outFile <- paste(output_folder,
                 "/integrated_mandible_E9.5-E11.5_hindlimb_E10.5-E13.5.UMAP.goi_spatial.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 25, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi_spatial) {
  if (gene %in% rownames(so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_mandible_E9.5_E10.5_E11.5_hindlimb_E10.5_E11.5_E13.5_integrated, 
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


################################################################################
###################### INTEGRATION OF DATASETS ON WYNTON #######################

# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)

# Set output path
output_folder <- "/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/hindlimb+midface_integration/analysis/"

# Load Seurat objects from each pre-analysis folder
so_hindlimb_E10.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE10.5/so_hindlimb_E10.5.rds")
so_hindlimb_E11.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE11.5/so_hindlimb_E11.5.rds")
so_hindlimb_E12.5_autopod <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE12.5_autopod/so_autopod_E12.5.rds")
so_hindlimb_E12.5_stylopod <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE12.5_stylopod/so_stylopod_E12.5.rds")
so_hindlimb_E12.5_stylopod_zeugopod <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE12.5_stylopod_zeugopod/so_stylopod_zeugopod_E12.5.rds")
so_hindlimb_E13.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE13.5/so_hindlimb_E13.5.rds")
so_hindlimb_E13.5_autopod <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE13.5_autopod/so_autopod_E13.5.rds")
so_hindlimb_E15.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE15.5/so_hindlimb_E15.5.rds")
so_hindlimb_E18.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE18.5/so_hindlimb_E18.5.rds")

#so_midface_E9.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre_analysis_midfaceE9.5/so_midface_E9.5.rds")
#so_midface_E10.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre_analysis_midfaceE10.5/so_midface_E10.5.rds")
#so_midface_E11.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre_analysis_midfaceE11.5/so_midface_E11.5.rds")

so_mandible_E9.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE9.5/so_mandible_E9.5.rds")
so_mandible_E10.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE10.5/so_mandible_E10.5.rds")
so_mandible_E11.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_mandibleE11.5/so_mandible_E11.5.rds")

# Change the default assay to "SCT"
DefaultAssay(so_hindlimb_E10.5) <- "SCT"
DefaultAssay(so_hindlimb_E11.5) <- "SCT"
DefaultAssay(so_hindlimb_E12.5_autopod) <- "SCT"
DefaultAssay(so_hindlimb_E12.5_stylopod) <- "SCT"
DefaultAssay(so_hindlimb_E12.5_stylopod_zeugopod) <- "SCT"
DefaultAssay(so_hindlimb_E13.5) <- "SCT"
DefaultAssay(so_hindlimb_E13.5_autopod) <- "SCT"
DefaultAssay(so_hindlimb_E15.5) <- "SCT"
DefaultAssay(so_hindlimb_E18.5) <- "SCT"

#DefaultAssay(so_midface_E9.5) <- "SCT"
#DefaultAssay(so_midface_E10.5) <- "SCT"
#DefaultAssay(so_midface_E11.5) <- "SCT"

DefaultAssay(so_mandible_E9.5) <- "SCT"
DefaultAssay(so_mandible_E10.5) <- "SCT"
DefaultAssay(so_mandible_E11.5) <- "SCT"

# Check memory usage before and after merge (due to Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Identify variable features for all datasets
so_hindlimb_E10.5 <- FindVariableFeatures(so_hindlimb_E10.5, selection.method = "vst", nfeatures = 2000)
so_hindlimb_E11.5 <- FindVariableFeatures(so_hindlimb_E11.5, selection.method = "vst", nfeatures = 2000)
so_hindlimb_E12.5_autopod <- FindVariableFeatures(so_hindlimb_E12.5_autopod, selection.method = "vst", nfeatures = 2000)
so_hindlimb_E12.5_stylopod <- FindVariableFeatures(so_hindlimb_E12.5_stylopod, selection.method = "vst", nfeatures = 2000)
so_hindlimb_E12.5_stylopod_zeugopod <- FindVariableFeatures(so_hindlimb_E12.5_stylopod_zeugopod, selection.method = "vst", nfeatures = 2000)
so_hindlimb_E13.5 <- FindVariableFeatures(so_hindlimb_E13.5, selection.method = "vst", nfeatures = 2000)
so_hindlimb_E13.5_autopod <- FindVariableFeatures(so_hindlimb_E13.5_autopod, selection.method = "vst", nfeatures = 2000)
so_hindlimb_E15.5 <- FindVariableFeatures(so_hindlimb_E15.5, selection.method = "vst", nfeatures = 2000)
so_hindlimb_E18.5 <- FindVariableFeatures(so_hindlimb_E18.5, selection.method = "vst", nfeatures = 2000)

#so_midface_E9.5 <- FindVariableFeatures(so_midface_E9.5, selection.method = "vst", nfeatures = 2000)
#so_midface_E10.5 <- FindVariableFeatures(so_midface_E10.5, selection.method = "vst", nfeatures = 2000)
#so_midface_E11.5 <- FindVariableFeatures(so_midface_E11.5, selection.method = "vst", nfeatures = 2000)

so_mandible_E9.5 <- FindVariableFeatures(so_mandible_E9.5, selection.method = "vst", nfeatures = 2000)
so_mandible_E10.5 <- FindVariableFeatures(so_mandible_E10.5, selection.method = "vst", nfeatures = 2000)
so_mandible_E11.5 <- FindVariableFeatures(so_mandible_E11.5, selection.method = "vst", nfeatures = 2000)

# Select integration features
SelectIntegrationFeatures(object.list = list(so_hindlimb_E10.5,
                                             so_hindlimb_E11.5,
                                             so_hindlimb_E12.5_autopod,
                                             so_hindlimb_E12.5_stylopod,
                                             so_hindlimb_E12.5_stylopod_zeugopod,
                                             so_hindlimb_E13.5,
                                             so_hindlimb_E13.5_autopod,
                                             so_hindlimb_E15.5,
                                             so_hindlimb_E18.5,
                                             #so_midface_E9.5,
                                             #so_midface_E10.5,
                                             #so_midface_E11.5,
                                             so_mandible_E9.5,
                                             so_mandible_E10.5,
                                             so_mandible_E11.5),
                          nfeatures = 2000,
                          verbose = TRUE)

# Step 1: Get the variable features for all datasets
var_features_hindlimb_E10.5 <- so_hindlimb_E10.5@assays[["SCT"]]@var.features
var_features_hindlimb_E11.5 <- so_hindlimb_E11.5@assays[["SCT"]]@var.features
var_features_hindlimb_E12.5_autopod <- so_hindlimb_E12.5_autopod@assays[["SCT"]]@var.features
var_features_hindlimb_E12.5_stylopod <- so_hindlimb_E12.5_stylopod@assays[["SCT"]]@var.features
var_features_hindlimb_E12.5_stylopod_zeugopod <- so_hindlimb_E12.5_stylopod_zeugopod@assays[["SCT"]]@var.features
var_features_hindlimb_E13.5 <- so_hindlimb_E13.5@assays[["SCT"]]@var.features
var_features_hindlimb_E13.5_autopod <- so_hindlimb_E13.5_autopod@assays[["SCT"]]@var.features
var_features_hindlimb_E15.5 <- so_hindlimb_E15.5@assays[["SCT"]]@var.features
var_features_hindlimb_E18.5 <- so_hindlimb_E18.5@assays[["SCT"]]@var.features
var_features_mandible_E9.5 <- so_mandible_E9.5@assays[["SCT"]]@var.features
var_features_mandible_E10.5 <- so_mandible_E10.5@assays[["SCT"]]@var.features
var_features_mandible_E11.5 <- so_mandible_E11.5@assays[["SCT"]]@var.features

# Step 2: Find the common variable features between all datasets
common_var_features <- Reduce(intersect, list(var_features_hindlimb_E10.5,
                                              var_features_hindlimb_E11.5,
                                              var_features_hindlimb_E12.5_autopod,
                                              var_features_hindlimb_E12.5_stylopod,
                                              var_features_hindlimb_E12.5_stylopod_zeugopod,
                                              var_features_hindlimb_E13.5,
                                              var_features_hindlimb_E13.5_autopod,
                                              var_features_hindlimb_E15.5,
                                              var_features_hindlimb_E18.5,
                                              var_features_mandible_E9.5,
                                              var_features_mandible_E10.5,
                                              var_features_mandible_E11.5))

# Step 3: Prepare the objects for integration using the common features
objects <- list(so_hindlimb_E10.5,
                so_hindlimb_E11.5,
                so_hindlimb_E12.5_autopod,
                so_hindlimb_E12.5_stylopod,
                so_hindlimb_E12.5_stylopod_zeugopod,
                so_hindlimb_E13.5,
                so_hindlimb_E13.5_autopod,
                so_hindlimb_E15.5,
                so_hindlimb_E18.5,
                so_mandible_E9.5,
                so_mandible_E10.5,
                so_mandible_E11.5)

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

# Step 5: Integrate the datasets using the found anchors
so_mandible_hindlimb_integrated <- IntegrateData(anchorset = anchors, 
                                                 normalization.method = "SCT",
                                                 dims = 1:10)

# Step 6: Perform scaling and PCA on the integrated data
so_mandible_hindlimb_integrated <- ScaleData(so_mandible_hindlimb_integrated)
so_mandible_hindlimb_integrated <- RunPCA(so_mandible_hindlimb_integrated, verbose = FALSE)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_mandible_hindlimb_integrated, ndims = 20)
out_path <- paste(output_folder, "integrated_mandible_hindlimb.data.qc.ellbowplot.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

# Store number of principle components in new variable (to be used later)
pca_dim_sel <- 8

# Perform UMAP on integrated data
so_mandible_hindlimb_integrated <- RunUMAP(so_mandible_hindlimb_integrated,
                                           dims = 1:pca_dim_sel, 
                                           reduction = "pca", 
                                           reduction.name = "umap.integrated")

# Change the default assay to "SCT" (normalized dataset)
DefaultAssay(so_mandible_hindlimb_integrated) <- "SCT"

# Clustering (Leiden) - Seurat v5 should work similarly
so_mandible_hindlimb_integrated <- FindNeighbors(so_mandible_hindlimb_integrated,
                                                 dims = 1:pca_dim_sel)
so_mandible_hindlimb_integrated <- FindClusters(so_mandible_hindlimb_integrated,
                                                resolution = 0.3,
                                                algorithm = 4,
                                                graph.name = "integrated_snn")

# Visualize datasets as UMAP after Integration
p <- DimPlot(so_mandible_hindlimb_integrated, 
             reduction = "umap.integrated", 
             group.by = c("orig.ident", "seurat_clusters"))
outFile <- paste(output_folder, "/integrated_mandible_hindlimb.UMAP.orig.ident.clusters.pdf", sep = "")
pdf(outFile, width = 12, height = 5)
plot(p)
dev.off()

p <- DimPlot(object = so_mandible_hindlimb_integrated,
             reduction = "umap.integrated",
             group.by = 'seurat_clusters',
             split.by = 'orig.ident',
             label = TRUE)
outFile <- paste(output_folder, "/integrated_mandible_hindlimb.UMAP.orig.ident.clusters_split.pdf", sep = "")
pdf(outFile, width = 50, height = 5)
plot(p)
dev.off()


# Save the Seurat object
outFile <- paste(output_folder, "so_mandible_hindlimb_integrated.rds", sep = "")
saveRDS(so_mandible_hindlimb_integrated, file = outFile)

# Load data
so_mandible_hindlimb_integrated <- readRDS("/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_hindlimb+mandible_integration/analysis/so_mandible_hindlimb_integrated.rds")


######## Visualizations

### GOI for TALE-HD/Hand2/Hox cluster, cell type and spatial distribution

# TALE-HD, Hand2 and Hox expression
goi <- c("Pbx1", "Pbx2", "Pbx3", "Hand2", "Irx3","Irx5", "Meis1", "Meis2",      # TALE-HD and Hand2
         "Hoxa1", "Hoxa2", "Hoxa3", "Hoxa4", "Hoxa5", "Hoxa6", "Hoxa7",         # Hox cluster genes
         "Hoxa9", "Hoxa10", "Hoxa11", "Hoxa13", "Hoxd13")

# Celltype-clusters via marker genes
goi_celltype <- c("Lin28a", "Sall4", "Fzd7",                                    # undifferentiated cells; Fernandez-Guerrero et al., 2021
                  "Bmp2", "Col2a1", "Sox9",                                     # chondrogenic differentiation and limb skeletal, digit and joint morphogenesis; Fernandez-Guerrero et al., 2021
                  "Runx2", "Dlx5", "Sp7",                                       # osteoblast progenitor and differentiation marker 
                  "Alx3", "Alx4",                                               # proximal anterior mesoderm; Fernandez-Guerrero et al., 2021                                             # autopod-associated genes; Fernandez-Guerrero et al., 2021
                  "Asb4",                                                       # distal anterior mesoderm; Fernandez-Guerrero et al., 2021 
                  "Krt8", "Krt14", "Krt15",                                     # epidermis/ epidermal keratins; Kelly et al., 2020; Fernandez-Guerrero et al., 2021 
                  "Cdh5",                                                       # vasculature @E11.5; Kelly et al., 2020
                  "Lyz2",                                                       # blood @E11.5; Kelly et al., 2020
                  "Bcan", "Dsg2", "Esrp1"                                       # epithelial markers; Fernandez-Guerrero et al., 2021 
)

# Spatial distribution via marker genes
goi_spatial <- c("Hoxa9", "Hoxd9", "Shox2",                                     # stylopod-associated genes; Fernandez-Guerrero et al., 2021
                 "Hoxa11",                                                      # zeugopod-associated genes; Fernandez-Guerrero et al., 2021
                 "Sall1", "Sall3", "Hoxa13",                                    # autopod-associated genes; Fernandez-Guerrero et al., 2021
                 "Hoxd11", "Hoxd12", "Hoxd13",                                  # digit patterning; Fernandez-Guerrero et al., 2021
                 "Fgf8", "Wwc1", "Pdgfa"                                        # AER; Fernandez-Guerrero et al., 2021
)


# Change the default assay to "SCT"
DefaultAssay(so_mandible_hindlimb_integrated) <- "SCT"

# Set the desired order of orig.ident
so_mandible_hindlimb_integrated$orig.ident <- factor(
  so_mandible_hindlimb_integrated$orig.ident,
  levels = c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5", 
             "hindlimb_E10.5", "hindlimb_E11.5", "stylopod_zeugopod_E12.5",
             "stylopod_E12.5", "autopod_E12.5", 
             "hindlimb_E13.5", "autopod_E13.5", 
             "hindlimb_E15.5", "hindlimb_E18.5")
)

### UMAP Plots for goi
outFile <- paste(output_folder,
                 "/integrated_mandible_hindlimb.UMAP.goi.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 50, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi) {
  if (gene %in% rownames(so_mandible_hindlimb_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_mandible_hindlimb_integrated, 
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

# Violin Plot for goi
so <- so_mandible_hindlimb_integrated

# Ensure clusters are identities
Idents(so) <- "seurat_clusters"

# Create composite cluster + sample identity
so$cluster_identity <- paste0("Cluster_", Idents(so), "_", so$orig.ident)

# Get all unique combinations, sorted by desired sample order and cluster
sample_order <- c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5", 
                  "hindlimb_E10.5", "hindlimb_E11.5", "stylopod_zeugopod_E12.5",
                  "stylopod_E12.5", "autopod_E12.5", 
                  "hindlimb_E13.5", "autopod_E13.5", 
                  "hindlimb_E15.5", "hindlimb_E18.5")

# Build desired factor levels
all_combos <- expand.grid(
  cluster = sort(unique(Idents(so))),
  sample = sample_order
)
cluster_levels <- paste0("Cluster_", all_combos$cluster, "_", all_combos$sample)

# Set the cluster_identity as a factor with desired order
so$cluster_identity <- factor(so$cluster_identity, levels = cluster_levels)

outFile <- paste(output_folder,
                 "/integrated_mandible_hindlimb.Violin.goi.by.cluster.and.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 50, height = 6)

for (gene in goi) {
  if (gene %in% rownames(so)) {
    p <- VlnPlot(
      so,
      features = gene,
      group.by = "cluster_identity",
      pt.size = 0.1
    ) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(p)
  } else {
    message(paste("Gene not found in data:", gene))
  }
}

dev.off()


### UMAP Plots for goi_celltype
outFile <- paste(output_folder,
                 "/integrated_mandible_hindlimb.UMAP.goi_celltype.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 50, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi_celltype) {
  if (gene %in% rownames(so_mandible_hindlimb_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_mandible_hindlimb_integrated, 
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

### UMAP Plots for goi_spatial
outFile <- paste(output_folder,
                 "/integrated_mandible_hindlimb.UMAP.goi_spatial.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 50, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi_spatial) {
  if (gene %in% rownames(so_mandible_hindlimb_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_mandible_hindlimb_integrated, 
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

################################################################################
################## INTEGRATION OF HINDLIMB DATASETS ON WYNTON ##################

# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)

# Set output path
output_folder <- "/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/hindlimb+midface_integration/analysis/"

# Load Seurat objects from each pre-analysis folder
so_hindlimb_E10.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE10.5/so_hindlimb_E10.5.rds")
so_hindlimb_E11.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE11.5/so_hindlimb_E11.5.rds")
so_hindlimb_E12.5_autopod <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE12.5_autopod/so_autopod_E12.5.rds")
so_hindlimb_E12.5_stylopod <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE12.5_stylopod/so_stylopod_E12.5.rds")
so_hindlimb_E12.5_stylopod_zeugopod <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE12.5_stylopod_zeugopod/so_stylopod_zeugopod_E12.5.rds")
so_hindlimb_E13.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE13.5/so_hindlimb_E13.5.rds")
so_hindlimb_E13.5_autopod <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE13.5_autopod/so_autopod_E13.5.rds")
so_hindlimb_E15.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE15.5/so_hindlimb_E15.5.rds")
so_hindlimb_E18.5 <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/pre-analysis_hindlimbE18.5/so_hindlimb_E18.5.rds")


# Change the default assay to "SCT"
DefaultAssay(so_hindlimb_E10.5) <- "SCT"
DefaultAssay(so_hindlimb_E11.5) <- "SCT"
DefaultAssay(so_hindlimb_E12.5_autopod) <- "SCT"
DefaultAssay(so_hindlimb_E12.5_stylopod) <- "SCT"
DefaultAssay(so_hindlimb_E12.5_stylopod_zeugopod) <- "SCT"
DefaultAssay(so_hindlimb_E13.5) <- "SCT"
DefaultAssay(so_hindlimb_E13.5_autopod) <- "SCT"
DefaultAssay(so_hindlimb_E15.5) <- "SCT"
DefaultAssay(so_hindlimb_E18.5) <- "SCT"


# Check memory usage before and after merge (due to Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Identify variable features for all datasets
so_hindlimb_E10.5 <- FindVariableFeatures(so_hindlimb_E10.5, selection.method = "vst", nfeatures = 2000)
so_hindlimb_E11.5 <- FindVariableFeatures(so_hindlimb_E11.5, selection.method = "vst", nfeatures = 2000)
so_hindlimb_E12.5_autopod <- FindVariableFeatures(so_hindlimb_E12.5_autopod, selection.method = "vst", nfeatures = 2000)
so_hindlimb_E12.5_stylopod <- FindVariableFeatures(so_hindlimb_E12.5_stylopod, selection.method = "vst", nfeatures = 2000)
so_hindlimb_E12.5_stylopod_zeugopod <- FindVariableFeatures(so_hindlimb_E12.5_stylopod_zeugopod, selection.method = "vst", nfeatures = 2000)
so_hindlimb_E13.5 <- FindVariableFeatures(so_hindlimb_E13.5, selection.method = "vst", nfeatures = 2000)
so_hindlimb_E13.5_autopod <- FindVariableFeatures(so_hindlimb_E13.5_autopod, selection.method = "vst", nfeatures = 2000)
so_hindlimb_E15.5 <- FindVariableFeatures(so_hindlimb_E15.5, selection.method = "vst", nfeatures = 2000)
so_hindlimb_E18.5 <- FindVariableFeatures(so_hindlimb_E18.5, selection.method = "vst", nfeatures = 2000)


# Select integration features
SelectIntegrationFeatures(object.list = list(so_hindlimb_E10.5,
                                             so_hindlimb_E11.5,
                                             so_hindlimb_E12.5_autopod,
                                             so_hindlimb_E12.5_stylopod,
                                             so_hindlimb_E12.5_stylopod_zeugopod,
                                             so_hindlimb_E13.5,
                                             so_hindlimb_E13.5_autopod,
                                             so_hindlimb_E15.5,
                                             so_hindlimb_E18.5),
                          nfeatures = 2000,
                          verbose = TRUE)

# Step 1: Get the variable features for all datasets
var_features_hindlimb_E10.5 <- so_hindlimb_E10.5@assays[["SCT"]]@var.features
var_features_hindlimb_E11.5 <- so_hindlimb_E11.5@assays[["SCT"]]@var.features
var_features_hindlimb_E12.5_autopod <- so_hindlimb_E12.5_autopod@assays[["SCT"]]@var.features
var_features_hindlimb_E12.5_stylopod <- so_hindlimb_E12.5_stylopod@assays[["SCT"]]@var.features
var_features_hindlimb_E12.5_stylopod_zeugopod <- so_hindlimb_E12.5_stylopod_zeugopod@assays[["SCT"]]@var.features
var_features_hindlimb_E13.5 <- so_hindlimb_E13.5@assays[["SCT"]]@var.features
var_features_hindlimb_E13.5_autopod <- so_hindlimb_E13.5_autopod@assays[["SCT"]]@var.features
var_features_hindlimb_E15.5 <- so_hindlimb_E15.5@assays[["SCT"]]@var.features
var_features_hindlimb_E18.5 <- so_hindlimb_E18.5@assays[["SCT"]]@var.features

# Step 2: Find the common variable features between all datasets
common_var_features <- Reduce(intersect, list(var_features_hindlimb_E10.5,
                                              var_features_hindlimb_E11.5,
                                              var_features_hindlimb_E12.5_autopod,
                                              var_features_hindlimb_E12.5_stylopod,
                                              var_features_hindlimb_E12.5_stylopod_zeugopod,
                                              var_features_hindlimb_E13.5,
                                              var_features_hindlimb_E13.5_autopod,
                                              var_features_hindlimb_E15.5,
                                              var_features_hindlimb_E18.5))

# Step 3: Prepare the objects for integration using the common features
objects <- list(so_hindlimb_E10.5,
                so_hindlimb_E11.5,
                so_hindlimb_E12.5_autopod,
                so_hindlimb_E12.5_stylopod,
                so_hindlimb_E12.5_stylopod_zeugopod,
                so_hindlimb_E13.5,
                so_hindlimb_E13.5_autopod,
                so_hindlimb_E15.5,
                so_hindlimb_E18.5)

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

# Step 5: Integrate the datasets using the found anchors
so_hindlimb_integrated <- IntegrateData(anchorset = anchors, 
                                        normalization.method = "SCT",
                                        dims = 1:10)

# Step 6: Perform scaling and PCA on the integrated data
so_hindlimb_integrated <- ScaleData(so_hindlimb_integrated)
so_hindlimb_integrated <- RunPCA(so_hindlimb_integrated, verbose = FALSE)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_hindlimb_integrated, ndims = 20)
out_path <- paste(output_folder, "integrated_hindlimb.data.qc.ellbowplot.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

# Store number of principle components in new variable (to be used later)
pca_dim_sel <- 8

# Perform UMAP on integrated data
so_hindlimb_integrated <- RunUMAP(so_hindlimb_integrated,
                                  dims = 1:pca_dim_sel, 
                                  reduction = "pca", 
                                  reduction.name = "umap.integrated")

# Change the default assay to "SCT" (normalized dataset)
DefaultAssay(so_hindlimb_integrated) <- "SCT"

# Clustering (Leiden) - Seurat v5 should work similarly
so_hindlimb_integrated <- FindNeighbors(so_hindlimb_integrated,
                                        dims = 1:pca_dim_sel)
so_hindlimb_integrated <- FindClusters(so_hindlimb_integrated,
                                       resolution = 0.3,
                                       algorithm = 4,
                                       graph.name = "integrated_snn")

# Visualize datasets as UMAP after Integration
p <- DimPlot(so_hindlimb_integrated, 
             reduction = "umap.integrated", 
             group.by = c("orig.ident", "seurat_clusters"))
outFile <- paste(output_folder, "/integrated_hindlimb.UMAP.orig.ident.clusters.pdf", sep = "")
pdf(outFile, width = 12, height = 5)
plot(p)
dev.off()

p <- DimPlot(object = so_hindlimb_integrated,
             reduction = "umap.integrated",
             group.by = 'seurat_clusters',
             split.by = 'orig.ident',
             label = TRUE)
outFile <- paste(output_folder, "/integrated_hindlimb.UMAP.orig.ident.clusters_split.pdf", sep = "")
pdf(outFile, width = 30, height = 5)
plot(p)
dev.off()

# Save the Seurat object
outFile <- paste(output_folder, "so_hindlimb_integrated.rds", sep = "")
saveRDS(so_hindlimb_integrated, file = outFile)


######## Visualizations
# Load data
so_hindlimb_integrated <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/hindlimb+midface_integration/analysis/so_hindlimb_integrated.rds")

### GOI for TALE-HD/bHLH/Hox cluster, cell type and spatial distribution
# TALE-HD, Hand2 and Hox expression
goi <- c("Pbx1", "Pbx2", "Pbx3", "Pbx4",                                        # TALE-HD, PBC domain
         "Meis1", "Meis2", "Pknox1", "Pknox2",                                  # TALE-HD, MEINOX domain
         "Tcf3", "Tcf4", "Tcf12",                                               # bHLH class I (E-proteins)
         "Hand1", "Hand2",                                                      # bHLH, class II (Hand1/2)
         "Hoxa1", "Hoxa2", "Hoxa3", "Hoxa4", "Hoxa5", "Hoxa6", "Hoxa7",         # Hox cluster genes
         "Hoxa9", "Hoxa10", "Hoxa11", 
         "Hoxa13", "Hoxb13", "Hoxc13", "Hoxd13")

# Celltype-clusters via marker genes
goi_celltype <- c("Lin28a", "Sall4", "Fzd7",                                    # undifferentiated cells; Fernandez-Guerrero et al., 2021
                  "Bmp2", "Col2a1", "Sox9",                                     # chondrogenic differentiation and limb skeletal, digit and joint morphogenesis; Fernandez-Guerrero et al., 2021
                  "Runx2", "Dlx5", "Sp7",                                       # osteoblast progenitor and differentiation marker 
                  "Alx3", "Alx4",                                               # proximal anterior mesoderm; Fernandez-Guerrero et al., 2021                                             # autopod-associated genes; Fernandez-Guerrero et al., 2021
                  "Asb4",                                                       # distal anterior mesoderm; Fernandez-Guerrero et al., 2021 
                  "Krt8", "Krt14", "Krt15",                                     # epidermis/ epidermal keratins; Kelly et al., 2020; Fernandez-Guerrero et al., 2021 
                  "Cdh5",                                                       # vasculature @E11.5; Kelly et al., 2020
                  "Lyz2",                                                       # blood @E11.5; Kelly et al., 2020
                  "Bcan", "Dsg2", "Esrp1"                                       # epithelial markers; Fernandez-Guerrero et al., 2021 
)

# Spatial distribution via marker genes
goi_spatial <- c("Hoxa9", "Hoxd9", "Shox2",                                     # stylopod-associated genes; Fernandez-Guerrero et al., 2021
                 "Hoxa11",                                                      # zeugopod-associated genes; Fernandez-Guerrero et al., 2021
                 "Sall1", "Sall3", "Hoxa13",                                    # autopod-associated genes; Fernandez-Guerrero et al., 2021
                 "Hoxd11", "Hoxd12", "Hoxd13",                                  # digit patterning; Fernandez-Guerrero et al., 2021
                 "Fgf8", "Wwc1", "Pdgfa"                                        # AER; Fernandez-Guerrero et al., 2021
)


# Change the default assay to "SCT"
DefaultAssay(so_hindlimb_integrated) <- "SCT"

# Set the desired order of orig.ident
so_hindlimb_integrated$orig.ident <- factor(
  so_hindlimb_integrated$orig.ident,
  levels = c("hindlimb_E10.5", "hindlimb_E11.5", "stylopod_E12.5",
             "stylopod_zeugopod_E12.5", "autopod_E12.5", 
             "hindlimb_E13.5", "autopod_E13.5", 
             "hindlimb_E15.5", "hindlimb_E18.5")
)

### UMAP Plots for goi
outFile <- paste(output_folder,
                 "/integrated_hindlimb.UMAP.goi.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 30, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi) {
  if (gene %in% rownames(so_hindlimb_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_hindlimb_integrated, 
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

# Violin Plot for goi
so <- so_hindlimb_integrated

# Ensure clusters are identities
Idents(so) <- "seurat_clusters"

# Create composite cluster + sample identity
so$cluster_identity <- paste0("Cluster_", Idents(so), "_", so$orig.ident)

# Get all unique combinations, sorted by desired sample order and cluster
sample_order <- c("hindlimb_E10.5", "hindlimb_E11.5", "stylopod_E12.5",
                  "stylopod_zeugopod_E12.5", "autopod_E12.5", 
                  "hindlimb_E13.5", "autopod_E13.5", 
                  "hindlimb_E15.5", "hindlimb_E18.5")

# Build desired factor levels
all_combos <- expand.grid(
  cluster = sort(unique(Idents(so))),
  sample = sample_order
)
cluster_levels <- paste0("Cluster_", all_combos$cluster, "_", all_combos$sample)

# Set the cluster_identity as a factor with desired order
so$cluster_identity <- factor(so$cluster_identity, levels = cluster_levels)

outFile <- paste(output_folder,
                 "/integrated_hindlimb.Violin.goi.by.cluster.and.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 30, height = 6)

for (gene in goi) {
  if (gene %in% rownames(so)) {
    p <- VlnPlot(
      so,
      features = gene,
      group.by = "cluster_identity",
      pt.size = 0.1
    ) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(p)
  } else {
    message(paste("Gene not found in data:", gene))
  }
}

dev.off()


### UMAP Plots for goi_celltype
outFile <- paste(output_folder,
                 "/integrated_hindlimb.UMAP.goi_celltype.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 30, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi_celltype) {
  if (gene %in% rownames(so_hindlimb_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_hindlimb_integrated, 
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

### UMAP Plots for goi_spatial
outFile <- paste(output_folder,
                 "/integrated_hindlimb.UMAP.goi_spatial.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 30, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi_spatial) {
  if (gene %in% rownames(so_hindlimb_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_hindlimb_integrated, 
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


### ADDITIONAL VISUALIZATOINS

### Upregulated genes in bulk RNA-seq E10.5 of Pbx1/2 or Hand2 mutant
# Genes derived from DESeq2 analysis 05/20/2025
goi_pbx_hand2_mutUP <- c("Gfra2", "LOC102636514", "Mme", "Opcml", "BC068157",
                         "Scara5", "Abca9", "Maob", "C130026L21Rik", "Zfp811", 
                         "Pclo", "Stat2", "Nrk", "Dab1", "Pdzd7", "Sox5", 
                         "Nipal3", "Pcdh19", "Fbln5", "Kctd12", "Cdon", "Rin2", 
                         "Gdf6", "Pgap1", "Srcin1", "Col24a1", "Tex15", 
                         "Agtrap", "Kif16b", "Sepp1", "Pml", "Mxd4", "Map1a", 
                         "Pard3b", "Hif3a", "Ildr2", "Snai2", "Heatr5b", 
                         "Prrx1", "Dab2", "Krt19", "Lamb2", "Npnt", "Wdr19", 
                         "Spats2l", "Crat", "Smarca1", "Ddr1", "Magee1", "Lcp1", 
                         "Pcmtd2", "Kcna6", "Speg", "Mgat3", "Pls3", "Wdr86", 
                         "Epb4.1l3", "Rfx5", "Gprasp2", "Tenm3", "Stox2", 
                         "Nusap1", "Bcorl1", "1600014C10Rik", "Col15a1", 
                         "Igdcc4", "Ptpru", "Wrn", "Kmt2e", "Cxxc5", "C2cd3", 
                         "Sort1", "Wdtc1", "Wars2", "Fut10", "Inadl", "Lhfpl2", 
                         "Sesn1", "Irx5", "Etnk1", "Efnb2", "Cic", "Ift140", 
                         "Mlh3", "Sgol2", "Marf1", "Zfp219", "Ncoa2", "Capn2", 
                         "Stat5b", "Tbc1d9b", "Nedd4", "Erlin2"
)

### UMAP Plots for goi_pbx_hand2_mutUP
outFile <- paste(output_folder,
                 "/integrated_hindlimb.UMAP.goi_pbx_hand2_mutUP.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 30, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi_pbx_hand2_mutUP) {
  if (gene %in% rownames(so_hindlimb_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_hindlimb_integrated, 
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

### Violin Plots for goi_pbx_hand2_mutUP
outFile <- paste(output_folder,
                 "/integrated_hindlimb.Violin.goi_pbx_hand2_mutUP.by.cluster.and.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 30, height = 6)
for (gene in goi_pbx_hand2_mutUP) {
  if (gene %in% rownames(so)) {
    p <- VlnPlot(
      so,
      features = gene,
      group.by = "cluster_identity",
      pt.size = 0.1
    ) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(p)
  } else {
    message(paste("Gene not found in data:", gene))
  }
}
dev.off()

### Downregulated genes in bulk RNA-seq E10.5 of Pbx1/2 or Hand2 mutant
# Genes derived from DESeq2 analysis 05/20/2025
goi_pbx_hand2_mutDOWN <- c("Ewsr1", "Dnajc8", "Pim1", "Rybp", "Rbmx", "Eif1", 
                           "Anapc15", "Cep76", "Prkrip1", "Tipin", "Chuk", 
                           "Rpl10a", "Trim27", "Aif1l", "Fam204a", "Mesdc1", 
                           "Igf2os", "Chsy1", "Strap", "Aars", "Cwf19l1", 
                           "Fjx1", "Prps1", "Guk1", "Adrm1", "Stk39", "Uqcrq", 
                           "2410002F23Rik", "Crcp", "Rpl36", "Atf4", "Rassf1", 
                           "Hnrnpdl", "Herpud1", "Bdh1", "Polr3d", 
                           "1500012F01Rik", "Enc1", "Rbm3", "Snrpa1", "Gas5", 
                           "Snhg8", "Slc25a33", "Snai1", "Mthfd2", "Ifrd1", 
                           "Mat2a", "Wnt11", "Pcolce2", "Snhg5", 
                           "D430041D05Rik", "Dusp4", "Hey1", "Ets2", "Timm17a", 
                           "Trappc5", "Gmip", "Apln", "Serinc5", "Asns", 
                           "Snhg3", "Galnt16", "Gpr4", "Snhg4", "Tspan7", 
                           "Mir17hg", "Mri1", "Mamdc2", "Bves", "Nkx1-1", 
                           "Tmem132e", "Ccser1", "Tulp2", "Clgn", "Leprel1", 
                           "Dlg2", "Dok4", "Slc35d3", "Pdlim3", "Olfm1", 
                           "Krt23", "Ano1", "Nkx1-2", "Gja3", "Cdh22", 
                           "Trim55", "Fgf15")

### UMAP Plots for goi_pbx_hand2_mutDOWN
outFile <- paste(output_folder,
                 "/integrated_hindlimb.UMAP.goi_pbx_hand2_mutDOWN.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 30, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi_pbx_hand2_mutDOWN) {
  if (gene %in% rownames(so_hindlimb_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_hindlimb_integrated, 
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

### Violin Plots for goi_pbx_hand2_mutDOWN
outFile <- paste(output_folder,
                 "/integrated_hindlimb.Violin.goi_pbx_hand2_mutDOWN.by.cluster.and.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 30, height = 6)
for (gene in goi_pbx_hand2_mutDOWN) {
  if (gene %in% rownames(so)) {
    p <- VlnPlot(
      so,
      features = gene,
      group.by = "cluster_identity",
      pt.size = 0.1
    ) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(p)
  } else {
    message(paste("Gene not found in data:", gene))
  }
}
dev.off()


### Genes UPregulated in Pbx1/2 mutants and DOWNregulated in Hand2 mutants
# Genes derived from DESeq2 analysis 05/20/2025
goi_pbxmutUP_hand2mutDOWN <- c("Fto", "Kif26a", "Habp4", "Mycl", "Arrb1", "Prex1", 
                               "Micall2", "Lrrc17", "Crnde", "Cited1", "Lypd6", 
                               "Stom", "Rspo3", "Arhgap8", "Grm4", "Ptch1", "Cntfr", 
                               "Tmem37", "Postn", "Aadat", "Kcne3", "2210039B01Rik", 
                               "Gli1", "Ptch2", "B3gnt3", "Hand2")

### UMAP Plots for goi_pbxmutUP_hand2mutDOWN
outFile <- paste(output_folder,
                 "/integrated_hindlimb.UMAP.goi_pbxmutUP_hand2mutDOWN.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 30, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi_pbxmutUP_hand2mutDOWN) {
  if (gene %in% rownames(so_hindlimb_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_hindlimb_integrated, 
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

### Violin Plots for goi_pbxmutUP_hand2mutDOWN
outFile <- paste(output_folder,
                 "/integrated_hindlimb.Violin.goi_pbxmutUP_hand2mutDOWN.by.cluster.and.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 30, height = 6)
for (gene in goi_pbxmutUP_hand2mutDOWN) {
  if (gene %in% rownames(so)) {
    p <- VlnPlot(
      so,
      features = gene,
      group.by = "cluster_identity",
      pt.size = 0.1
    ) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(p)
  } else {
    message(paste("Gene not found in data:", gene))
  }
}
dev.off()

### Genes DOWNregulated in Pbx1/2 mutants and UPregulated in Hand2 mutants
# Genes derived from DESeq2 analysis 05/20/2025
goi_pbxmutDOWN_hand2mutUP <- c("Alcam", "Map6", "Wif1", "Cux2", "Zic3", "Alx3", 
                               "Lhx9", "Dach2", "Epha4", "Reep2", "Camk2n1", 
                               "Msx1", "Aff3", "Tspan13", "Scn2b", "Rnpep", 
                               "Snhg6", "Snrpb2", "Ehd3", "Pitpnm1", "Usp29", 
                               "Disp1", "Hmgcr", "Vars", "Zhx2", "Nek7")

### UMAP Plots for goi_pbxmutDOWN_hand2mutUP
outFile <- paste(output_folder,
                 "/integrated_hindlimb.UMAP.goi_pbxmutDOWN_hand2mutUP.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 30, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi_pbxmutDOWN_hand2mutUP) {
  if (gene %in% rownames(so_hindlimb_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_hindlimb_integrated, 
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

### Violin Plots for goi_pbxmutDOWN_hand2mutUP
outFile <- paste(output_folder,
                 "/integrated_hindlimb.Violin.goi_pbxmutDOWN_hand2mutUP.by.cluster.and.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 30, height = 6)
for (gene in goi_pbxmutDOWN_hand2mutUP) {
  if (gene %in% rownames(so)) {
    p <- VlnPlot(
      so,
      features = gene,
      group.by = "cluster_identity",
      pt.size = 0.1
    ) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(p)
  } else {
    message(paste("Gene not found in data:", gene))
  }
}
dev.off()

################################################################################
######### Analysis of cell trajectories in integrated hindlimb dataset #########

# Load data
so_hindlimb_integrated <- readRDS("/wynton/home/selleri/veritasnondatur/scRNA-seq/integrated_hindlimbE10.5-18.5_midfaceE9.5-E11.5/hindlimb+midface_integration/analysis/so_hindlimb_integrated.rds")

### Minimal code (ChatGPT)
# Load required libraries
library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(dplyr)

# Your Seurat object
seurat_obj <- so_hindlimb_integrated

#------------------------------------------
# STEP 1: Convert Seurat to Monocle3 CDS
#------------------------------------------

# Use SeuratWrappers to convert
cds <- as.cell_data_set(seurat_obj)

# Optional: copy Seurat cluster info to CDS
colData(cds)$seurat_clusters <- seurat_obj$seurat_clusters

#------------------------------------------
# STEP 2: Add UMAP embeddings from Seurat (optional, recommended)
#------------------------------------------

reducedDims(cds)$UMAP <- Embeddings(seurat_obj, reduction = "umap.integrated")

#------------------------------------------
# STEP 3: Preprocess CDS
#------------------------------------------

# Skip this if Seurat already did PCA and you're using its UMAP
#cds <- preprocess_cds(cds, num_dim = 50)

# Optional: run dimension reduction (skip if reusing Seurat UMAP)
#cds <- reduce_dimension(cds, reduction_method = "UMAP")

#------------------------------------------
# STEP 4: Cluster Cells (optional)
#------------------------------------------

# If you want Monocle to assign clusters:
cds <- cluster_cells(cds)

# Or use Seurat clusters directly (already in colData)

#------------------------------------------
# STEP 5: Learn the Trajectory Graph
#------------------------------------------

cds <- learn_graph(cds)

#------------------------------------------
# STEP 6: Order Cells (pseudotime)
#------------------------------------------

# Interactive mode: manually click root node in plot
# cds <- order_cells(cds)

# OR manually define root cells (e.g., based on known starting cluster)
# Example: use 10 cells from cluster "0" as root
root_cells <- colnames(seurat_obj)[seurat_obj$seurat_clusters == "2"][1:10]
cds <- order_cells(cds, root_cells = root_cells)

# Ensure cell names are matching
common_cells <- intersect(colnames(cds), colnames(seurat_obj))

# Subset both objects to only include common cells
cds <- cds[, common_cells]
seurat_obj <- subset(seurat_obj, cells = common_cells)

# Now transfer Seurat clusters to CDS metadata
cds@colData$seurat_clusters <- seurat_obj$seurat_clusters[common_cells]

# Confirm it's added
table(cds@colData$seurat_clusters)

#------------------------------------------
# STEP 7: Plot the Trajectory
#------------------------------------------

# Open PDF device
outFile <- paste(output_folder,
                 "/integrated_hindlimb.UMAP.trajectory_plot.seurat_clusters.pdf", 
                 sep = "")
pdf(outFile, width = 20, height = 15)  # Adjust size as needed
# Plot the cells
plot_cells(cds, 
           color_cells_by = "seurat_clusters",
           label_groups_by_cluster = TRUE,
           label_leaves = TRUE,
           label_branch_points = TRUE,
           graph_label_size = 6)  # Increase this value for larger text)
# Close the PDF device
dev.off()

# Open PDF device
outFile <- paste(output_folder,
                 "/integrated_hindlimb.UMAP.trajectory_plot.pseudotime.pdf", 
                 sep = "")
pdf(outFile, width = 20, height = 15)  # Adjust size as needed
# Plot colored by pseudotime
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = TRUE,
           label_leaves = TRUE,
           label_branch_points = TRUE,
           graph_label_size = 6)  # Increase this value for larger text
# Close the PDF device
dev.off()

#------------------------------------------
# STEP 8: Add Pseudotime Back to Seurat
#------------------------------------------

seurat_obj$pseudotime <- pseudotime(cds)

# Now pseudotime is available as metadata for downstream Seurat plots
FeaturePlot(seurat_obj, features = "pseudotime", reduction = "umap.integrated")



### Workflow according to Monocle3 documentation 
#[for C.elegans data, needs to be adapted to my context]

# Pre-process the data
mandible_hindlimb_int_trajectory <- new_cell_data_set(so_mandible_hindlimb_integrated)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")

# Reduce dimensionality and visualize the results
cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cell.type")

ciliated_genes <- c("che-1",
                    "hlh-17",
                    "nhr-6",
                    "dmd-6",
                    "ceh-36",
                    "ham-1")

plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

# Cluster your cells
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

# Learn the trajectory graph
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "cell.type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

# Order the cells in pseudotime
plot_cells(cds,
           color_cells_by = "embryo.time.bin",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="130-170"){
  cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

# Subset cells by branch
cds_sub <- choose_graph_segments(cds)

# Working with 3D trajectories
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))

cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")



########################## MAKE PLOTS FOR TALKS ETC. ##########################
# Load data
so_mandible_hindlimb_integrated <- readRDS("/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_hindlimb+mandible_integration/analysis/so_mandible_hindlimb_integrated.rds")

# Set output folder
output_folder <- "/Users/veralaub/Documents/postdoc/presentations/2025-09-018_QDB_talk/plots/"

### For QDB talk 2025-09-18

# Change the default assay to "SCT"
DefaultAssay(so_mandible_hindlimb_integrated) <- "SCT"

# Set the desired order of orig.ident
so_mandible_hindlimb_integrated$orig.ident <- factor(
  so_mandible_hindlimb_integrated$orig.ident,
  levels = c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5", 
             "hindlimb_E10.5", "hindlimb_E11.5", "stylopod_zeugopod_E12.5",
             "stylopod_E12.5", "autopod_E12.5", 
             "hindlimb_E13.5", "autopod_E13.5", 
             "hindlimb_E15.5", "hindlimb_E18.5")
)

#Dimplot
p <- DimPlot(object = so_mandible_hindlimb_integrated,
             reduction = "umap.integrated",
             group.by = 'seurat_clusters',
             split.by = 'orig.ident',
             label = TRUE)
outFile <- paste(output_folder, "/integrated_mandible_hindlimb.UMAP.orig.ident.clusters_split.pdf", sep = "")
pdf(outFile, width = 50, height = 5)
plot(p)
dev.off()

# Set genes to plot
goi <- c("Pbx1", "Pbx2", "Pbx3", "Pbx4", "Hand2",                               # TALE-HD and Hand2
         "Hoxa1", "Hoxa2", "Hoxa3", "Hoxa4", "Hoxa5", "Hoxa6", "Hoxa7",         # Hox cluster genes
         "Hoxa9", "Hoxa10", "Hoxa11", "Hoxa13", "Hoxd13")

### UMAP Plots for goi
outFile <- paste(output_folder,
                 "/integrated_mandible_hindlimb.UMAP.goi_Pbx1-4,Hand2,Hoxa.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 50, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi) {
  if (gene %in% rownames(so_mandible_hindlimb_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_mandible_hindlimb_integrated, 
                     features = gene,
                     reduction = "umap.integrated",
                     split.by = "orig.ident",
                     order = TRUE,
                     pt.size = 2)
    plot(p)
  } else {
    # Print a message for missing genes (optional)
    message(paste("Gene not found in data: ", gene))
  }
}

dev.off()

# FeaturePlots for goi (all integrated datasets together), with contrast raster + axes + improved colors
# Adapted FeaturePlots with contrast raster style
library(ggplot2)

so <- so_mandible_hindlimb_integrated  # convenience short name

# define color palettes and suffixes
palettes <- list(
  blue    = c("lightskyblue1", "steelblue", "darkblue"),
  purple  = c("violet", "mediumorchid", "purple4"),
  green   = c("palegoldenrod", "seagreen3", "darkgreen")
)

# parameters
top_quantile <- 0.90    # top 10% of positives
nonexp_size <- 0.35     
low_size <- 0.7         
high_size <- 1.6        
grey_col <- "grey80"

# explicitly define orig_levels before the loop
orig_levels <- c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5", 
                 "hindlimb_E10.5", "hindlimb_E11.5", "stylopod_zeugopod_E12.5",
                 "stylopod_E12.5", "autopod_E12.5", 
                 "hindlimb_E13.5", "autopod_E13.5", 
                 "hindlimb_E15.5", "hindlimb_E18.5")

for (pal_name in names(palettes)) {
  palette_cols <- palettes[[pal_name]]
  
  outFile <- paste0(
    output_folder, "/integrated_mandible_hindlimb.UMAP.goi_Pbx1-4,Hand2,Hoxa.orig.ident_ContrastRaster_", 
    pal_name, ".pdf"
  )
  
  pdf(outFile, width = 50, height = 5)
  
  for (gene in goi) {
    if (!(gene %in% rownames(so))) {
      message("Gene not found in data: ", gene)
      next
    }
    
    emb <- Embeddings(so, "umap.integrated")
    fetch <- FetchData(so, vars = c(gene, "orig.ident"))
    
    df <- data.frame(
      cell = rownames(fetch),
      UMAP_1 = emb[rownames(fetch), 1],
      UMAP_2 = emb[rownames(fetch), 2],
      expr = as.numeric(fetch[, gene]),
      orig.ident = factor(fetch[, "orig.ident"], levels = orig_levels),
      stringsAsFactors = FALSE
    )
    
    # split non / low / high
    non_idx <- which(df$expr == 0 | is.na(df$expr))
    pos_idx <- which(df$expr > 0)
    if (length(pos_idx) == 0) {
      qtop <- NA
    } else {
      qtop <- as.numeric(quantile(df$expr[pos_idx], probs = top_quantile, na.rm = TRUE))
    }
    
    if (is.na(qtop) || qtop == 0) {
      low_idx <- integer(0)
      high_idx <- pos_idx
    } else {
      low_idx <- which(df$expr > 0 & df$expr <= qtop)
      high_idx <- which(df$expr > qtop)
      if (length(high_idx) == 0 && length(pos_idx) >= 1) {
        qtop2 <- as.numeric(quantile(df$expr[pos_idx], probs = 0.75, na.rm = TRUE))
        low_idx <- which(df$expr > 0 & df$expr <= qtop2)
        high_idx <- which(df$expr > qtop2)
      }
    }
    
    df_non  <- df[non_idx, , drop = FALSE]
    df_low  <- df[low_idx, , drop = FALSE]
    df_high <- df[high_idx, , drop = FALSE]
    
    # plotting
    p <- ggplot() +
      # Non-expressors
      (if (nrow(df_non) > 0) {
        if (has_ggrastr) ggrastr::geom_point_rast(
          data = df_non, aes(x = UMAP_1, y = UMAP_2),
          color = grey_col, size = nonexp_size, raster.dpi = 150
        ) else geom_point(
          data = df_non, aes(x = UMAP_1, y = UMAP_2),
          color = grey_col, size = nonexp_size
        )
      } else NULL) +
      # Low expressors
      (if (nrow(df_low) > 0) {
        if (has_ggrastr) ggrastr::geom_point_rast(
          data = df_low, aes(x = UMAP_1, y = UMAP_2, color = expr),
          size = low_size, raster.dpi = 150
        ) else geom_point(
          data = df_low, aes(x = UMAP_1, y = UMAP_2, color = expr),
          size = low_size
        )
      } else NULL) +
      # High expressors
      (if (nrow(df_high) > 0) {
        if (has_ggrastr) ggrastr::geom_point_rast(
          data = df_high, aes(x = UMAP_1, y = UMAP_2, color = expr),
          size = high_size, raster.dpi = 150
        ) else geom_point(
          data = df_high, aes(x = UMAP_1, y = UMAP_2, color = expr),
          size = high_size
        )
      } else NULL) +
      facet_wrap(~ orig.ident, nrow = 1) +
      scale_color_gradientn(colors = palette_cols, na.value = palette_cols[1]) +
      ggtitle(gene) +
      theme_minimal() +
      labs(x = "umapintegrated_1", y = "umapintegrated_2") +
      theme(
        strip.text = element_text(size = 12),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid = element_blank()
      ) +
      guides(color = guide_colorbar(title = "Expression"))
    
    print(p)
  }
  
  dev.off()
}

### Plots for subset of tissues
library(ggplot2)
has_ggrastr <- requireNamespace("ggrastr", quietly = TRUE)

# subset of interest
subset_levels <- c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5", 
                   "hindlimb_E10.5", "hindlimb_E11.5", "autopod_E12.5",
                   "hindlimb_E13.5", "autopod_E13.5", 
                   "hindlimb_E15.5", "hindlimb_E18.5")

# compute a PDF width that scales with number of samples (adjust multiplier as needed)
n_samples <- length(subset_levels)
pdf_width <- max(10, n_samples * 3.5)  # ~3.5 inches per sample; adjust to taste
pdf_height <- 5

for (pal_name in names(palettes)) {
  palette_cols <- palettes[[pal_name]]
  
  outFile <- file.path(
    output_folder,
    paste0("integrated_mandible_hindlimb.SUBSET.UMAP.goi_ContrastRaster_", pal_name, ".pdf")
  )
  
  pdf(outFile, width = pdf_width, height = pdf_height)
  
  for (gene in goi) {
    if (!(gene %in% rownames(so))) {
      message("Gene not found in data: ", gene)
      next
    }
    
    emb <- Embeddings(so, "umap.integrated")
    fetch <- FetchData(so, vars = c(gene, "orig.ident"))
    
    # build df, keep orig.ident as character first, then filter, then set factor with desired order
    df <- data.frame(
      cell = rownames(fetch),
      UMAP_1 = emb[rownames(fetch), 1],
      UMAP_2 = emb[rownames(fetch), 2],
      expr = as.numeric(fetch[, gene]),
      orig.ident = as.character(fetch[, "orig.ident"]),
      stringsAsFactors = FALSE
    )
    
    # keep only requested subset_levels
    df <- df[df$orig.ident %in% subset_levels, , drop = FALSE]
    
    # if no cells for this gene in subset, skip with a message
    if (nrow(df) == 0) {
      message("No cells for gene ", gene, " in chosen subset. Skipping.")
      next
    }
    
    # enforce factor order (this will keep facets in the order subset_levels,
    # and drop any levels from subset_levels that aren't present in df automatically)
    df$orig.ident <- factor(df$orig.ident, levels = subset_levels)
    
    # split non / low / high (recomputed per gene)
    non_idx <- which(df$expr == 0 | is.na(df$expr))
    pos_idx <- which(df$expr > 0)
    
    if (length(pos_idx) == 0) {
      qtop <- NA
    } else {
      qtop <- as.numeric(quantile(df$expr[pos_idx], probs = top_quantile, na.rm = TRUE))
    }
    
    if (is.na(qtop) || qtop == 0) {
      low_idx <- integer(0)
      high_idx <- pos_idx
    } else {
      low_idx <- which(df$expr > 0 & df$expr <= qtop)
      high_idx <- which(df$expr > qtop)
      if (length(high_idx) == 0 && length(pos_idx) >= 1) {
        qtop2 <- as.numeric(quantile(df$expr[pos_idx], probs = 0.75, na.rm = TRUE))
        low_idx <- which(df$expr > 0 & df$expr <= qtop2)
        high_idx <- which(df$expr > qtop2)
      }
    }
    
    # always produce these df pieces (may be zero-row data.frames)
    df_non  <- df[non_idx, , drop = FALSE]
    df_low  <- df[low_idx, , drop = FALSE]
    df_high <- df[high_idx, , drop = FALSE]
    
    # build plot (use df as base to ensure facet variable is present somewhere)
    p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
      # Non-expressors (grey)
      (if (nrow(df_non) > 0) {
        if (has_ggrastr) ggrastr::geom_point_rast(
          data = df_non, aes(x = UMAP_1, y = UMAP_2),
          color = grey_col, size = nonexp_size, raster.dpi = 150
        ) else geom_point(
          data = df_non, aes(x = UMAP_1, y = UMAP_2),
          color = grey_col, size = nonexp_size
        )
      } else NULL) +
      # Low expressors (color by expr)
      (if (nrow(df_low) > 0) {
        if (has_ggrastr) ggrastr::geom_point_rast(
          data = df_low, aes(x = UMAP_1, y = UMAP_2, color = expr),
          size = low_size, raster.dpi = 150
        ) else geom_point(
          data = df_low, aes(x = UMAP_1, y = UMAP_2, color = expr),
          size = low_size
        )
      } else NULL) +
      # High expressors (color by expr, plotted last to be on top)
      (if (nrow(df_high) > 0) {
        if (has_ggrastr) ggrastr::geom_point_rast(
          data = df_high, aes(x = UMAP_1, y = UMAP_2, color = expr),
          size = high_size, raster.dpi = 150
        ) else geom_point(
          data = df_high, aes(x = UMAP_1, y = UMAP_2, color = expr),
          size = high_size
        )
      } else NULL) +
      facet_wrap(~ orig.ident, nrow = 1, drop = TRUE) +
      scale_color_gradientn(colors = palette_cols, na.value = palette_cols[1]) +
      ggtitle(gene) +
      theme_minimal() +
      labs(x = "UMAP_1", y = "UMAP_2") +
      theme(
        strip.text = element_text(size = 12),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid = element_blank()
      ) +
      guides(color = guide_colorbar(title = "Expression"))
    
    print(p)
  } # end gene loop
  
  dev.off()
} # end palette loop


### Plots for cell markers
# Celltype-clusters via marker genes
goi_celltype <- c("Lin28a", "Sall4", "Fzd7",                                    # undifferentiated cells; Fernandez-Guerrero et al., 2021
                  "Bmp2", "Col2a1", "Sox9", "Acan",                             # chondrogenic differentiation and limb skeletal, digit and joint morphogenesis; Fernandez-Guerrero et al., 2021
                  "Runx2", "Dlx5", "Sp7",                                       # osteoblast progenitor and differentiation marker 
                  "Alx3", "Alx4",                                               # proximal anterior mesoderm; Fernandez-Guerrero et al., 2021                                             # autopod-associated genes; Fernandez-Guerrero et al., 2021
                  "Asb4",                                                       # distal anterior mesoderm; Fernandez-Guerrero et al., 2021 
                  "Krt8", "Krt14", "Krt15",                                     # epidermis/ epidermal keratins; Kelly et al., 2020; Fernandez-Guerrero et al., 2021 
                  "Cdh5",                                                       # vasculature @E11.5; Kelly et al., 2020
                  "Lyz2",                                                       # blood @E11.5; Kelly et al., 2020
                  "Bcan", "Dsg2", "Esrp1", "Epcam",                             # epithelial markers; Fernandez-Guerrero et al., 2021 
                  "Pax3", "Pax7", "Myf5", "Myog",                               # muscle (progenitor) markers
                  "Twist1", "Ctcf"
)

### Plots for subset of tissues
library(ggplot2)
has_ggrastr <- requireNamespace("ggrastr", quietly = TRUE)

# define color palettes and suffixes
palettes <- list(
  red    = c("darksalmon", "red", "brown4"),
  brown  = c("beige", "darkgoldenrod3", "coral3"),
  turquoise   = c("aquamarine", "aquamarine4", "darkslategrey"),
  pink   = c("azure3", "deeppink1", "deeppink4")
)

# parameters
top_quantile <- 0.90    # top 10% of positives
nonexp_size <- 0.35     
low_size <- 0.7         
high_size <- 1.6        
grey_col <- "grey80"

# subset of interest
subset_levels <- c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5", 
                   "hindlimb_E10.5", "hindlimb_E11.5", "autopod_E12.5",
                   "hindlimb_E13.5", "autopod_E13.5", 
                   "hindlimb_E15.5", "hindlimb_E18.5")

# compute a PDF width that scales with number of samples (adjust multiplier as needed)
n_samples <- length(subset_levels)
pdf_width <- max(10, n_samples * 3.5)  # ~3.5 inches per sample; adjust to taste
pdf_height <- 5

for (pal_name in names(palettes)) {
  palette_cols <- palettes[[pal_name]]
  
  outFile <- file.path(
    output_folder,
    paste0("integrated_mandible_hindlimb.SUBSET.UMAP.goi_celltype_ContrastRaster_", pal_name, ".pdf")
  )
  
  pdf(outFile, width = pdf_width, height = pdf_height)
  
  for (gene in goi_celltype) {
    if (!(gene %in% rownames(so))) {
      message("Gene not found in data: ", gene)
      next
    }
    
    emb <- Embeddings(so, "umap.integrated")
    fetch <- FetchData(so, vars = c(gene, "orig.ident"))
    
    # build df, keep orig.ident as character first, then filter, then set factor with desired order
    df <- data.frame(
      cell = rownames(fetch),
      UMAP_1 = emb[rownames(fetch), 1],
      UMAP_2 = emb[rownames(fetch), 2],
      expr = as.numeric(fetch[, gene]),
      orig.ident = as.character(fetch[, "orig.ident"]),
      stringsAsFactors = FALSE
    )
    
    # keep only requested subset_levels
    df <- df[df$orig.ident %in% subset_levels, , drop = FALSE]
    
    # if no cells for this gene in subset, skip with a message
    if (nrow(df) == 0) {
      message("No cells for gene ", gene, " in chosen subset. Skipping.")
      next
    }
    
    # enforce factor order (this will keep facets in the order subset_levels,
    # and drop any levels from subset_levels that aren't present in df automatically)
    df$orig.ident <- factor(df$orig.ident, levels = subset_levels)
    
    # split non / low / high (recomputed per gene)
    non_idx <- which(df$expr == 0 | is.na(df$expr))
    pos_idx <- which(df$expr > 0)
    
    if (length(pos_idx) == 0) {
      qtop <- NA
    } else {
      qtop <- as.numeric(quantile(df$expr[pos_idx], probs = top_quantile, na.rm = TRUE))
    }
    
    if (is.na(qtop) || qtop == 0) {
      low_idx <- integer(0)
      high_idx <- pos_idx
    } else {
      low_idx <- which(df$expr > 0 & df$expr <= qtop)
      high_idx <- which(df$expr > qtop)
      if (length(high_idx) == 0 && length(pos_idx) >= 1) {
        qtop2 <- as.numeric(quantile(df$expr[pos_idx], probs = 0.75, na.rm = TRUE))
        low_idx <- which(df$expr > 0 & df$expr <= qtop2)
        high_idx <- which(df$expr > qtop2)
      }
    }
    
    # always produce these df pieces (may be zero-row data.frames)
    df_non  <- df[non_idx, , drop = FALSE]
    df_low  <- df[low_idx, , drop = FALSE]
    df_high <- df[high_idx, , drop = FALSE]
    
    # build plot (use df as base to ensure facet variable is present somewhere)
    p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
      # Non-expressors (grey)
      (if (nrow(df_non) > 0) {
        if (has_ggrastr) ggrastr::geom_point_rast(
          data = df_non, aes(x = UMAP_1, y = UMAP_2),
          color = grey_col, size = nonexp_size, raster.dpi = 150
        ) else geom_point(
          data = df_non, aes(x = UMAP_1, y = UMAP_2),
          color = grey_col, size = nonexp_size
        )
      } else NULL) +
      # Low expressors (color by expr)
      (if (nrow(df_low) > 0) {
        if (has_ggrastr) ggrastr::geom_point_rast(
          data = df_low, aes(x = UMAP_1, y = UMAP_2, color = expr),
          size = low_size, raster.dpi = 150
        ) else geom_point(
          data = df_low, aes(x = UMAP_1, y = UMAP_2, color = expr),
          size = low_size
        )
      } else NULL) +
      # High expressors (color by expr, plotted last to be on top)
      (if (nrow(df_high) > 0) {
        if (has_ggrastr) ggrastr::geom_point_rast(
          data = df_high, aes(x = UMAP_1, y = UMAP_2, color = expr),
          size = high_size, raster.dpi = 150
        ) else geom_point(
          data = df_high, aes(x = UMAP_1, y = UMAP_2, color = expr),
          size = high_size
        )
      } else NULL) +
      facet_wrap(~ orig.ident, nrow = 1, drop = TRUE) +
      scale_color_gradientn(colors = palette_cols, na.value = palette_cols[1]) +
      ggtitle(gene) +
      theme_minimal() +
      labs(x = "UMAP_1", y = "UMAP_2") +
      theme(
        strip.text = element_text(size = 12),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid = element_blank()
      ) +
      guides(color = guide_colorbar(title = "Expression"))
    
    print(p)
  } # end gene loop
  
  dev.off()
} # end palette loop
