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

