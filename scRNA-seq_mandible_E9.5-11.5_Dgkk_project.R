################################################################################
# Individual and integrated analysis of scRNA-seq mandible E9.5, E10.5 and E11.5

# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)

################################################################################
################### Pre-analysis of mandible_E9.5 dataset  #####################

# Set out folder (to store results)
output_folder <- "~/Documents/postdoc/collaboration/Pauline/pre-analysis_mandibleE9.5/"

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
saveRDS(so_mandible_E9.5, file = "~/Documents/postdoc/collaboration/Pauline/pre-analysis_mandibleE9.5/so_mandible_E9.5.rds")


################################################################################
################## Pre-analysis of mandible_E10.5 dataset  #####################

# Set out folder (to store results)
output_folder <- "~/Documents/postdoc/collaboration/Pauline/pre-analysis_mandibleE10.5/"

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
             pt.size = 0.1,
             label = TRUE)
out_path <- paste(output_folder, "/mandible_E10.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 7, height = 5)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_mandible_E10.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/mandible_E10.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 7, height = 5)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_mandible_E10.5, file = "~/Documents/postdoc/collaboration/Pauline/pre-analysis_mandibleE10.5/so_mandible_E10.5.rds")


###################################  Plots  ####################################

# Define genes of interest (goi)
goi <- c("Dgkk", "Gpr50", "Hand2", "Dlx6")


################################ Violin plots  #################################
# Violin Plot for goi
outFile <- paste(output_folder,
                 "/mandible_E10.5.Violin.goi.seurat_clusters.pdf", 
                 sep = "")
pdf(outFile, width = 20, height = 5)
for (gene in goi) {
  if (gene %in% rownames(so_mandible_E10.5)) {
    p <- VlnPlot(so_mandible_E10.5,
                 features = gene,
                 group.by = "seurat_clusters",
                 pt.size = 0.1
    ) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(p)
  } else {
    message(paste("Gene not found in data:", gene))
  }
}
dev.off()


### UMAP Plots for goi
outFile <- paste0(output_folder, "/mandible_E10.5.UMAP.goi.orig.ident.pdf")
pdf(outFile, width = 7, height = 5)

for (gene in goi) {
  if (gene %in% rownames(so_mandible_E10.5)) {
    p <- FeaturePlot(
      so_mandible_E10.5, 
      features = gene,
      reduction = "umap",
      split.by = "orig.ident",
      pt.size = 0.8,
      order = TRUE,                              # put expressors on top
      min.cutoff = "q05",                        # shrink influence of 0’s
      cols = c("grey90", "purple")                  # light grey → red gradient
    )
    print(p)
  } else {
    message(paste("Gene not found in data: ", gene))
  }
}

dev.off()


################################################################################
################### Pre-analysis of mandible_E11.5 dataset  ####################

# Set out folder (to store results)
output_folder <- "~/Documents/postdoc/collaboration/Pauline/pre-analysis_mandibleE11.5/"

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
             pt.size = 0.8,
             label = TRUE)
out_path <- paste(output_folder, "/mandible_E11.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 7, height = 5)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_mandible_E11.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/mandible_E11.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 7, height = 5)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_mandible_E11.5, file = "~/Documents/postdoc/collaboration/Pauline/pre-analysis_mandibleE11.5/so_mandible_E11.5.rds")


################################################################################
####################### INTEGRATION OF MANDIBLE DATASETS #######################

# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/collaboration/Pauline/mandible_E9.5-11.5_integration/"

# Load raw data (pre-processed Seurat objects produced above)
so_mandible_E9.5 <- readRDS("~/Documents/postdoc/collaboration/Pauline/pre-analysis_mandibleE9.5/so_mandible_E9.5.rds")
so_mandible_E10.5 <- readRDS("~/Documents/postdoc/collaboration/Pauline/pre-analysis_mandibleE10.5/so_mandible_E10.5.rds")
so_mandible_E11.5 <- readRDS("~/Documents/postdoc/collaboration/Pauline/pre-analysis_mandibleE11.5/so_mandible_E11.5.rds")


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
             group.by = c("orig.ident", "seurat_clusters"),
             pt.size = 0.8)
outFile <- paste(output_folder, "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.orig.ident.clusters.pdf", sep = "")
pdf(outFile, width = 10, height = 5)
plot(p)
dev.off()

p <- DimPlot(object = so_mandible_E9.5_E10.5_E11.5_integrated,
             reduction = "umap.integrated",
             group.by = 'seurat_clusters',
             split.by = 'orig.ident',
             label = TRUE,
             pt.size = 0.8)
outFile <- paste(output_folder, "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.orig.ident.clusters_split.pdf", sep = "")
pdf(outFile, width = 15, height = 5)
plot(p)
dev.off()


# Save the Seurat object
outFile <- paste(output_folder, "so_mandible_E9.5_E10.5_E11.5_integrated.rds", sep = "")
saveRDS(so_mandible_E9.5_E10.5_E11.5_integrated, file = outFile)



##############################################################################
################################# Plots ######################################

# Load integrated Dataset
so_mandible_E9.5_E10.5_E11.5_integrated <- readRDS("~/Documents/postdoc/collaboration/Pauline/mandible_E9.5-11.5_integration/so_mandible_E9.5_E10.5_E11.5_integrated.rds")

# Define genes of interest (goi)
goi <- c("Dgkk", "Dlx6", "Hand2")

# Define output folder for publication
output_folder <- "~/Documents/postdoc/collaboration/Pauline/Dgkk_project/results_for_publication/"

########################### Violin plots per dataset ###########################
# Violin Plot for goi
outFile <- paste(output_folder,
                 "/mandible_E9.5_E10.5_E11.5_integrated.Violin.goi.seurat_clusters.pdf", 
                 sep = "")
pdf(outFile, width = 20, height = 5)
for (gene in goi) {
  if (gene %in% rownames(so_mandible_E9.5_E10.5_E11.5_integrated)) {
    p <- VlnPlot(so_mandible_E9.5_E10.5_E11.5_integrated,
                 features = gene,
                 group.by = "seurat_clusters",
                 pt.size = 0.1
    ) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(p)
  } else {
    message(paste("Gene not found in data:", gene))
  }
}
dev.off()

################## Violin plots per dataset and seurat cluster #################
so <- so_mandible_E9.5_E10.5_E11.5_integrated

# Ensure clusters are identities
Idents(so) <- "seurat_clusters"

# Create composite cluster + sample identity
so$cluster_identity <- paste0("Cluster_", Idents(so), "_", so$orig.ident)

# Get all unique combinations, sorted by desired sample order and cluster
sample_order <- c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5")

# Build desired factor levels
all_combos <- expand.grid(
  cluster = sort(unique(Idents(so))),
  sample = sample_order
)
cluster_levels <- paste0("Cluster_", all_combos$cluster, "_", all_combos$sample)

# Set the cluster_identity as a factor with desired order
so$cluster_identity <- factor(so$cluster_identity, levels = cluster_levels)

outFile <- paste(output_folder,
                 "/mandible_E9.5_E10.5_E11.5_integrated.Violin.goi.by.cluster.and.orig.ident.pdf", 
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


#################################### DOTPLOTS ##################################
################################################################################

# Exploratory genes of interest (goi) according to (Xu et al., 2019)
goi <- c("Dgkk", 
         "Dlx5", "Dlx6", 
         "Hand1",                         # (Xu et al., 2019) 
         "Hand2",                         # distal domain marker gene (Xu et al., 2019, NC1/2 Subgroups 0 and 1)
         "Pou3f3", "Emx2", "Barx1",       # proximal domain marker genes (Xu et al., 2019, NC1/2 Subgroup 2)
         "Lhx6", "Lhx8",                  # rostral domain marker genes (Xu et al., 2019, NC1/2 Subgroup 4)
         "Gsc",                           # caudal domain marker gene  (Xu et al., 2019, NC1/2 Subgroup 3)
         "Epcam",                         # Epithelial cells (Xu et al., 2019)
         "Pecam1",                        # Endothelial cells (Xu et al., 2019)
         "Msx1", "Msx2", "Bambi"          # BMP signaling targets (Xu et al., 2019, NC1/2 Subgroup)
         )

# Exploratory genes of interest (goi)  (Yuan et al., 2020)
goi <- c("Dgkk", "Hand2", "Dlx6",             # Genes of interest for this study
         "Hand1", "Alx1", "Alx3",             # distal domain markers (Yuan et al., 2020)
         "Gsc",                               # aboral domain marker (Yuan et al., 2020)
         "Foxc1", "Dclk1",                    # proximo-aboral (Yuan et al., 2020 cluster 1), marker cluster 2+4 in our analysis
         "Barx1", "Nbl1",                     # proximal domain marker gene Yuan et al., 2020), broadly expressed clusters 2-6 in our analysis
         "Pitx1",                             # general mesenchyme (Yuan et al., 2020)
         "Notch2",                            #  (Yuan et al., 2020)
         "Lhx6", "Lhx8", "Etv4",              # oral domain marker genes (Yuan et al., 2020)
         "Cited1", "Mecom", "Crlf1",          #  (Yuan et al., 2020)
         "Epcam",                        # Epithelial marker
         # What is cluster 8 -> marker gene!
         "Car2",                         #  (Yuan et al., 2020), marker cluster 9 in our analysis
         "Cgnl1",                        #  (Yuan et al., 2020), marker cluster 10 in our analysis
         "Maf"                           #  (Yuan et al., 2020), marker cluster 11 in our analysis
         )

# Define genes of interest (goi)  (Xu et al., 2019; Yuan et al., 2020)
goi <- c("Dgkk", "Hand2", "Dlx6",        # Genes of interest for this study
         "Hand1", "Alx1", "Alx3",        # distal domain markers (Yuan et al., 2020)
         "Gsc",                          # aboral domain marker (Yuan et al., 2020)
         "Foxc1", "Dclk1",               # proximo-aboral (Yuan et al., 2020 cluster 1), marker cluster 2+4 in our analysis
         "Barx1", "Nbl1",                # proximal domain marker gene Yuan et al., 2020), broadly expressed clusters 2-6 in our analysis
         "Pitx1",                        # general mesenchyme (Yuan et al., 2020)
         "Notch2",                       # (Yuan et al., 2020)
         "Lhx6", "Lhx8", "Etv4",         # oral domain marker genes (Yuan et al., 2020)
         "Cited1", "Mecom", "Crlf1",     # (Yuan et al., 2020)
         "Epcam",                        # Epithelial marker
         "Malat1", "Xist",               # Dataset artefact, most likely, from incomplete lysis of a subset of single cells, likely "The NC3 cluster, consisting of 484 cells, is distinguished from the other clusters primarily by overall lower levels of transcripts and dramatically lower number of genes detected per cell [...]. In addition, the NC3 cluster cells showed very low levels or absence of several ubiquitously expressed, nuclear localized long noncoding RNAs, including Malat1, Xist, Meg3, and Kcnq1ot1" (Xu et al., 2019)
         "Car2",                         # (Yuan et al., 2020), marker cluster 9 in our analysis
         "Cgnl1",                        # (Yuan et al., 2020), marker cluster 10 in our analysis
         "Maf"                           # (Yuan et al., 2020), marker cluster 11 in our analysis
)

# Define genes of interest (goi), final list for revision (Xu et al., 2019; Yuan et al., 2020; Paulines candidates)
goi <- c("Dgkk", "Hand1", "Hand2", "Alx1", "Alx3",        # Genes of interest for this study + distal domain markers (Yuan et al., 2020)
         "Foxf1", "Lhx8", "Etv4",               # (distal-)oral markers (Pauline; Yuan et al., 2020)
         "Cgnl1", "Maf",                                  # aboral (Pauline; Yuan et al., 2020), markers cluster 10/11 in our analysis
         "Foxc1",                                         # proximo-aboral (Yuan et al., 2020 cluster 1), marker cluster 2+4 in our analysis
         "Nbl1", "Pou3f3", "Car2", "Crlf1",               # proximal domain marker genes (Yuan et al., 2020; Pauline); marker cluster 9 in our analysis
         "Barx1", "Pitx1", "Dlx5", "Dlx6", "Notch2",               # ectomesenchyme (Pauline; Yuan et al., 2020)
         "Epcam",                                         # Epithelial marker
         "Bmp4", "Edn1"                                  # markers to address reviewer comment (Pauline))
)

#################### DotPlot: GOI to define clusters (revision) ###################
# ================================================
# FILTER GENES PRESENT IN OBJECT
# ================================================
goi_present <- goi[goi %in% rownames(so_mandible_E9.5_E10.5_E11.5_integrated)]

missing_genes <- setdiff(goi, goi_present)
if (length(missing_genes) > 0) {
  message("Missing genes: ", paste(missing_genes, collapse = ", "))
}

# ================================================
# DOTPLOT
# ================================================
outFile <- paste(
  output_folder,
  "/mandible_E9.5_E10.5_E11.5_integrated.DotPlot.goi.seurat_clusters.pdf",
  sep = ""
)

pdf(outFile, width = 10, height = 5)

p <- DotPlot(
  so_mandible_E9.5_E10.5_E11.5_integrated,
  features = goi_present,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  ggtitle("Gene expression (DotPlot)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p)

dev.off()


# Define genes of interest (goi) for Dgk genes to assess expression across clusters
goi <- c("Dgka", "Dgkb", "Dgkd", "Dgke", "Dgkeos", "Dgkg", "Dgkh", "Dgki", 
         "Dgkk", "Dgkq", "Dgkz")

#################### DotPlot: GOI to define clusters (revision) ###################
# ================================================
# FILTER GENES PRESENT IN OBJECT
# ================================================
goi_present <- goi[goi %in% rownames(so_mandible_E9.5_E10.5_E11.5_integrated)]

missing_genes <- setdiff(goi, goi_present)
if (length(missing_genes) > 0) {
  message("Missing genes: ", paste(missing_genes, collapse = ", "))
}

# ================================================
# DOTPLOT
# ================================================
outFile <- paste(
  output_folder,
  "/mandible_E9.5_E10.5_E11.5_integrated.DotPlot.Dgk-genes.seurat_clusters.pdf",
  sep = ""
)

pdf(outFile, width = 10, height = 5)

p <- DotPlot(
  so_mandible_E9.5_E10.5_E11.5_integrated,
  features = goi_present,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  ggtitle("Gene expression (DotPlot)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p)

dev.off()

#################### Heatmap: Dgk genes across clusters ###################

# ================================================
# FILTER GENES PRESENT IN OBJECT
# ================================================
goi_present <- goi[goi %in% rownames(so_mandible_E9.5_E10.5_E11.5_integrated)]

missing_genes <- setdiff(goi, goi_present)
if (length(missing_genes) > 0) {
  message("Missing genes: ", paste(missing_genes, collapse = ", "))
}

# ================================================
# OPTIONAL: SET CLUSTER ORDER
# ================================================
Idents(so_mandible_E9.5_E10.5_E11.5_integrated) <- "seurat_clusters"

# Example:
# so_mandible_E9.5_E10.5_E11.5_integrated$seurat_clusters <- factor(
#   so_mandible_E9.5_E10.5_E11.5_integrated$seurat_clusters,
#   levels = c("0","1","2","3","4","5")
# )

# ================================================
# HEATMAP
# ================================================
outFile <- paste(
  output_folder,
  "/mandible_E9.5_E10.5_E11.5_integrated.Heatmap.Dgk-genes.seurat_clusters.pdf",
  sep = ""
)

pdf(outFile, width = 10, height = 8)

p <- DoHeatmap(
  object = so_mandible_E9.5_E10.5_E11.5_integrated,
  features = goi_present,
  group.by = "seurat_clusters",
  slot = "scale.data"
) +
  ggtitle("Dgk gene expression across clusters") +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p)

dev.off()

#################### DotPlot: all GOI, timepoints top → bottom ###################

library(Seurat)
library(viridis)

# Seurat object
so <- so_mandible_E9.5_E10.5_E11.5_integrated

# Ensure clusters are identities
Idents(so) <- "seurat_clusters"

# Define timepoint order
sample_order <- c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5")

# Create composite identity: Cluster_X + sample
so$cluster_identity <- paste0(
  "Cluster_",
  Idents(so),
  "_",
  so$orig.ident
)

# Get unique clusters
clusters <- sort(as.numeric(levels(Idents(so))))

# Cluster-first ordering:
# Cluster1_E9.5, Cluster2_E9.5, ..., Cluster1_E10.5, ...
cluster_levels <- unlist(
  lapply(sample_order, function(sample) {
    paste0("Cluster_", clusters, "_", sample)
  })
)

# IMPORTANT: reverse levels so Cluster 1 appears at the TOP
so$cluster_identity <- factor(
  so$cluster_identity,
  levels = rev(cluster_levels)
)

# Keep only genes present in the object
valid_goi <- goi[goi %in% rownames(so)]
if (length(valid_goi) == 0) {
  stop("None of the genes in 'goi' were found in the Seurat object.")
}

# Create DotPlot
p <- DotPlot(
  so,
  features = valid_goi,
  group.by = "cluster_identity",
  scale = TRUE
) +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90")
  ) +
  labs(
    x = "Gene",
    y = "Cluster × timepoint",
    color = "Avg. expression (scaled)",
    size = "% expressing"
  )

# Save to PDF
outFile <- paste(
  output_folder,
  "/mandible_E9.5_E10.5_E11.5_integrated.DotPlot.goi_all.timepoint_top_to_bottom.pdf",      # mandible_E9.5_E10.5_E11.5_integrated.DotPlot.Dgk-genes.timepoint_top_to_bottom.pdf"
  sep = ""
)

pdf(
  outFile,
  width = 12,
  height = 0.5 * length(valid_goi) + 2
)

print(p)
dev.off()



############### DotPlot: all GOI, ordered by cluster top → bottom ##############
library(Seurat)
library(viridis)

# Seurat object
so <- so_mandible_E9.5_E10.5_E11.5_integrated

# Ensure cluster identities are set
Idents(so) <- "seurat_clusters"

# Define timepoint order
sample_order <- c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5")

# Create composite cluster × timepoint identity
so$cluster_identity <- paste0(
  "Cluster_",
  Idents(so),
  "_",
  so$orig.ident
)

# Build desired ordering:
# Cluster 1 (E9.5 → E10.5 → E11.5), Cluster 2, ...
clusters <- sort(as.numeric(levels(Idents(so))))

cluster_levels <- unlist(
  lapply(clusters, function(cl) {
    paste0("Cluster_", cl, "_", sample_order)
  })
)

# IMPORTANT: reverse levels so Cluster 1 appears at the TOP
so$cluster_identity <- factor(
  so$cluster_identity,
  levels = rev(cluster_levels)
)

# Keep only genes present in the object
valid_goi <- goi[goi %in% rownames(so)]
if (length(valid_goi) == 0) {
  stop("None of the genes in 'goi' were found in the Seurat object.")
}

# Create DotPlot
p <- DotPlot(
  so,
  features = valid_goi,
  group.by = "cluster_identity",
  scale = TRUE
) +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90")
  ) +
  labs(
    x = "Gene",
    y = "Cluster × timepoint",
    color = "Avg. expression (scaled)",
    size = "% expressing"
  )

# Save to PDF
outFile <- paste(
  output_folder,
  "/mandible_E9.5_E10.5_E11.5_integrated.DotPlot.all_goi.cluster_top_to_bottom.pdf",
  sep = ""
)

pdf(
  outFile,
  width = 12,
  height = 0.5 * length(valid_goi) + 2
)

print(p)
dev.off()


################ DotPlot: all GOI, ordered by cluster left → right #############

library(Seurat)
library(viridis)

# Seurat object
so <- so_mandible_E9.5_E10.5_E11.5_integrated

# Ensure clusters are identities
Idents(so) <- "seurat_clusters"

# Define timepoint order
sample_order <- c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5")

# Create composite identity: Cluster_X + sample
so$cluster_identity <- paste0(
  "Cluster_",
  Idents(so),
  "_",
  so$orig.ident
)

# Get unique clusters (numeric order)
clusters <- sort(as.numeric(levels(Idents(so))))

# CLUSTER-FIRST ordering:
# Cluster1_E9.5, Cluster1_E10.5, Cluster1_E11.5,
# Cluster2_E9.5, Cluster2_E10.5, Cluster2_E11.5, ...
cluster_levels <- unlist(
  lapply(clusters, function(cl) {
    paste0("Cluster_", cl, "_", sample_order)
  })
)

# Apply ordering (Cluster 1 block will be far LEFT after coord_flip())
so$cluster_identity <- factor(
  so$cluster_identity,
  levels = cluster_levels
)

# Keep only genes present
valid_goi <- goi[goi %in% rownames(so)]
if (length(valid_goi) == 0) {
  stop("None of the genes in 'goi' were found in the Seurat object.")
}

# Create DotPlot
p <- DotPlot(
  so,
  features = valid_goi,
  group.by = "cluster_identity",
  scale = TRUE
) +
  coord_flip() +
  scale_color_viridis_c() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 9
    ),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90")
  ) +
  labs(
    x = "Gene",
    y = "Cluster × timepoint",
    color = "Avg. expression (scaled)",
    size = "% expressing"
  )

# Save to PDF
outFile <- paste(
  output_folder,
  "/mandible_E9.5_E10.5_E11.5_integrated.DotPlot.all_goi.cluster_left_to_right.pdf",
  sep = ""
)

pdf(
  outFile,
  width = 0.5 * length(valid_goi) + 2,
  height = 12
)

print(p)
dev.off()


############################ UMAP plots for goi ################################

# Set the desired order of orig.ident
so_mandible_E9.5_E10.5_E11.5_integrated$orig.ident <- factor(
  so_mandible_E9.5_E10.5_E11.5_integrated$orig.ident,
  levels = c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5")
)

# Visualize clusters of integrated dataset as Dimplot
p <- DimPlot(object = so_mandible_E9.5_E10.5_E11.5_integrated,
             reduction = "umap.integrated",
             group.by = 'seurat_clusters',
             split.by = "orig.ident",
             pt.size = 0.8,
             label = TRUE)
out_path <- paste(output_folder, "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.Dimplot.orig.ident.pdf", sep = "")
pdf(out_path, width = 30, height = 10)
plot(p)
dev.off()

# Visualize clusters of integrated dataset as Dimplot
p <- DimPlot(object = so_mandible_E9.5_E10.5_E11.5_integrated,
             reduction = "umap.integrated",
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()


### UMAP Plots for goi
outFile <- paste(output_folder,
                 "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.goi.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 30, height = 10)
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

### UMAP Plots for goi
outFile <- paste(output_folder,
                 "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.goi.pdf", 
                 sep = "")
pdf(outFile, width = 15, height = 10)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi) {
  if (gene %in% rownames(so_mandible_E9.5_E10.5_E11.5_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_mandible_E9.5_E10.5_E11.5_integrated, 
                     features = gene,
                     reduction = "umap.integrated")
    plot(p)
  } else {
    # Print a message for missing genes (optional)
    message(paste("Gene not found in data: ", gene))
  }
}

dev.off()


### Combinatory plots
# Hand2 + Dgkk
p1 <- FeaturePlot(
  so_mandible_E9.5_E10.5_E11.5_integrated,
  features = c("Hand2", "Dgkk"),
  reduction = "umap.integrated",
  blend = TRUE,                  # blend colors
  cols = c("grey90", "blue", "magenta"),
  order = TRUE                   # plot higher expressers on top
)

# Dlx6 + Dgkk
p2 <- FeaturePlot(
  so_mandible_E9.5_E10.5_E11.5_integrated,
  features = c("Dlx6", "Dgkk"),
  reduction = "umap.integrated",
  blend = TRUE,
  cols = c("grey90", "blue", "magenta"),
  order = TRUE
)

# Save to PDF
outFile <- paste0(output_folder, "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.Hand2_Dlx6_Dgkk.pdf")
pdf(outFile, width = 25, height = 7)
print(p1)
print(p2)
dev.off()


### Overlaps at each timepoint separated
# Define pairs of interest
gene_pairs <- list(
  c("Hand2", "Dgkk"),
  c("Dlx6", "Dgkk")
)

# Output file
outFile <- paste0(output_folder,
                  "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.Hand2_Dlx6_Dgkk.orig.ident.pdf")
pdf(outFile, width = 20, height = 10)

for (pair in gene_pairs) {
  if (all(pair %in% rownames(so_mandible_E9.5_E10.5_E11.5_integrated))) {
    p <- FeaturePlot(
      so_mandible_E9.5_E10.5_E11.5_integrated,
      features = pair,
      reduction = "umap.integrated",
      blend = TRUE,
      cols = c("grey90", "blue", "magenta"),
      order = TRUE,
      split.by = "orig.ident"
    )
    print(p)
  } else {
    message(paste("One or both genes not found in data: ", paste(pair, collapse = ", ")))
  }
}

dev.off()


## Downsized UMAP plots for goi, increased dot size and plotting +cells last ###

# Set the desired order of orig.ident
so_mandible_E9.5_E10.5_E11.5_integrated$orig.ident <- factor(
  so_mandible_E9.5_E10.5_E11.5_integrated$orig.ident,
  levels = c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5")
)

# 1. DimPlot split by orig.ident
out_path <- paste0(output_folder, "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.Dimplot.orig.ident_downsize.pdf")
pdf(out_path, width = 15, height = 5)
p <- DimPlot(
  so_mandible_E9.5_E10.5_E11.5_integrated,
  reduction = "umap.integrated",
  group.by = "seurat_clusters",
  split.by = "orig.ident",
  label = TRUE,
  pt.size = 3,
  raster = TRUE
)
print(p)
dev.off()

# 2. DimPlot clusters (all together)
out_path <- paste0(output_folder, "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.Dimplot.clusters_downsize.pdf")
pdf(out_path, width = 10, height = 7)
p <- DimPlot(
  so_mandible_E9.5_E10.5_E11.5_integrated,
  reduction = "umap.integrated",
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 2,
  raster = TRUE
)
print(p)
dev.off()

# 3.0 FeaturePlots for goi split by orig.ident (positives on top)
outFile <- paste0(output_folder, "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.goi.orig.ident_downsize.pdf")
pdf(outFile, width = 15, height = 5)

for (gene in goi) {
  if (gene %in% rownames(so_mandible_E9.5_E10.5_E11.5_integrated)) {
    p <- FeaturePlot(
      so_mandible_E9.5_E10.5_E11.5_integrated,
      features = gene,
      reduction = "umap.integrated",
      split.by = "orig.ident",
      pt.size = 3,
      raster = TRUE,
      order = TRUE,                            # bring positive cells forward
      min.cutoff = 0.3,                      # reduce noise from low values
      cols = c("grey90", "purple")                # darker grey for negatives
    )
    print(p)
  } else {
    message(paste("Gene not found in data: ", gene))
  }
}

dev.off()

# 3.1 FeaturePlots for goi split by orig.ident with contrast raster + axes + improved colors
library(ggplot2)

# Desired order
orig_levels <- c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5")

so <- so_mandible_E9.5_E10.5_E11.5_integrated  # convenience short name

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

for (pal_name in names(palettes)) {
  palette_cols <- palettes[[pal_name]]
  
  outFile <- paste0(
    output_folder, "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.goi.orig.ident_ContrastRaster_", 
    pal_name, "_downsize.pdf"
  )
  
  pdf(outFile, width = 5, height = 12)
  
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
      facet_wrap(~ orig.ident, ncol = 1, scales = "fixed") +
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





# 4.0 FeaturePlots for goi (all integrated datasets together)
outFile <- paste0(output_folder, "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.goi_downsize.pdf")
pdf(outFile, width = 10, height = 7)
for (gene in goi) {
  if (gene %in% rownames(so_mandible_E9.5_E10.5_E11.5_integrated)) {
    p <- FeaturePlot(
      so_mandible_E9.5_E10.5_E11.5_integrated,
      features = gene,
      reduction = "umap.integrated",
      pt.size = 2,
      raster = TRUE,
      order = TRUE,
      min.cutoff = 0.05,                      # reduce noise from low values
    )
    print(p)
  } else {
    message(paste("Gene not found in data: ", gene))
  }
}
dev.off()

# 4.1 FeaturePlots for goi (all integrated datasets together), with contrast raster + axes + improved colors
# Adapted FeaturePlots with contrast raster style
library(ggplot2)

so <- so_mandible_E9.5_E10.5_E11.5_integrated  # convenience short name
orig_levels <- c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5")

outFile <- paste0(output_folder, "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.goi_ContrastRaster_downsize.pdf")

# parameters
top_quantile <- 0.90
nonexp_size <- 0.35
low_size <- 0.7
high_size <- 1.6
grey_col <- "grey80"
palette_cols <- c("lightblue", "steelblue", "purple")  # gradient for expression

pdf(outFile, width = 7, height = 5)

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
        size = low_size, raster.dpi = 300
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
    scale_color_gradientn(colors = palette_cols, na.value = palette_cols[1]) +
    ggtitle(gene) +
    theme_minimal() +
    labs(x = "umapintegrated_1", y = "umapintegrated_2") +
    theme(
      strip.text = element_text(size = 20),             # facet titles twice as big
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.line = element_line(color = "black", linewidth = 0.6),
      panel.grid = element_blank()
    ) +
    guides(color = guide_colorbar(title = "Expression"))
  
  print(p)
}

dev.off()


# 5. Combinatory plots Hand2+Dgkk / Dlx6+Dgkk
outFile <- paste0(output_folder, "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.Hand2_Dlx6_Dgkk_downsize.pdf")
pdf(outFile, width = 20, height = 5)
p1 <- FeaturePlot(
  so_mandible_E9.5_E10.5_E11.5_integrated,
  features = c("Hand2", "Dgkk"),
  reduction = "umap.integrated",
  blend = TRUE,
  cols = c("grey90", "blue", "magenta"),
  order = TRUE,
  pt.size = 3,
  raster = TRUE,
  #min.cutoff = 0.05                      # reduce noise from low values
)
p2 <- FeaturePlot(
  so_mandible_E9.5_E10.5_E11.5_integrated,
  features = c("Dlx6", "Dgkk"),
  reduction = "umap.integrated",
  blend = TRUE,
  cols = c("grey90", "blue", "magenta"),
  order = TRUE,
  pt.size = 3,
  raster = TRUE,
  min.cutoff = 0.001                      # reduce noise from low values
)
print(p1)
print(p2)
dev.off()


# 6. Overlaps at each timepoint
gene_pairs <- list(
  c("Hand2", "Dgkk"),
  c("Dlx6", "Dgkk")
)

outFile <- paste0(output_folder, "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.Hand2_Dlx6_Dgkk.orig.ident_downsize.pdf")
pdf(outFile, width = 15, height = 10)
for (pair in gene_pairs) {
  if (all(pair %in% rownames(so_mandible_E9.5_E10.5_E11.5_integrated))) {
    p <- FeaturePlot(
      so_mandible_E9.5_E10.5_E11.5_integrated,
      features = pair,
      reduction = "umap.integrated",
      blend = TRUE,
      cols = c("grey90", "blue", "magenta"),
      order = TRUE,
      split.by = "orig.ident",
      pt.size = 3,
      raster = TRUE
    )
    print(p)
  } else {
    message(paste("One or both genes not found in data: ", paste(pair, collapse = ", ")))
  }
}
dev.off()


################################################################################
############################### MARKER ANALYSIS  ###############################

# Change the default assay to "SCT"
DefaultAssay(so_mandible_E9.5_E10.5_E11.5_integrated) <- "SCT"

## Identify top 100 markers + all per cluster
# Correcting SCT counts before running FindAllMarkers
so_mandible_E9.5_E10.5_E11.5_integrated <- PrepSCTFindMarkers(so_mandible_E9.5_E10.5_E11.5_integrated,
                                                              assay = "SCT", 
                                                              verbose = TRUE)
markers <- FindAllMarkers(so_mandible_E9.5_E10.5_E11.5_integrated,
                          min.pct = 0.1,
                          test.use = "wilcox")

# View markers for all clusters
marker_table <- table(markers$cluster)  # Shows the number of markers for each cluster

# Show markers for the first few clusters
cluster_ids <- unique(so_mandible_E9.5_E10.5_E11.5_integrated$seurat_clusters)   # Get unique cluster identities
num_clusters <- length(cluster_ids)   # Count the number of unique clusters

# Loop over each cluster to extract and print top 25 markers
for (cluster in cluster_ids) {
  # Extract the top 25 markers for this cluster
  top_markers <- head(markers[markers$cluster == cluster, ], 100)  # Get top 100 markers for the current cluster
  # Print the top 20 markers for the current cluster
  cat("Top 100 markers for cluster ", cluster, " are: \n", sep = "")
  # Print the gene names (marker genes) for the current cluster
  print(top_markers$gene)   # Assuming 'gene' is the column containing marker gene names
  cat("\n")  # Add a line break between clusters
}

# Save the marker list to a CSV file
write.csv(markers, file = "~/Documents/postdoc/collaboration/Pauline/Dgkk_project/results_for_publication/mandible_E9.5_E10.5_E11.5_integrated_markers_by_cluster.csv", row.names = TRUE)


#################### GO term analysis using clusterProfiler ####################

# Get unique cluster IDs
cluster_ids <- unique(so_mandible_E9.5_E10.5_E11.5_integrated$seurat_clusters)

# Initialize an empty list to store top 100 markers per cluster
top100_markers_list <- list()

# Loop over clusters
for (cluster in cluster_ids) {
  
  # Extract top 100 markers for the current cluster
  top_markers <- head(markers[markers$cluster == cluster, ], 100)
  
  # Store only the gene names in the list, using cluster ID as name
  top100_markers_list[[as.character(cluster)]] <- top_markers$gene
  
  # Optional: print them for verification
  cat("Top 100 markers for cluster", cluster, ":\n")
  print(top_markers$gene)
  cat("\n")
}

# Use clusterProfiler for GO term analysis/cluster definition
library(clusterProfiler)
library(org.Mm.eg.db)

for (cluster in names(top100_markers_list)) {
  genes <- top100_markers_list[[cluster]]
  
  # Convert SYMBOL to ENTREZID
  gene_df <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  entrez_genes <- gene_df$ENTREZID
  
  # GO enrichment
  ego <- enrichGO(
    gene = entrez_genes,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  # Optional: save results to a file per cluster
  write.csv(as.data.frame(ego), file = paste0("~/Documents/postdoc/collaboration/Pauline/Gpr50_project/results/GO_cluster_", cluster, ".csv"), row.names = FALSE)
}


################################################################################
# Within-stage differential expression analysis, comparing Dgkk⁺ vs Dgkk⁻ cells 
# separately for E9.5, E10.5, and E11.5

# Prepare the Seurat object
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Make new Seurat object for plots
so <- so_mandible_E9.5_E10.5_E11.5_integrated

DefaultAssay(so) <- "SCT"
so <- PrepSCTFindMarkers(so, verbose = TRUE)


# Define Dgkk positivity
so$Dgkk_status <- ifelse(
  GetAssayData(so, layer = "data")["Dgkk", ] > 0,
  "Dgkk_pos",
  "Dgkk_neg"
)

so$Dgkk_status <- factor(
  so$Dgkk_status,
  levels = c("Dgkk_neg", "Dgkk_pos")
)


# Differential expression per developmental stage
run_Dgkk_DE <- function(seurat_obj, stage) {
  
  message("Running DE for ", stage)
  
  so_sub <- subset(
    seurat_obj,
    subset = orig.ident == stage
  )
  
  Idents(so_sub) <- "Dgkk_status"
  
  markers <- FindMarkers(
    so_sub,
    ident.1 = "Dgkk_pos",
    ident.2 = "Dgkk_neg",
    test.use = "wilcox",
    min.pct = 0.1,
    logfc.threshold = 0.25
  )
  
  markers$gene <- rownames(markers)
  markers$stage <- stage
  
  return(markers)
}


# Run for all stages
stages <- c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5")

Dgkk_markers <- lapply(stages, function(s) run_Dgkk_DE(so, s))
Dgkk_markers <- bind_rows(Dgkk_markers)


# Save DE results (logFC + FDR)
write.csv(
  Dgkk_markers,
  file = paste0(
    output_folder,
    "/Dgkk_pos_vs_neg.DE.by_stage.csv"
  ),
  row.names = FALSE
)


### Heatmap of top enriched genes per stage
# Heatmap, stratified by Dgkk status and stage

# Ensure correct assay
DefaultAssay(so) <- "SCT"

# Ensure stage order (left → right)
so$orig.ident <- factor(
  so$orig.ident,
  levels = c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5")
)

# Sanity check: Dgkk_status must already exist
table(so$Dgkk_status)
# should show Dgkk_pos / Dgkk_neg

# Create composite grouping: stage × Dgkk status
so$stage_Dgkk <- paste0(
  so$orig.ident, "_", so$Dgkk_status
)

# Explicit ordering: stage first, then Dgkk− → Dgkk+
so$stage_Dgkk <- factor(
  so$stage_Dgkk,
  levels = c(
    "mandible_E9.5_Dgkk_neg", "mandible_E9.5_Dgkk_pos",
    "mandible_E10.5_Dgkk_neg", "mandible_E10.5_Dgkk_pos",
    "mandible_E11.5_Dgkk_neg", "mandible_E11.5_Dgkk_pos"
  )
)

# Select top enriched genes per stage (from your DE results)
top_genes <- Dgkk_markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  group_by(stage) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  pull(gene) %>%
  unique()

# Keep only genes present
top_genes <- top_genes[top_genes %in% rownames(so)]

# Plot heatmap
pdf(
  paste0(
    output_folder,
    "/Dgkk_top_enriched_genes.Heatmap.by_stage_and_Dgkk_status.pdf"
  ),
  width = 14,   # wider to avoid clipping
  height = 0.35 * length(top_genes) + 4
)

DoHeatmap(
  so,
  features = top_genes,
  group.by = "stage_Dgkk"
) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 10
    ),
    plot.margin = margin(
      t = 10,
      r = 80,   # 👈 critical: prevents last label cutoff
      b = 10,
      l = 10
    )
  ) +
  labs(
    title = "Top genes enriched in Dgkk⁺ cells\nstratified by developmental stage",
    fill = "Scaled\nexpression"   # legend title
  )

dev.off()


### DotPlot (stage × Dgkk status)
# Explicit ordering: stage first, then Dgkk− → Dgkk+
so$stage_Dgkk <- factor(
  so$stage_Dgkk,
  levels = c(
    "mandible_E9.5_Dgkk_neg", "mandible_E9.5_Dgkk_pos",
    "mandible_E10.5_Dgkk_neg", "mandible_E10.5_Dgkk_pos",
    "mandible_E11.5_Dgkk_neg", "mandible_E11.5_Dgkk_pos"
  )
)

Idents(so) <- "stage_Dgkk"

pdf(
  paste0(output_folder, "/Dgkk_top_enriched_genes.DotPlot.by_stage_and_Dgkk_status.pdf"),
  width = 10,
  height = 20
)

DotPlot(
  so,
  features = top_genes,
  scale = TRUE
) +
  coord_flip() +
  scale_color_viridis_c(option = "cividis") +
  labs(
    x = "Gene",
    y = "Stage × Dgkk status",
    color = "Avg. expression",
    size = "% expressing"
  ) +
  theme(axis.text.x = element_text(angle = 90))

dev.off()

########## Violin plots
# Sanity-check violin plots: Dgkk+ vs Dgkk−, split by stage
# Select top genes to visualize
genes_to_plot <- head(top_genes, 50)

outFile <- paste0(
  output_folder,
  "/Dgkk_sanitycheck.Violin.top50_genes.by_stage.pdf"
)

pdf(
  outFile,
  width  = 18,
  height = 1 * length(genes_to_plot) + 4
)

p <- VlnPlot(
  so,
  features = genes_to_plot,
  group.by = "Dgkk_status",
  split.by = "orig.ident",
  pt.size = 0,
  stack = TRUE,
  flip = TRUE
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10)
  ) +
  labs(
    x = "Scaled expression",
    y = "Gene",
    title = "Top Dgkk⁺-enriched genes (sanity check)"
  )

print(p)
dev.off()

## ===============================
## Violin plots per stage (E9.5 / E10.5 / E11.5)
## ===============================

library(Seurat)
library(dplyr)
library(ggplot2)

## -------------------------------
## Inputs you already have
## -------------------------------
# so_mandible_E9.5_E10.5_E11.5_integrated
# Dgkk_markers
# output_folder

## -------------------------------
## Load integrated object ONCE
## -------------------------------
so <- so_mandible_E9.5_E10.5_E11.5_integrated
DefaultAssay(so) <- "SCT"

## -------------------------------
## Ensure Dgkk_status exists
## (safe even if it already does)
## -------------------------------
if (!"Dgkk_status" %in% colnames(so@meta.data)) {
  message("Adding Dgkk_status metadata")
  
  so$Dgkk_status <- ifelse(
    GetAssayData(so, assay = "SCT", layer = "data")["Dgkk", ] > 0,
    "Dgkk_pos",
    "Dgkk_neg"
  )
  
  so$Dgkk_status <- factor(
    so$Dgkk_status,
    levels = c("Dgkk_neg", "Dgkk_pos")
  )
}

## Quick sanity check
print(table(so$Dgkk_status, useNA = "ifany"))

## -------------------------------
## Stages in chronological order
## -------------------------------
stages <- c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5")

## -------------------------------
## Open PDF
## -------------------------------
pdf(
  file   = paste0(output_folder,
                  "/Dgkk_sanitycheck.Violin.top50_genes.by_stage.pdf"),
  width  = 18,
  height = 60   # room for 50 stacked violins
)

## -------------------------------
## Loop over stages
## -------------------------------
for (stage_i in stages) {
  
  message("Plotting stage: ", stage_i)
  
  ## ---- subset Seurat object to this dataset ONLY
  so_stage <- subset(
    so,
    subset = orig.ident == stage_i
  )
  
  ## ---- select top 50 DE genes for this stage
  top_genes_stage <- Dgkk_markers %>%
    filter(
      stage == stage_i,
      p_val_adj < 0.05,
      avg_log2FC > 0.5
    ) %>%
    slice_max(order_by = avg_log2FC, n = 50) %>%
    pull(gene)
  
  ## ---- keep genes present in the object
  top_genes_stage <- intersect(
    top_genes_stage,
    rownames(so_stage)
  )
  
  if (length(top_genes_stage) == 0) {
    message("No genes to plot for ", stage_i)
    next
  }
  
  ## ---- violin plot
  p <- VlnPlot(
    so_stage,
    features = top_genes_stage,
    group.by = "Dgkk_status",
    pt.size  = 0,
    stack    = TRUE,
    flip     = TRUE
  ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9),
      plot.title  = element_text(face = "bold")
    ) +
    labs(
      x = "SCT-normalized expression",
      y = "Gene",
      title = paste0(
        "Top 50 enriched genes in Dgkk+ cells in ",
        stage_i
      )
    )
  
  print(p)
}

dev.off()


## ===============================
## Violin plots: top50 genes per stage
## shown across ALL stages
## ===============================

library(Seurat)
library(dplyr)
library(ggplot2)

## -------------------------------
## Input object
## -------------------------------
so <- so_mandible_E9.5_E10.5_E11.5_integrated
DefaultAssay(so) <- 'SCT'

## -------------------------------
## Ensure stage order
## -------------------------------
so$orig.ident <- factor(
  so$orig.ident,
  levels = c('mandible_E9.5', 'mandible_E10.5', 'mandible_E11.5')
)

## -------------------------------
## Ensure Dgkk_status exists
## -------------------------------
if (!'Dgkk_status' %in% colnames(so@meta.data)) {
  so$Dgkk_status <- ifelse(
    GetAssayData(so, assay = 'SCT', layer = 'data')['Dgkk', ] > 0,
    'Dgkk_pos',
    'Dgkk_neg'
  )
  so$Dgkk_status <- factor(
    so$Dgkk_status,
    levels = c('Dgkk_neg', 'Dgkk_pos')
  )
}

## -------------------------------
## Stages in chronological order
## -------------------------------
stages <- c('mandible_E9.5', 'mandible_E10.5', 'mandible_E11.5')

## -------------------------------
## Output PDF
## -------------------------------
pdf(
  file   = paste0(output_folder,
                  '/Dgkk_sanitycheck.Violin.top50_genes.by_stage_across_time.pdf'),
  width  = 18,
  height = 60
)

## -------------------------------
## Loop over stages
## -------------------------------
for (stage_i in stages) {
  
  message('Plotting genes enriched in: ', stage_i)
  
  ## ---- select top 50 genes ENRICHED in this stage
  top_genes_stage <- Dgkk_markers %>%
    filter(
      stage == stage_i,
      p_val_adj < 0.05,
      avg_log2FC > 0.5
    ) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 50) %>%
    pull(gene)
  
  ## ---- keep genes present in object
  top_genes_stage <- intersect(
    top_genes_stage,
    rownames(so)
  )
  
  if (length(top_genes_stage) == 0) {
    message('No genes to plot for ', stage_i)
    next
  }
  
  ## ---- violin plot across ALL stages
  p <- VlnPlot(
    so,
    features = top_genes_stage,
    group.by = 'Dgkk_status',
    split.by = 'orig.ident',
    pt.size  = 0,
    stack    = TRUE,
    flip     = TRUE
  ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9),
      plot.title  = element_text(face = 'bold')
    ) +
    labs(
      x = 'SCT-normalized expression',
      y = 'Gene',
      title = paste0(
        'Top 50 enriched genes in Dgkk+ cells ',
        stage_i,
        ' (shown across all stages)'
      )
    )
  
  print(p)
}

dev.off()


# =====================================================
# Dgkk GO ANALYSIS — BY STAGE (MANDIBLE E9.5–E11.5)
# =====================================================

library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

# =========================
# OUTPUT FOLDER
# =========================
output_folder <- "~/Documents/postdoc/collaboration/Pauline/Dgkk_project/results_for_publication/"
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# =========================
# LOAD DATA
# =========================
so <- readRDS("~/Documents/postdoc/collaboration/Pauline/mandible_E9.5-11.5_integration/so_mandible_E9.5_E10.5_E11.5_integrated.rds")

DefaultAssay(so) <- "SCT"

# Ensure SCT consistency once globally
so <- PrepSCTFindMarkers(so)

# =========================
# DEFINE STAGES
# =========================
stages <- c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5")

# =========================
# OUTPUT PDF (6 PANELS TOTAL)
# =========================
pdf(file.path(output_folder, "GO_Dgkk_by_stage.pdf"),
    width = 8, height = 12)

# =====================================================
# LOOP OVER STAGES
# =====================================================
for (st in stages) {
  
  cat("\n============================\n")
  cat("Processing:", st, "\n")
  
  # -------------------------
  # SUBSET STAGE
  # -------------------------
  so_stage <- subset(
    so,
    subset = orig.ident == st
  )
  
  cat("Cells in stage:\n")
  print(table(so_stage$orig.ident))
  
  if (ncol(so_stage) < 50) {
    cat("Skipping stage (too few cells)\n")
    next
  }
  
  # -------------------------
  # REBUILD SCT (CRITICAL FIX)
  # -------------------------
  so_stage <- SCTransform(so_stage, verbose = FALSE)
  so_stage <- PrepSCTFindMarkers(so_stage)
  
  # -------------------------
  # DEFINE Dgkk STATUS
  # -------------------------
  if (!"Dgkk" %in% rownames(so_stage)) {
    stop("Dgkk not found in dataset")
  }
  
  dgkk_vals <- FetchData(so_stage, vars = "Dgkk")[,1]
  
  cutoff <- quantile(dgkk_vals, 0.9, na.rm = TRUE)
  
  so_stage$Dgkk_status <- ifelse(
    dgkk_vals > cutoff,
    "Dgkk_pos",
    "Dgkk_neg"
  )
  
  cat("Dgkk distribution:\n")
  print(table(so_stage$Dgkk_status))
  
  if (length(unique(so_stage$Dgkk_status)) < 2) {
    cat("Skipping stage (no variation)\n")
    next
  }
  
  # -------------------------
  # DE ANALYSIS
  # -------------------------
  Idents(so_stage) <- "Dgkk_status"
  
  markers <- FindMarkers(
    so_stage,
    ident.1 = "Dgkk_pos",
    ident.2 = "Dgkk_neg",
    min.pct = 0.1,
    logfc.threshold = 0.25
  )
  
  # -------------------------
  # FILTER GENES
  # -------------------------
  genes_up <- rownames(markers %>%
                         filter(p_val_adj < 0.05, avg_log2FC > 0.5))
  
  genes_down <- rownames(markers %>%
                           filter(p_val_adj < 0.05, avg_log2FC < -0.5))
  
  # -------------------------
  # SKIP EMPTY
  # -------------------------
  if (length(genes_up) == 0 & length(genes_down) == 0) next
  
  # -------------------------
  # ENTREZ CONVERSION
  # -------------------------
  gene_up_entrez <- bitr(
    genes_up,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
  
  gene_down_entrez <- bitr(
    genes_down,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
  
  # -------------------------
  # GO ENRICHMENT
  # -------------------------
  ego_up <- NULL
  ego_down <- NULL
  
  if (nrow(gene_up_entrez) > 0) {
    ego_up <- enrichGO(
      gene = gene_up_entrez$ENTREZID,
      OrgDb = org.Mm.eg.db,
      ont = "BP",
      readable = TRUE
    )
  }
  
  if (nrow(gene_down_entrez) > 0) {
    ego_down <- enrichGO(
      gene = gene_down_entrez$ENTREZID,
      OrgDb = org.Mm.eg.db,
      ont = "BP",
      readable = TRUE
    )
  }
  
  # -------------------------
  # PLOTTING (6 PANELS TOTAL)
  # -------------------------
  if (!is.null(ego_up) && nrow(ego_up) > 0) {
    print(
      dotplot(ego_up, showCategory = 20) +
        ggtitle(paste0("Dgkk+ ", st))
    )
  }
  
  if (!is.null(ego_down) && nrow(ego_down) > 0) {
    print(
      dotplot(ego_down, showCategory = 20) +
        ggtitle(paste0("Dgkk- ", st))
    )
  }
}

dev.off()

cat("\nDONE — 6-panel GO figure generated.\n")
