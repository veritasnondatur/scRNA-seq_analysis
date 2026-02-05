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

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/collaboration/Pauline/Gpr50_project/results/"

# Define genes of interest (goi)
goi <- c("Gpr50", "Hand2", "Dlx6", "Pcsk1n")


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
# Hand2 + Gpr50
p1 <- FeaturePlot(
  so_mandible_E9.5_E10.5_E11.5_integrated,
  features = c("Hand2", "Gpr50"),
  reduction = "umap.integrated",
  blend = TRUE,                  # blend colors
  cols = c("grey90", "blue", "darkgreen"),
  order = TRUE                   # plot higher expressers on top
)

# Dlx6 + Gpr50
p2 <- FeaturePlot(
  so_mandible_E9.5_E10.5_E11.5_integrated,
  features = c("Dlx6", "Gpr50"),
  reduction = "umap.integrated",
  blend = TRUE,
  cols = c("grey90", "blue", "darkgreen"),
  order = TRUE
)

# Save to PDF
outFile <- paste0(output_folder, "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.Hand2_Dlx6_Gpr50.pdf")
pdf(outFile, width = 25, height = 7)
print(p1)
print(p2)
dev.off()


### Overlaps at each timepoint separated
# Define pairs of interest
gene_pairs <- list(
  c("Hand2", "Gpr50"),
  c("Dlx6", "Gpr50")
)

# Output file
outFile <- paste0(output_folder,
                  "/mandible_E9.5_E10.5_E11.5_integrated.UMAP.Hand2_Dlx6_Gpr50.orig.ident.pdf")
pdf(outFile, width = 20, height = 10)

for (pair in gene_pairs) {
  if (all(pair %in% rownames(so_mandible_E9.5_E10.5_E11.5_integrated))) {
    p <- FeaturePlot(
      so_mandible_E9.5_E10.5_E11.5_integrated,
      features = pair,
      reduction = "umap.integrated",
      blend = TRUE,
      cols = c("grey90", "blue", "darkgreen"),
      order = TRUE,
      split.by = "orig.ident"
    )
    print(p)
  } else {
    message(paste("One or both genes not found in data: ", paste(pair, collapse = ", ")))
  }
}

dev.off()

################################################################################
################# Exploratory analysis after meeting 11/20/2025 ################

######## MARKER ANALYSIS 

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
  top_markers <- head(markers[markers$cluster == cluster, ], 100)  # Get top 25 markers for the current cluster
  # Print the top 20 markers for the current cluster
  cat("Top 100 markers for cluster ", cluster, " are: \n", sep = "")
  # Print the gene names (marker genes) for the current cluster
  print(top_markers$gene)   # Assuming 'gene' is the column containing marker gene names
  cat("\n")  # Add a line break between clusters
}

# Save the marker list to a CSV file
write.csv(markers, file = "~/Documents/postdoc/collaboration/Pauline/Gpr50_project/results/mandible_E9.5_E10.5_E11.5_integrated_markers_by_cluster.csv", row.names = TRUE)


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

###### Renaming clusters according to GO term analysis of marker gene expression
### Still ongoing (as of 02/03/2026)

# Example mapping: cluster ID -> new name
new_cluster_names <- c(
  "1" = "",
  "2" = "",
  "3" = "",
  "4" = "",
  "5" = "",
  "6" = "",
  "7" = "",
  "8" = "",
  "9" = "",
  "10" = "",
  "11" = ""
)

################################################################################
# Within-stage differential expression analysis, comparing Gpr50⁺ vs Gpr50⁻ cells 
# separately for E9.5, E10.5, and E11.5

# Load integrated Dataset
so_mandible_E9.5_E10.5_E11.5_integrated <- readRDS("~/Documents/postdoc/collaboration/Pauline/mandible_E9.5-11.5_integration/so_mandible_E9.5_E10.5_E11.5_integrated.rds")

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/collaboration/Pauline/Gpr50_project/results/"

# Prepare the Seurat object
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Make new Seurat object for plots
so <- so_mandible_E9.5_E10.5_E11.5_integrated

DefaultAssay(so) <- "SCT"
so <- PrepSCTFindMarkers(so, verbose = TRUE)


# Define Gpr50 positivity
so$Gpr50_status <- ifelse(
  GetAssayData(so, layer = "data")["Gpr50", ] > 0,
  "Gpr50_pos",
  "Gpr50_neg"
)

so$Gpr50_status <- factor(
  so$Gpr50_status,
  levels = c("Gpr50_neg", "Gpr50_pos")
)


# Differential expression per developmental stage
run_gpr50_DE <- function(seurat_obj, stage) {
  
  message("Running DE for ", stage)
  
  so_sub <- subset(
    seurat_obj,
    subset = orig.ident == stage
  )
  
  Idents(so_sub) <- "Gpr50_status"
  
  markers <- FindMarkers(
    so_sub,
    ident.1 = "Gpr50_pos",
    ident.2 = "Gpr50_neg",
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

gpr50_markers <- lapply(stages, function(s) run_gpr50_DE(so, s))
gpr50_markers <- bind_rows(gpr50_markers)


# Save DE results (logFC + FDR)
write.csv(
  gpr50_markers,
  file = paste0(
    output_folder,
    "/Gpr50_pos_vs_neg.DE.by_stage.csv"
  ),
  row.names = FALSE
)


### Heatmap of top enriched genes per stage
# Heatmap, stratified by Gpr50 status and stage

# Ensure correct assay
DefaultAssay(so) <- "SCT"

# Ensure stage order (left → right)
so$orig.ident <- factor(
  so$orig.ident,
  levels = c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5")
)

# Sanity check: Gpr50_status must already exist
table(so$Gpr50_status)
# should show Gpr50_pos / Gpr50_neg

# Create composite grouping: stage × Gpr50 status
so$stage_Gpr50 <- paste0(
  so$orig.ident, "_", so$Gpr50_status
)

# Explicit ordering: stage first, then Gpr50− → Gpr50+
so$stage_Gpr50 <- factor(
  so$stage_Gpr50,
  levels = c(
    "mandible_E9.5_Gpr50_neg", "mandible_E9.5_Gpr50_pos",
    "mandible_E10.5_Gpr50_neg", "mandible_E10.5_Gpr50_pos",
    "mandible_E11.5_Gpr50_neg", "mandible_E11.5_Gpr50_pos"
  )
)

# Select top enriched genes per stage (from your DE results)
top_genes <- gpr50_markers %>%
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
    "/Gpr50_top_enriched_genes.Heatmap.by_stage_and_Gpr50_status.pdf"
  ),
  width = 14,   # wider to avoid clipping
  height = 0.35 * length(top_genes) + 4
)

DoHeatmap(
  so,
  features = top_genes,
  group.by = "stage_Gpr50"
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
    title = "Top genes enriched in Gpr50⁺ cells\nstratified by developmental stage",
    fill = "Scaled\nexpression"   # legend title
  )

dev.off()


### DotPlot (stage × Gpr50 status)
# Explicit ordering: stage first, then Gpr50− → Gpr50+
so$stage_Gpr50 <- factor(
  so$stage_Gpr50,
  levels = c(
    "mandible_E9.5_Gpr50_neg", "mandible_E9.5_Gpr50_pos",
    "mandible_E10.5_Gpr50_neg", "mandible_E10.5_Gpr50_pos",
    "mandible_E11.5_Gpr50_neg", "mandible_E11.5_Gpr50_pos"
  )
)

Idents(so) <- "stage_Gpr50"

pdf(
  paste0(output_folder, "/Gpr50_top_enriched_genes.DotPlot.by_stage_and_Gpr50_status.pdf"),
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
    y = "Stage × Gpr50 status",
    color = "Avg. expression",
    size = "% expressing"
  ) +
  theme(axis.text.x = element_text(angle = 90))

dev.off()

########## Violin plots
# Sanity-check violin plots: Gpr50+ vs Gpr50−, split by stage
# Select top genes to visualize
genes_to_plot <- head(top_genes, 50)

outFile <- paste0(
  output_folder,
  "/Gpr50_sanitycheck.Violin.top50_genes.by_stage.pdf"
)

pdf(
  outFile,
  width  = 18,
  height = 1 * length(genes_to_plot) + 4
)

p <- VlnPlot(
  so,
  features = genes_to_plot,
  group.by = "Gpr50_status",
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
    title = "Top Gpr50⁺-enriched genes (sanity check)"
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
# gpr50_markers
# output_folder

## -------------------------------
## Load integrated object ONCE
## -------------------------------
so <- so_mandible_E9.5_E10.5_E11.5_integrated
DefaultAssay(so) <- "SCT"

## -------------------------------
## Ensure Gpr50_status exists
## (safe even if it already does)
## -------------------------------
if (!"Gpr50_status" %in% colnames(so@meta.data)) {
  message("Adding Gpr50_status metadata")
  
  so$Gpr50_status <- ifelse(
    GetAssayData(so, assay = "SCT", layer = "data")["Gpr50", ] > 0,
    "Gpr50_pos",
    "Gpr50_neg"
  )
  
  so$Gpr50_status <- factor(
    so$Gpr50_status,
    levels = c("Gpr50_neg", "Gpr50_pos")
  )
}

## Quick sanity check
print(table(so$Gpr50_status, useNA = "ifany"))

## -------------------------------
## Stages in chronological order
## -------------------------------
stages <- c("mandible_E9.5", "mandible_E10.5", "mandible_E11.5")

## -------------------------------
## Open PDF
## -------------------------------
pdf(
  file   = paste0(output_folder,
                  "/Gpr50_sanitycheck.Violin.top50_genes.by_stage.pdf"),
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
  top_genes_stage <- gpr50_markers %>%
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
    group.by = "Gpr50_status",
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
        "Top 50 Gpr50+ enriched genes — ",
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
## Ensure Gpr50_status exists
## -------------------------------
if (!'Gpr50_status' %in% colnames(so@meta.data)) {
  so$Gpr50_status <- ifelse(
    GetAssayData(so, assay = 'SCT', layer = 'data')['Gpr50', ] > 0,
    'Gpr50_pos',
    'Gpr50_neg'
  )
  so$Gpr50_status <- factor(
    so$Gpr50_status,
    levels = c('Gpr50_neg', 'Gpr50_pos')
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
                  '/Gpr50_sanitycheck.Violin.top50_genes.by_stage_across_time.pdf'),
  width  = 18,
  height = 60
)

## -------------------------------
## Loop over stages
## -------------------------------
for (stage_i in stages) {
  
  message('Plotting genes enriched in: ', stage_i)
  
  ## ---- select top 50 genes ENRICHED in this stage
  top_genes_stage <- gpr50_markers %>%
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
    group.by = 'Gpr50_status',
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
        'Top 50 Gpr50+ enriched genes in ',
        stage_i,
        ' (shown across all stages)'
      )
    )
  
  print(p)
}

dev.off()
