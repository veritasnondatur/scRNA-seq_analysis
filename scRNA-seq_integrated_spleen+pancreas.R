#### Comparative analysis of integrated pancreas_E14.5 and spleen_E15.5 datasets

# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/collaboration/Maurizio/WIP_scRNA-seq_integrated_spleen+pancreas/results/"

################################################################################
## Visualize integrated dataset with high resolution for clusteridentification

# Load Seurat object
so_spleenE15.5_pancreasE14.5_integrated_old <- readRDS("/Users/veralaub/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/integrated_analysis_spleen+pancreas/int_appr3.so_spleenE15.5_pancreasE14.5_integrated.rds")

# Set the active assay to SCT (to plot normalized expression)
DefaultAssay(so_spleenE15.5_pancreasE14.5_integrated_old) <- "SCT"

# Definition of GOIs
goi <- c("Tlx1", "Barx1", "Klf4", "Nr2f2", "Nkx2-5",
         "Fgf8", "Fgf9", "Fgfr1", "Fgfr2", "Fgfr3", 
         "Cdk1", "Cdkn1c", "Mki67",
         "Epcam", "Upk3b", "Lum")
outFile <- paste(output_folder,
                 "/int_appr3.UMAP.spleenE15.5_pancreasE14.5_integrated_old.goi.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 12, height = 5)

# Loop through each gene and check if it exists in the Seurat object
for (gene in goi) {
  if (gene %in% rownames(so_spleenE15.5_pancreasE14.5_integrated_old)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_spleenE15.5_pancreasE14.5_integrated_old, 
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
######### Re-analysis of spleen_E15.5 dataset with different resolution ########

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Load E15.5 spleen data (filtered)
d_spleenE15.5 <- Read10X(data.dir ="/Users/veralaub/Documents/postdoc/collaboration/Maurizio/E15.5_spleen/MR3_filtered_feature_bc_matrix/")

# Create a Seurat object
so_spleenE15.5 <- CreateSeuratObject(counts = d_spleenE15.5,
                                     project = "spleen_E15.5",
                                     min.cells = 3, min.features = 200)

# Add mito fraction
so_spleenE15.5[["percent.mt"]] <- PercentageFeatureSet(so_spleenE15.5,
                                                       pattern = "^mt-")

# QC stats before filtering
p <- RidgePlot(so_spleenE15.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "spleenE15.5.data.qc.unfiltered.pdf", sep = "")
pdf(out_path, width = 4, height = 4)
plot(p)
dev.off()

# Define thresholds for filtering cells (can be adapted ad gusto)
percent.mt_max <- 5  # maximum percentage of mitochondrial genes (adjust as needed)
nFeature_RNA_min <- 1000  # minimum number of features per cell
nFeature_RNA_max <- 10000 # minimum number of features per cell
nCount_RNA_min <- 100  # minimum number of RNA counts per cell
nCount_RNA_max <- 50000  # maximum number of RNA counts per cell

# Subset the Seurat object (filter based on thresholds above)
so_spleenE15.5_filtered <- subset(so_spleenE15.5,
                                  subset = percent.mt <= percent.mt_max &
                                  nFeature_RNA >= nFeature_RNA_min &
                                  nFeature_RNA <= nFeature_RNA_max &
                                  nCount_RNA >= nCount_RNA_min &
                                  nCount_RNA <= nCount_RNA_max)

# QC stats after filtering
p <- RidgePlot(so_spleenE15.5_filtered,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "spleenE15.5.data.qc.filtered.pdf", sep = "")
pdf(out_path, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_spleenE15.5_filtered[["RNA"]]), invert = TRUE)

# Exclude mito genes
filtered_counts <- GetAssayData(so_spleenE15.5_filtered[["RNA"]], layer = "counts")[keep, ]

# Create Seurat object with filtered genes
# Access counts using GetAssayData instead of directly from counts slot
filtered_counts <- GetAssayData(so_spleenE15.5_filtered[["RNA"]], layer = "counts")[keep, ]

# Create a new Seurat object after filtering
so_spleenE15.5_filtered <- CreateSeuratObject(counts = filtered_counts, 
                                              project = "spleen_E15.5", 
                                              meta.data = so_spleenE15.5_filtered@meta.data)

# Save the Seurat object
#so_path_spleenE15.5 <- paste(output_folder, "so_spleenE15.5_beforeNormalization.rds", sep = "")
#saveRDS(so_spleenE15.5_filtered, file = so_path_spleenE15.5)

# Normalization with SCTransform (has not changed in Seurat v5)
n_features <- 2000
so_spleenE15.5_filtered_norm <- SCTransform(so_spleenE15.5_filtered,
                                            verbose = FALSE,
                                            variable.features.n = n_features)

## Cell Cycle regression
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
so_spleenE15.5_filtered_norm <- CellCycleScoring(so_spleenE15.5_filtered_norm,
                                                  s.features = cc.genes_corrected$s.genes,   # S phase genes
                                                  g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                                  set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling
so_spleenE15.5_filtered_norm <- ScaleData(so_spleenE15.5_filtered_norm,
                                          features = rownames(so_spleenE15.5_filtered_norm),
                                          vars.to.regress = c("S.Score", "G2M.Score"),
                                          verbose = TRUE
)

# PCA (no changes for PCA)
so_spleenE15.5_filtered_norm <- RunPCA(so_spleenE15.5_filtered_norm,
                                       verbose = FALSE, npcs = 20)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_spleenE15.5_filtered_norm)
out_path <- paste(output_folder, "spleenE15.5.data.qc.ellbowplot.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

# Store number of principle components in new variable (to be used later)
pca_dim_sel <- 8

# UMAP
so_spleenE15.5_filtered_norm <- RunUMAP(so_spleenE15.5_filtered_norm,
                                        dims = 1:pca_dim_sel)

# Save the Seurat object
so_path_spleenE15.5 <- paste(output_folder, "so_spleenE15.5_norm_CCR.rds", sep = "")
saveRDS(so_spleenE15.5_filtered_norm, file = so_path_spleenE15.5)

# Load E15.5 spleen data (filtered)
so_spleenE15.5_filtered_norm <- readRDS("/Users/veralaub/Documents/postdoc/collaboration/Maurizio/WIP_scRNA-seq_integrated_spleen+pancreas/results/so_spleenE15.5_norm_CCR.rds")

# Clustering (Leiden algorithm)
so_spleenE15.5_filtered_norm <- FindNeighbors(so_spleenE15.5_filtered_norm,
                                              dims = 1:pca_dim_sel)
so_spleenE15.5_filtered_norm <- FindClusters(so_spleenE15.5_filtered_norm,
                                             resolution = 0.3,
                                             algorithm = 4)

# Visualize clusters as Dimplot
p <- DimPlot(object = so_spleenE15.5_filtered_norm,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/spleenE15.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_spleenE15.5_filtered_norm, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/spleenE15.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save the Seurat object
so_path_spleenE15.5 <- paste(output_folder, "so_spleenE15.5.rds", sep = "")
saveRDS(so_spleenE15.5_filtered_norm, file = so_path_spleenE15.5)


################################################################################
######## Re-analysis of pancreas_E14.5 dataset with different resolution #######

# Load pancreas data (WT only)
d_pancreasE14.5 <- Read10X(data.dir ="/Users/veralaub/Documents/postdoc/collaboration/Maurizio/Fgf9null_datasets/scRNA-seq/pancreas_Fgf9null/E14.5_pancreas_Fgf9WT_GSM6434043/")

# Create a Seurat object
so_pancreasE14.5 <- CreateSeuratObject(counts = d_pancreasE14.5, 
                                         project = "pancreas_E14.5",
                                         min.cells = 3, min.features = 200)

# Add mito fraction (no change here)
so_pancreasE14.5[["percent.mt"]] <- PercentageFeatureSet(so_pancreasE14.5,
                                                         pattern = "^mt-")

# QC stats before filtering
p <- RidgePlot(so_pancreasE14.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "pancreasE14.5.data.qc.unfiltered.pdf", sep = "")
pdf(out_path, width = 4, height = 4)
plot(p)
dev.off()

# Define thresholds for filtering cells (can be adapted ad gusto)
percent.mt_max <- 5  # maximum percentage of mitochondrial genes (adjust as needed)
nFeature_RNA_min <- 2500  # minimum number of features per cell
nFeature_RNA_max <- 15000 # minimum number of features per cell
nCount_RNA_min <- 100  # minimum number of RNA counts per cell
nCount_RNA_max <- 50000  # maximum number of RNA counts per cell

# Subset the Seurat object (filter based on thresholds above)
so_pancreasE14.5_filtered <- subset(so_pancreasE14.5,
                                    subset = percent.mt <= percent.mt_max &
                                    nFeature_RNA >= nFeature_RNA_min &
                                    nFeature_RNA <= nFeature_RNA_max &
                                    nCount_RNA >= nCount_RNA_min &
                                    nCount_RNA <= nCount_RNA_max)

# QC stats after filtering
p <- RidgePlot(so_pancreasE14.5_filtered,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "pancreasE14.5.data.qc.filtered.pdf", sep = "")
pdf(out_path, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_pancreasE14.5_filtered[["RNA"]]), invert = TRUE)

# Exclude mito genes
filtered_counts <- GetAssayData(so_pancreasE14.5_filtered[["RNA"]], layer = "counts")[keep, ]

# Create Seurat object with filtered genes
# Access counts using GetAssayData instead of directly from counts slot
filtered_counts <- GetAssayData(so_pancreasE14.5_filtered[["RNA"]], layer = "counts")[keep, ]

# Create a new Seurat object after filtering
so_pancreasE14.5_filtered <- CreateSeuratObject(counts = filtered_counts,
                                                project = "pancreas_E14.5", 
                                                meta.data = so_pancreasE14.5_filtered@meta.data)

# Save the Seurat object
#so_path_pancreasE14.5 <- paste(output_folder, "so_pancreasE14.5_beforeNormalization.rds", sep = "")
#saveRDS(so_pancreasE14.5_filtered, file = so_path_pancreasE14.5)

# Normalization with SCTransform (has not changed in Seurat v5)
n_features <- 2000
so_pancreasE14.5_filtered_norm <- SCTransform(so_pancreasE14.5_filtered,
                                              verbose = TRUE,
                                              variable.features.n = n_features)

## Cell Cycle regression
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
so_pancreasE14.5_filtered_norm <- CellCycleScoring(so_pancreasE14.5_filtered_norm,
                                                   s.features = cc.genes_corrected$s.genes,   # S phase genes
                                                   g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                                   set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling
so_pancreasE14.5_filtered_norm <- ScaleData(so_pancreasE14.5_filtered_norm,
                                            features = rownames(so_pancreasE14.5_filtered_norm),
                                            vars.to.regress = c("S.Score", "G2M.Score"),
                                            verbose = TRUE
)

# PCA
so_pancreasE14.5_filtered_norm <- RunPCA(so_pancreasE14.5_filtered_norm,
                                         verbose = FALSE, npcs = 20)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_pancreasE14.5_filtered_norm, ndims = 20)
out_path <- paste(output_folder, "pancreasE14.5.data.qc.ellbowplot.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

# Store number of principle components in new variable (to be used later)
pca_dim_sel <- 8

# UMAP
so_pancreasE14.5_filtered_norm <- RunUMAP(so_pancreasE14.5_filtered_norm,
                                          dims = 1:pca_dim_sel)

# Save the Seurat object
so_path_pancreasE14.5 <- paste(output_folder, "so_pancreasE14.5_norm_CCR.rds", sep = "")
saveRDS(so_pancreasE14.5_filtered_norm, file = so_path_pancreasE14.5)

# Load E15.5 spleen data (filtered)
so_pancreasE14.5_filtered_norm <- readRDS("/Users/veralaub/Documents/postdoc/collaboration/Maurizio/WIP_scRNA-seq_integrated_spleen+pancreas/results/so_pancreasE14.5_norm_CCR.rds")

# Clustering (Leiden) - Seurat v5 should work similarly
so_pancreasE14.5_filtered_norm <- FindNeighbors(so_pancreasE14.5_filtered_norm,
                                                dims = 1:pca_dim_sel)
so_pancreasE14.5_filtered_norm <- FindClusters(so_pancreasE14.5_filtered_norm,
                                                resolution = 0.3,
                                                algorithm = 4)

# Visualize clusters as Dimplot
p <- DimPlot(object = so_pancreasE14.5_filtered_norm,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/pancreasE14.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_pancreasE14.5_filtered_norm, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/pancreasE14.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save the Seurat object
so_path_pancreasE14.5 <- paste(output_folder, "so_pancreasE14.5.rds", sep = "")
saveRDS(so_pancreasE14.5_filtered_norm, file = so_path_pancreasE14.5)


################################################################################
######################## INTEGRATION OF THE TWO DATASETS #######################

# Load raw data (pre-processed Seurat objects produced above)
so_spleenE15.5 <- so_spleenE15.5_filtered_norm
so_spleenE15.5 <- readRDS("/Users/veralaub/Documents/postdoc/collaboration/Maurizio/WIP_scRNA-seq_integrated_spleen+pancreas/results/so_spleenE15.5.rds")

so_pancreasE14.5 <- so_pancreasE14.5_filtered_norm
so_pancreasE14.5 <- readRDS("/Users/veralaub/Documents/postdoc/collaboration/Maurizio/WIP_scRNA-seq_integrated_spleen+pancreas/results/so_pancreasE14.5.rds")

# Change the default assay to "SCT"
DefaultAssay(so_spleenE15.5) <- "SCT"
DefaultAssay(so_pancreasE14.5) <- "SCT"

# Check memory usage before and after merge (due to Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Identify variable features for both datasets
so_spleenE15.5 <- FindVariableFeatures(so_spleenE15.5, 
                                       selection.method = "vst", 
                                       nfeatures = 2000)
so_pancreasE14.5 <- FindVariableFeatures(so_pancreasE14.5, 
                                         selection.method = "vst", 
                                         nfeatures = 2000)

# Select integration features
SelectIntegrationFeatures(object.list = list(so_spleenE15.5, so_pancreasE14.5),
                          nfeatures = 2000,
                          verbose = TRUE)

# Step 1: Get the variable features for both datasets
var_features_spleen <- so_spleenE15.5@assays[["SCT"]]@var.features
var_features_pancreas <- so_pancreasE14.5@assays[["SCT"]]@var.features

# Save Seurat objects (interim step to retrieve later [necessary due to Memory capacity error])
outFile <- paste(output_folder, "so_spleenE15.5_integrationWIP.rds", sep = "")
saveRDS(so_spleenE15.5, file = outFile)
outFile <- paste(output_folder, "so_pancreasE14.5_integrationWIP.rds", sep = "")
saveRDS(so_pancreasE14.5, file = outFile)

# Retrieve Seurat objects (interim step to retrieve later [necessary due to Memory capacity error])
#so_spleenE15.5 <- readRDS("/Users/veralaub/Documents/postdoc/collaboration/Maurizio/WIP_scRNA-seq_integrated_spleen+pancreas/results/so_spleenE15.5_integrationWIP.rds")
#so_pancreasE14.5 <- readRDS("/Users/veralaub/Documents/postdoc/collaboration/Maurizio/WIP_scRNA-seq_integrated_spleen+pancreas/results/so_pancreasE14.5_integrationWIP.rds")

# Step 2: Find the common variable features between the two datasets
common_var_features <- intersect(var_features_spleen, var_features_pancreas)

# Step 3: Prepare the objects for integration using the common features
objects <- list(so_spleenE15.5, so_pancreasE14.5)

# Prepare objects for integration
objects <- PrepSCTIntegration(object.list = objects, 
                              anchor.features = common_var_features, 
                              verbose = TRUE)

# Step 4: Find integration anchors - make sure to specify the common features
anchors <- FindIntegrationAnchors(object.list = objects, 
                                  normalization.method = "SCT", 
                                  dims = 1:30, 
                                  anchor.features = common_var_features,  # explicitly specify the features
                                  verbose = TRUE)

# Save IntegrationAnchorSet (interim step to retrieve later)
# [necessary due to Memory capacity error]
outFile <- paste(output_folder, "IntegrationAnchorSet_spleenE15.5_pancreasE14.5.rds", sep = "")
saveRDS(anchors, file = outFile)

# Clear objects from workspace or clear all/close R and retrieve saved IntegrationAnchorSet
# [necessary due to Memory capacity error]
rm(so_spleenE15.5)
rm(so_pancreasE14.5)
rm(objects)
rm(so_pancreasE14.5_filtered)
rm(so_pancreasE14.5_filtered_norm)
rm(so_spleenE15.5_filtered)
rm(so_spleenE15.5_filtered_norm)
anchors <- readRDS("/Users/veralaub/Documents/postdoc/collaboration/Maurizio/WIP_scRNA-seq_integrated_spleen+pancreas/results/IntegrationAnchorSet_spleenE15.5_pancreasE14.5.rds")

# Step 5: Integrate the datasets using the found anchors
so_spleenE15.5_pancreasE14.5_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

# Step 6: Perform scaling and PCA on the integrated data
so_spleenE15.5_pancreasE14.5_integrated <- ScaleData(so_spleenE15.5_pancreasE14.5_integrated)
so_spleenE15.5_pancreasE14.5_integrated <- RunPCA(so_spleenE15.5_pancreasE14.5_integrated, verbose = FALSE)

# Perform UMAP on integrated data
so_spleenE15.5_pancreasE14.5_integrated <- RunUMAP(so_spleenE15.5_pancreasE14.5_integrated, 
                                                   dims = 1:30, 
                                                   reduction = "pca", 
                                                   reduction.name = "umap.integrated")

# Change the 'orig.ident' metadata of E15.5 spleen (to match name)
so_spleenE15.5_pancreasE14.5_integrated$orig.ident <- gsub("^MR4", "spleen_E15.5", so_spleenE15.5_pancreasE14.5_integrated$orig.ident)

# Visualize datasets as UMAP after Integration
p <- DimPlot(so_spleenE15.5_pancreasE14.5_integrated, 
             reduction = "umap.integrated", 
             group.by = c("orig.ident", "seurat_clusters"))
outFile <- paste(output_folder, "/spleenE15.5_pancreasE14.5_integrated.UMAP.orig.ident.clusters.pdf", sep = "")
pdf(outFile, width = 12, height = 5)
plot(p)
dev.off()

# Save the Seurat object
outFile <- paste(output_folder, "so_spleenE15.5_pancreasE14.5_integrated.rds", sep = "")
saveRDS(so_spleenE15.5_pancreasE14.5_integrated, file = outFile)

### Visualization
# Change the default assay to "SCT"
DefaultAssay(so_spleenE15.5_pancreasE14.5_integrated) <- "SCT"

# Visualize as FeaturePlot
FeaturePlot(so_spleenE15.5_pancreasE14.5_integrated, 
            features = "Tlx1",
            reduction = "umap.integrated",
            split.by = "orig.ident")

goi <- c("Tlx1", "Barx1", "Klf4", "Nr2f2", "Nkx2-5",
         "Fgf8", "Fgf9", "Fgfr1", "Fgfr2", "Fgfr3", 
         "Cdk1", "Cdkn1c", "Mki67",
         "Epcam", "Upk3b", "Lum")
outFile <- paste(output_folder,
                 "/spleenE15.5_pancreasE14.5_integrated.UMAP.goi.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 12, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi) {
  if (gene %in% rownames(so_spleenE15.5_pancreasE14.5_integrated)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_spleenE15.5_pancreasE14.5_integrated, 
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
######## MARKER ANALYSIS (following meeting with Maurizio on 02/27/2025) #######

# Change the default assay to "SCT"
DefaultAssay(so_spleenE15.5_pancreasE14.5_integrated) <- "SCT"

## Identify top 25 markers + all per cluster
# Correcting SCT counts before running FindAllMarkers
so_spleenE15.5_pancreasE14.5_integrated <- PrepSCTFindMarkers(so_spleenE15.5_pancreasE14.5_integrated,
                                                              assay = "SCT", 
                                                              verbose = TRUE)
markers <- FindAllMarkers(so_spleenE15.5_pancreasE14.5_integrated,
                          min.pct = 0.1,
                          test.use = "wilcox")

# View markers for all clusters
marker_table <- table(markers$cluster)  # Shows the number of markers for each cluster

# Show markers for the first few clusters
cluster_ids <- unique(so_spleenE15.5_pancreasE14.5_integrated$seurat_clusters)   # Get unique cluster identities
num_clusters <- length(cluster_ids)   # Count the number of unique clusters

# Loop over each cluster to extract and print top 25 markers
for (cluster in cluster_ids) {
  # Extract the top 25 markers for this cluster
  top_markers <- head(markers[markers$cluster == cluster, ], 25)  # Get top 25 markers for the current cluster
  # Print the top 20 markers for the current cluster
  cat("Top 25 markers for cluster ", cluster, " are: \n", sep = "")
  # Print the gene names (marker genes) for the current cluster
  print(top_markers$gene)   # Assuming 'gene' is the column containing marker gene names
  cat("\n")  # Add a line break between clusters
}

# Save the marker list to a CSV file
write.csv(markers, file = "/Users/veralaub/Documents/postdoc/collaboration/Maurizio/WIP_scRNA-seq_integrated_spleen+pancreas/results/spleenE15.5_pancreasE14.5_integrated_markers_by_cluster.csv", row.names = TRUE)


################################################################################
########### Exploratory post-hoc analysis to explore integrated dataset #########
# Can be run from here without preloading any of the other datasets

# Load data
so_spleenE15.5_pancreasE14.5_integrated <- readRDS("/Users/veralaub/Documents/postdoc/collaboration/Maurizio/WIP_scRNA-seq_integrated_spleen+pancreas/results/so_spleenE15.5_pancreasE14.5_integrated.rds")

# Change the default assay to "SCT" (normalized dataset)
DefaultAssay(so_spleenE15.5_pancreasE14.5_integrated) <- "SCT"

# Visualize as FeaturePlot
FeaturePlot(so_spleenE15.5_pancreasE14.5_integrated, 
            features = "Tlx1",
            reduction = "umap.integrated",
            split.by = "orig.ident")

