#### Analysis of midface E10.5 datasets (Epi only and Mes+Epi)
# THIS IS ONLY A TESTVERSION AS WYNTON IS IN FUCKUP MODE ATM AND KEEPS CRASHING (08/21/2025)
# Does exclude a lot of data to aviod memory overload

# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)


################################################################################
############### Pre-analysis of midface_Epi_Mes_E10.5 dataset  #################

# Set out folder (to store results)
output_folder <- "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/MF_Epi+Mes-E10.5ctrl_ML7/analysis/"

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 16 * 1024 * 1024 * 1024)

# Read barcodes, features (genes) and matrix files
midface_Epi_Mes_E10.5 <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/MF_Epi+Mes-E10.5ctrl_ML7/filtered_feature_bc_matrix')

# Create Seurat object
so_midface_Epi_Mes_E10.5 <- CreateSeuratObject(counts = midface_Epi_Mes_E10.5, 
                                       project = "midface_Epi_Mes_E10.5", 
                                       min.cells = 3, min.features = 1000)      # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)


######## Standard pre-processing workflow
# Add mito fraction to object meta.data
so_midface_Epi_Mes_E10.5 <- PercentageFeatureSet(so_midface_Epi_Mes_E10.5, 
                                         pattern = "^mt-", 
                                         col.name = "percent.mt")

# QC stats before filtering
p <- RidgePlot(so_midface_Epi_Mes_E10.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "midface_Epi_Mes_E10.5.qc.unfiltered.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

### Selecting cells for further analysis
percent.mt_max <- 5          # maximum percentage of mitochondrial genes (adjust as needed)
nFeature_RNA_min <- 1000      # minimum number of features per cell
nFeature_RNA_max <- 10000     # minimum number of features per cell
nCount_RNA_min <- 1000        # minimum number of RNA counts per cell
nCount_RNA_max <- 50000      # maximum number of RNA counts per cell

# Subset the Seurat object (filter based on thresholds above)
so_midface_Epi_Mes_E10.5 <- subset(so_midface_Epi_Mes_E10.5,
                           subset = percent.mt <= percent.mt_max &
                             nFeature_RNA >= nFeature_RNA_min &
                             nFeature_RNA <= nFeature_RNA_max &
                             nCount_RNA >= nCount_RNA_min &
                             nCount_RNA <= nCount_RNA_max)

### Data normalization (using SCTransform)
# QC stats after filtering
p <- RidgePlot(so_midface_Epi_Mes_E10.5,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_file <- paste(output_folder, "midface_Epi_Mes_E10.5.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_midface_Epi_Mes_E10.5[["RNA"]]), invert = TRUE)
so_midface_Epi_Mes_E10.5 <- CreateSeuratObject(counts = GetAssayData(so_midface_Epi_Mes_E10.5[["RNA"]], layer = "counts")[keep, ], 
                                       project = "midface_Epi_Mes_E10.5", 
                                       meta.data = so_midface_Epi_Mes_E10.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Save Seurat object
saveRDS(so_midface_Epi_Mes_E10.5, file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/MF_Epi+Mes-E10.5ctrl_ML7/analysis/so_midface_Epi_Mes_E10.5.rds")

# Clear Workspace and re-load Seurat object (otherwise code keeps crashing during Cell Cycle regression)
so_midface_Epi_Mes_E10.5 <- readRDS("/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/MF_Epi+Mes-E10.5ctrl_ML7/analysis/so_midface_Epi_Mes_E10.5.rds")

# Normalization with SCTransform
n_features <- 2000
so_midface_Epi_Mes_E10.5 <- SCTransform(so_midface_Epi_Mes_E10.5,
                                verbose = TRUE,
                                variable.features.n = n_features)

# Save Seurat object
saveRDS(so_midface_Epi_Mes_E10.5, file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/MF_Epi+Mes-E10.5ctrl_ML7/analysis/so_midface_Epi_Mes_E10.5.rds")

# Clear Workspace and re-load Seurat object (otherwise code keeps crashing during Cell Cycle regression)
so_midface_Epi_Mes_E10.5 <- readRDS("/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/MF_Epi+Mes-E10.5ctrl_ML7/analysis/so_midface_Epi_Mes_E10.5.rds")

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
so_midface_Epi_Mes_E10.5 <- CellCycleScoring(so_midface_Epi_Mes_E10.5,
                                     s.features = cc.genes_corrected$s.genes,   # S phase genes
                                     g2m.features = cc.genes_corrected$g2m.genes,  # G2M phase genes
                                     set.ident = TRUE
)
# This adds two new metadata columns, `S.Score` and `G2M.Score`, for each cell indicating its level of expression in the S and G2M phases.

# 3. Regress out cell cycle scores during scaling (regress cell cycle information from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis)
so_midface_Epi_Mes_E10.5 <- ScaleData(so_midface_Epi_Mes_E10.5,
                              features = rownames(so_midface_Epi_Mes_E10.5),
                              vars.to.regress = c("S.Score", "G2M.Score"),
                              verbose = TRUE
)

# PCA
so_midface_Epi_Mes_E10.5 <- RunPCA(so_midface_Epi_Mes_E10.5,
                           verbose = FALSE, 
                           npcs = 100,                # Number of principal components to compute
                           ndims.print = 1:5,         # Print details for the first 5 PCs
                           nfeatures.print = 30       # Print details for the top 30 features
)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_midface_Epi_Mes_E10.5, ndims = 20)
out_file <- paste(output_folder, "midface_Epi_Mes_E10.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method
pca_dim_sel <- 8   

#UMAP
so_midface_Epi_Mes_E10.5 <- RunUMAP(so_midface_Epi_Mes_E10.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_midface_Epi_Mes_E10.5 <- FindNeighbors(so_midface_Epi_Mes_E10.5, dims = 1:pca_dim_sel)         
so_midface_Epi_Mes_E10.5 <- FindClusters(so_midface_Epi_Mes_E10.5, 
                                 resolution = 0.6,
                                 algorithm = 4)            

# Visualize clusters as Dimplot
p <- DimPlot(object = so_midface_Epi_Mes_E10.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_path <- paste(output_folder, "/midface_Epi_Mes_E10.5.UMAP.Dimplot.clusters.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Visualize the S and G2M scores on UMAP
p <- FeaturePlot(so_midface_Epi_Mes_E10.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_path <- paste(output_folder, "/midface_Epi_Mes_E10.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_path, width = 15, height = 10)
plot(p)
dev.off()

# Save Seurat object
saveRDS(so_midface_Epi_Mes_E10.5, file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/MF_Epi+Mes-E10.5ctrl_ML7/analysis/so_midface_Epi_Mes_E10.5.rds")


###################################  Plots  ####################################

# Define genes of interest (goi) to mark Epi+Mes and BMP-signalling pathway
goi <- c("Epcam", "Krt18",                                                      # Epithelial markers
         "Pdgfra", "Twist1", "Prrx1","Prrx2", "Snai1", "Alx4", "Col13a1",                              # Mesenchymal markers
         "Bmp4", "Bmpr1a", "Bmpr1b", "Bmpr2",                                   # BMP signalling
         "Smad1", "Smad2", "Smad4", "Smad5",
         "Id1", "Id2", "Id3",
         "Dcx",                                                                 # neurblast marker
         "Sox2",                                                                # neuro-epithelium
         "Sox10"                                                                # neural crest
         )


################################ Violin plots  #################################
# Violin Plot for goi
outFile <- paste(output_folder,
                 "/midface_Epi_Mes_E10.5.Violin.goi.seurat_clusters.pdf", 
                 sep = "")
pdf(outFile, width = 20, height = 5)
for (gene in goi) {
  if (gene %in% rownames(so_midface_Epi_Mes_E10.5)) {
    p <- VlnPlot(so_midface_Epi_Mes_E10.5,
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
outFile <- paste(output_folder,
                 "/midface_Epi_Mes_E10.5.UMAP.goi.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 15, height = 10)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi) {
  if (gene %in% rownames(so_midface_Epi_Mes_E10.5)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_midface_Epi_Mes_E10.5, 
                     features = gene,
                     reduction = "umap",
                     split.by = "orig.ident")
    plot(p)
  } else {
    # Print a message for missing genes (optional)
    message(paste("Gene not found in data: ", gene))
  }
}

dev.off()

### Determine marker genes for each cluster
so_midface_Epi_Mes_E10.5 <- PrepSCTFindMarkers(so_midface_Epi_Mes_E10.5)

markers <- FindAllMarkers(
  object = so_midface_Epi_Mes_E10.5,
  only.pos = TRUE,               # Return only positive markers (upregulated in the cluster)
  min.pct = 0.25,                # Gene expressed in at least 25% of cells in either group
  logfc.threshold = 0.25         # Minimum log fold change
)

# Save the marker list to a CSV file
write.csv(markers, file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/MF_Epi+Mes-E10.5ctrl_ML7/analysis/midface_Epi_Mes_E10.5.markers.csv")

## Very prominent markers
markers <- FindAllMarkers(
  object = so_midface_Epi_Mes_E10.5,
  only.pos = TRUE,               # Return only positive markers (upregulated in the cluster)
  min.pct = 0.25,                # Gene expressed in at least 25% of cells in either group
  logfc.threshold = 01         # Minimum log fold change
)

# Save the marker list to a CSV file
write.csv(markers, file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/MF_Epi+Mes-E10.5ctrl_ML7/analysis/midface_Epi_Mes_E10.5.prominent_markers.csv")


### Overlaps of expression
# Define pairs of interest
gene_pairs <- list(
  c("Bmpr1a", "Bmp4"),
  c("Bmpr1b", "Bmp4"),
  c("Bmpr2",  "Bmp4"),
  c("Smad1",  "Bmp4"),
  c("Smad2",  "Bmp4"),
  c("Smad4",  "Bmp4"),
  c("Smad5",  "Bmp4"),
  c("Id1",    "Bmp4"),
  c("Id2",    "Bmp4"),
  c("Id3",    "Bmp4")
)

# Output file
outFile <- paste0(output_folder,
                  "/midface_Epi_Mes_E10.5.UMAP.Bmp4_receptors+effectors.orig.ident.pdf")
pdf(outFile, width = 20, height = 5)

for (pair in gene_pairs) {
  if (all(pair %in% rownames(so_midface_Epi_Mes_E10.5))) {
    p <- FeaturePlot(
      so_midface_Epi_Mes_E10.5,
      features = pair,
      reduction = "umap",
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


### Plots for genes of cell adhesion molecules

# Alpha (ITGA) subunits (18)
alpha.integrins <- c(
  "Itga1","Itga2","Itga3","Itga4","Itga5","Itga6","Itga7",
  "Itga8","Itga9","Itga10","Itga11",
  "Itga2b","Itgav","Itgae","Itgal","Itgam","Itgax","Itgad"
)

# Beta (ITGB) subunits (8)
beta.integrins <- c(
  "Itgb1","Itgb2","Itgb3","Itgb4","Itgb5","Itgb6","Itgb7","Itgb8"
)

goi <-  c(alpha.integrins,
          beta.integrins,
          "Thbs1",   # thrombospondin-1
          "Tns1",    # tensin 1
          "Postn",   # periostin (osteoblast specific factor)
          "Fn1",     # fibronectin 1
          "Flrt2",   # fibronectin leucine rich transmambrane protein 2
          "Tgfbi",   # transforming growth factor, beta induced, 68kDA
          "Cdsn",    # corneodesmosin
          "Dsc2"     # desmocollin 2
)

################################ Violin plots  #################################
# Violin Plot for goi
outFile <- paste(output_folder,
                 "/midface_Epi_Mes_E10.5.Violin.goi_cell_adhesion.seurat_clusters.pdf", 
                 sep = "")
pdf(outFile, width = 20, height = 5)
for (gene in goi) {
  if (gene %in% rownames(so_midface_Epi_Mes_E10.5)) {
    p <- VlnPlot(so_midface_Epi_Mes_E10.5,
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
outFile <- paste(output_folder,
                 "/midface_Epi_Mes_E10.5.UMAP.goi_cell_adhesion.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 15, height = 10)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi) {
  if (gene %in% rownames(so_midface_Epi_Mes_E10.5)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_midface_Epi_Mes_E10.5, 
                     features = gene,
                     reduction = "umap",
                     split.by = "orig.ident")
    plot(p)
  } else {
    # Print a message for missing genes (optional)
    message(paste("Gene not found in data: ", gene))
  }
}

dev.off()
