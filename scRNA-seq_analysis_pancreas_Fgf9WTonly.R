# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)

# Set out folder (to store results)
outFolder_pancreas_WT <- "/Users/veralaub/Documents/postdoc/collaboration/Maurizio/Fgf9null_datasets/scRNA-seq/pancreas_Fgf9null/results_WT/"
so_path_pancreas_WT <- paste(outFolder_pancreas_WT, "so_pancreas_WT.rds", sep = "")

# Increase the maximum global size to 2 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# Load pancreas data (WT only)
d_pancreas_Fgf9WT <- Read10X(data.dir ="/Users/veralaub/Documents/postdoc/collaboration/Maurizio/Fgf9null_datasets/scRNA-seq/pancreas_Fgf9null/E14.5_pancreas_Fgf9WT_GSM6434043/")

# Create a Seurat object
so_pancreas_Fgf9WT <- CreateSeuratObject(counts = d_pancreas_Fgf9WT, 
                                         project = "pancreas_E14.5",
                                         min.cells = 3, min.features = 200)

# Add mito fraction (no change here)
so_pancreas_Fgf9WT[["percent.mt"]] <- PercentageFeatureSet(so_pancreas_Fgf9WT,
                                                           pattern = "^mt-")

# QC stats before filtering
p <- RidgePlot(so_pancreas_Fgf9WT, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(outFolder_pancreas_WT, "data.qc.unfiltered.pdf", sep = "")
pdf(out_path, width = 4, height = 4)
plot(p)
dev.off()

# Define thresholds for filtering cells (can be adapted ad gusto)
# Only such cells that pass these criteria are kept for further analysis
# The current threshold filters for rather lowly expressed genes
percent.mt_max <- 5  # maximum percentage of mitochondrial genes (adjust as needed)
nFeature_RNA_min <- 2500  # minimum number of features per cell [1. trial: 300]
nFeature_RNA_max <- 15000 # minimum number of features per cell [1. trial: 2500]
nCount_RNA_min <- 100  # minimum number of RNA counts per cell [1. trial: 100]
nCount_RNA_max <- 50000  # maximum number of RNA counts per cell [1. trial: 5000]

# Subset the Seurat object (filter based on thresholds above)
so_pancreas_Fgf9WT_filtered <- subset(so_pancreas_Fgf9WT,
                                      subset = percent.mt <= percent.mt_max &
                                        nFeature_RNA >= nFeature_RNA_min &
                                        nFeature_RNA <= nFeature_RNA_max &
                                        nCount_RNA >= nCount_RNA_min &
                                        nCount_RNA <= nCount_RNA_max)

# QC stats after filtering
p <- RidgePlot(so_pancreas_Fgf9WT_filtered,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(outFolder_pancreas_WT, "data.qc.filtered.pdf", sep = "")
pdf(out_path, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_pancreas_Fgf9WT_filtered[["RNA"]]), invert = TRUE)

# Exclude mito genes
filtered_counts <- GetAssayData(so_pancreas_Fgf9WT_filtered[["RNA"]], layer = "counts")[keep, ]

# Create Seurat object with filtered genes
# Access counts using GetAssayData instead of directly from counts slot
filtered_counts <- GetAssayData(so_pancreas_Fgf9WT_filtered[["RNA"]], layer = "counts")[keep, ]

# Create a new Seurat object after filtering
so_pancreas_Fgf9WT_filtered <- CreateSeuratObject(counts = filtered_counts, 
                                                  project = "pancreas_E14.5", 
                                                  meta.data = so_pancreas_Fgf9WT_filtered@meta.data)

# Check memory usage before and after merge (due to Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization (SCTransform has not changed)
n_features <- 100
so_pancreas_Fgf9WT_filtered_norm <- SCTransform(so_pancreas_Fgf9WT_filtered,
                                                verbose = TRUE,
                                                variable.features.n = n_features)

# PCA (no changes for PCA)
so_pancreas_Fgf9WT_filtered_norm <- RunPCA(so_pancreas_Fgf9WT_filtered_norm,
                                           verbose = FALSE, npcs = 30)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_pancreas_Fgf9WT_filtered_norm, ndims = 30)
out_path <- paste(outFolder_pancreas_WT, "data.qc.ellbowplot.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

# Store number of principle components in new variable (to be used later)
pca_dim_sel <- 21

# UMAP
so_pancreas_Fgf9WT_filtered_norm <- RunUMAP(so_pancreas_Fgf9WT_filtered_norm,
                                            dims = 1:pca_dim_sel)

# Clustering (Leiden) - Seurat v5 should work similarly
so_pancreas_Fgf9WT_filtered_norm <- FindNeighbors(so_pancreas_Fgf9WT_filtered_norm,
                                                  dims = 1:pca_dim_sel)
so_pancreas_Fgf9WT_filtered_norm <- FindClusters(so_pancreas_Fgf9WT_filtered_norm,
                                                 resolution = 0.4,
                                                 algorithm = 4)

# Save the Seurat object
saveRDS(so_pancreas_Fgf9WT_filtered_norm, file = so_path_pancreas_WT)


#UMAP plot (clusters)
p <- DimPlot(object = so_pancreas_Fgf9WT_filtered_norm,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
outFile_pancreas_WT <- paste(outFolder_pancreas_WT, "/UMAP.clusters.pdf", sep = "")
pdf(outFile_pancreas_WT, width = 7, height = 5)
plot(p)
dev.off()

#GOI: UMAP overlay
goi <- c("Tlx1", "Barx1", "Fgf9", "Fgf10", "Fgfr1", "Fgfr2", "Fgfr3", "Pdx1", "Prrx1",
         "Pbx1", "Pbx2", "Pbx3", "Pknox1", "Pknox2", "Meis1", "Meis2", "Meis3",
         "Cdk1", "Mki67", "Cdkn1c")
outFile_pancreas_WT <- paste(outFolder_pancreas_WT, "/UMAP.goi.pdf", sep = "")
pdf(outFile_pancreas_WT, width = 7, height = 5)
for (gene in goi) {
  p <- FeaturePlot(so_pancreas_Fgf9WT_filtered_norm, gene)
  plot(p)
}
dev.off()

#GOI: vioplot by cluster
outFile_pancreas_WT <- paste(outFolder_pancreas_WT, "/vioplot.goi.pdf", sep = "")
pdf(outFile_pancreas_WT, width = 7, height = 5)
for (gene in goi) {
  p <- VlnPlot(so_pancreas_Fgf9WT_filtered_norm, 
               features = gene, 
               pt.size = 0, 
               group.by = "seurat_clusters") +
    stat_summary(fun = median, geom='point', size = 10, colour = "black", shape = 95) +
    theme(legend.position = 'none')
  plot(p)
}
dev.off()

# Multipannel UMAP of goi
p <- FeaturePlot(so_pancreas_Fgf9WT_filtered_norm, features = goi,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
outFile_pancreas_WT <- paste(outFolder_pancreas_WT, "/UMAP.goi.multipanel.pdf", sep = "")
pdf(outFile_pancreas_WT, width = 25, height = 20)
plot(p)
dev.off()

## UMAP Coexpression
# Define a list of gene sets (co-expressed genes)
gene_set_1 <- list(Coexpression = c("Tlx1",
                                    "Barx1"))

# Add module scores to the Seurat object
data <- AddModuleScore(so_pancreas_Fgf9WT_filtered_norm, 
                       features = gene_set_1, name = "CoexpressionScore")

# Visualize the module score in UMAP
p <- FeaturePlot(data, features = "CoexpressionScore1",
                 pt.size = 0.5) +
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = median(data$CoexpressionScore1)) +
  labs(title = "Co-expression in E10.5 hindlimb scRNA-seq", 
       x = "UMAP 1", 
       y = "UMAP 2", 
       color = "Co-expression\n(Tlx1, Barx1)") +  # Custom legend title
  theme_minimal()
outFile_pancreas_WT <- paste(outFolder_pancreas_WT, "/UMAP.coexpression.pdf", sep = "")
pdf(outFile_pancreas_WT, width = 7, height = 5)
plot(p)
dev.off()

## Analysis of Tlx1 expressing cells
# Subset the Seurat object for Tlx1-expressing cells
tlx1_cells <- subset(so_pancreas_Fgf9WT_filtered_norm, subset = Tlx1 > 0)

# Find markers for Tlx1-expressing cells compared to the rest of the cells
markers_tlx1_vs_rest <- FindMarkers(
  so_pancreas_Fgf9WT_filtered_norm, 
  ident.1 = WhichCells(so_pancreas_Fgf9WT_filtered_norm, expression = Tlx1 > 0), 
  only.pos = TRUE
)

# View the top expressed genes in Tlx1+ cells
head(markers_tlx1_vs_rest)

# Save the marker list to a CSV file
write.csv(markers_tlx1_vs_rest, file = "/Users/veralaub/Documents/postdoc/collaboration/Maurizio/Fgf9null_datasets/scRNA-seq/pancreas_Fgf9null/results_WT/markers_tlx1+cells.csv", row.names = TRUE)

# Feature plot for some top genes (you can adjust the number)
top_genes <- rownames(markers_tlx1_vs_rest)[1:100]
p <- FeaturePlot(tlx1_cells, features = top_genes)
outFile_pancreas_WT <- paste(outFolder_pancreas_WT, "/UMAP.top10genes.Tlx1+cells.pdf", sep = "")
pdf(outFile_pancreas_WT, width = 15, height = 10)
plot(p)
dev.off()

# Feature Plot
goi <- top_genes
outFile_pancreas_WT <- paste(outFolder_pancreas_WT, "/UMAP.goi.tlx1+cells.pdf", sep = "")
pdf(outFile_pancreas_WT, width = 7, height = 5)
for (gene in goi) {
  p <- FeaturePlot(so_pancreas_Fgf9WT_filtered_norm, gene)
  plot(p)
}
dev.off()


############## ANALYSIS WITH WT DATA ONLY
# Load data
so_Fgf9WT_processed <- readRDS("/Users/veralaub/Documents/postdoc/collaboration/Maurizio/Fgf9null_datasets/scRNA-seq/pancreas_Fgf9null/results_WT/so_pancreas_WT.rds")
outFolder_pancreas_WT <- "/Users/veralaub/Documents/postdoc/collaboration/Maurizio/Fgf9null_datasets/scRNA-seq/pancreas_Fgf9null/results_WT/"

## Test for correlated expression of Tlx1 and Barx1 in UMAP
# Step 1: Calculate gene correlation
# Extract expression data for the two genes
gene_1_exprWT <- FetchData(so_Fgf9WT_processed, vars = "Tlx1")
gene_2_exprWT <- FetchData(so_Fgf9WT_processed, vars = "Barx1")

# Calculate the correlation between the two genes
correlation <- cor(gene_1_exprWT, gene_2_exprWT)

## Testing whether the expression of two genes, like Tlx1 and Barx1, is mutually exclusive 
# 1. Define gene expression thresholds
# Extract expression data for Tlx1 and Barx1
Tlx1_exprWT <- FetchData(so_Fgf9WT_processed, vars = "Tlx1")
Barx1_exprWT <- FetchData(so_Fgf9WT_processed, vars = "Barx1")

# Define a threshold for expression (e.g., greater than 0 means expressed)
Tlx1_expressedWT <- Tlx1_exprWT > 0
Barx1_expressedWT <- Barx1_exprWT > 0

# 2. Identify cells where only one gene is expressed
# Identify cells where Tlx1 is expressed but Barx1 is not, and vice versa
mutually_exclusive_Tlx1_WT <- Tlx1_expressedWT & !Barx1_expressedWT
mutually_exclusive_Barx1_WT <- Barx1_expressedWT & !Tlx1_expressedWT
both_expressed_WT <- Tlx1_expressedWT & Barx1_expressedWT

# Count the number of mutually exclusive cells for each gene
mutually_exclusive_Tlx1_cellsWT <- sum(mutually_exclusive_Tlx1_WT)
mutually_exclusive_Barx1_cellsWT <- sum(mutually_exclusive_Barx1_WT)

# Count the number of cells where both genes are expressed
both_expressed_cellsWT <- sum(both_expressed_WT)

# Print the results
cat("Number of cells (WT only) where Tlx1 is expressed but Barx1 is not: ", mutually_exclusive_Tlx1_cellsWT, "\n")
cat("Number of cells (WT only) where Barx1 is expressed but Tlx1 is not: ", mutually_exclusive_Barx1_cellsWT, "\n")
cat("Number of cells (WT only) where both Tlx1 and Barx1 are expressed: ", both_expressed_cellsWT, "\n")

# 3. Statistical test for mutual exclusivity
# Create a contingency table for the co-expression of Tlx1 and Barx1
contingency_table_WT <- table(Tlx1_expressedWT, Barx1_expressedWT)

# Perform Fisher's Exact Test (for small numbers) or Chi-square test
# Fisher's exact test is recommended for small sample sizes (less than 5 in any cell)
fisher_test_result_WT <- fisher.test(contingency_table_WT)

# Print the p-value from the Fisher's test
cat("P-value for mutual exclusivity (Fisher's Exact Test): ", fisher_test_result_WT$p.value, "\n")

# 4. Visualizing the results
# Create a new metadata column to label cells as mutually exclusive for each gene
so_Fgf9WT_processed$MutualExclusive_Tlx1 <- mutually_exclusive_Tlx1_WT
so_Fgf9WT_processed$MutualExclusive_Barx1 <- mutually_exclusive_Barx1_WT
so_Fgf9WT_processed$Both_Expressed <- both_expressed_WT

# Visualize mutual exclusivity of Tlx1 and Barx1 using UMAP
p <- FeaturePlot(so_Fgf9WT_processed, 
                 features = c("MutualExclusive_Tlx1", "MutualExclusive_Barx1", "Both_Expressed"))
outFile_pancreas_WT <- paste(outFolder_pancreas_WT, "/UMAP.mutual.exclusive.pdf", sep = "")
pdf(outFile_pancreas_WT, width = 10, height = 7)
plot(p)
dev.off()

# Visualize mutual exclusivity with a Venn diagram (with two overlapping circles)
library(grid)
library(futile.logger)
library(VennDiagram)

# Define the output file path for the Venn diagram
outFile_venn_WT <- paste(outFolder_pancreas_WT, "/VennDiagram_Tlx1_Barx1_expressed.pdf", sep = "")

# Open a PDF device to save the plot
pdf(outFile_venn_WT, width = 10, height = 5)

# Create the Venn diagram with two circles
venn.plot <- venn.diagram(
  x = list(
    "Tlx1 expressed" = which(Tlx1_expressedWT),
    "Barx1 expressed" = which(Barx1_expressedWT)
  ),
  category.names = c("Tlx1 expressed", "Barx1 expressed"),
  filename = NULL,  # We are using grid.draw() to plot, so no need for a file name here
  output = TRUE,
  lwd = 2,  # Line width of the circles
  fill = c("red", "blue"),  # Fill color for the circles
  alpha = c(0.5, 0.5),  # Transparency of the circles
  cex = 1.5,  # Text size for main title
  cat.cex = 1.5,  # Category text size
  cat.pos = 0,  # Category text position (0 is top-center)
  main = "Venn Diagram of Gene Expression, pancreas E14.5 WT",  # Main title
  family = "sans",  # Use a generic sans-serif font (usually Helvetica or Arial)
  cat.fontface = 1,  # Regular font style for category names
  fontface = 1  # Regular font style for main title
)

# Plot the Venn diagram
grid.draw(venn.plot)

# Close the PDF device to save the plot
dev.off()


########################### ANALYSIS WORK IN PROGRESS #########################
## Markers by cluster
# 12-16-2024: Currently working on this (not finished)
outFile_pancreas_WT <- paste(outFolder_pancreas_WT, "markers.txt", sep = "")
so_pancreas_Fgf9WT_filtered_norm <- PrepSCTFindMarkers(so_pancreas_Fgf9WT_filtered_norm,
                                                       assay = "SCT", 
                                                       verbose = TRUE)
markers <- FindAllMarkers(so_pancreas_Fgf9WT_filtered_norm,
                          min.pct = 0.1,
                          test.use = "wilcox")

# View markers for all clusters
marker_table <- table(markers$cluster)  # Shows the number of markers for each cluster

# Show markers for the first few clusters
# Assuming markers is already calculated using FindAllMarkers
cluster_ids <- unique(so_pancreas_Fgf9WT_filtered_norm$seurat_clusters)   # Get unique cluster identities
num_clusters <- length(cluster_ids)   # Count the number of unique clusters

# Loop over each cluster to extract and print top 10 markers
for (cluster in cluster_ids) {
  # Extract the top 20 markers for this cluster
  top_markers <- head(markers[markers$cluster == cluster, ], 20)  # Get top 10 markers for the current cluster
  # Print the top 20 markers for the current cluster
  cat("Top 20 markers for cluster ", cluster, " are: \n", sep = "")
  # Print the gene names (marker genes) for the current cluster
  print(top_markers$gene)   # Assuming 'gene' is the column containing marker gene names
  cat("\n")  # Add a line break between clusters
}


# UMAP plot (clusters with markers)
# NEEDS TO BE UPDATED
p <- DimPlot(object = so_pancreas_Fgf9WT_filtered_norm,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
outFile_pancreas_WT <- paste(outFolder_pancreas_WT, "/UMAP.clusters.markers.pdf", sep = "")
pdf(outFile_pancreas_WT, width = 7, height = 5)
plot(p)
dev.off()
