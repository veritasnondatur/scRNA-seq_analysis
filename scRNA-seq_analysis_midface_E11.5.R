## scRNA-seq data analysis of midface E11.5 with Seurat
# author: Vera Laub
# last edited: 2025-01-17

# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)

# to unzip downloaded files with terminal:
# cd ~/filepath
# tar -xvf GSE125416_RAW.tar 
# gunzip *.gz 

# Set out folder (to store results)
out_folder <- "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/MF-E11.5_rep1_ML8_071719/results/"

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
out_path <- paste(out_folder, "midface_E11.5.qc.unfiltered.pdf", sep = "")
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
out_file <- paste(out_folder, "midface_E11.5.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_midface_E11.5[["RNA"]]), invert = TRUE)
so_midface_E11.5 <- CreateSeuratObject(counts = GetAssayData(so_midface_E11.5[["RNA"]], layer = "counts")[keep, ], 
                                       project = "pancreas_E14.5", 
                                       meta.data = so_midface_E11.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 3000
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
pca_dim_sel <- 60                                      # Number of features selection by elbow method

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_midface_E11.5, ndims = 100)
out_file <- paste(out_folder, "midface_E11.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

#UMAP
so_midface_E11.5 <- RunUMAP(so_midface_E11.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_midface_E11.5 <- FindNeighbors(so_midface_E11.5, dims = 1:pca_dim_sel)         # "Cell clusters were identified by constructing a shared nearest neighbor graph followed by a modularity optimization-based clustering algorithm (Leiden algorithm) using the top 60 principal components as determined by PCA." (Losa et al., 2023)
so_midface_E11.5 <- FindClusters(so_midface_E11.5, resolution = 0.8)              # "Clustering was performed at multiple resolutions between 0.2 and 2, and optimal resolution was determined empirically based on the expression of known population markers (resolution = 0.8)." (Losa et al., 2023)

# Save Seurat object
saveRDS(so_midface_E11.5, file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/MF-E11.5_rep1_ML8_071719/results/so_midface_E11.5.rds")

# Read Seurat object (for post-hoc analysis)
so_midface_E11.5 <- readRDS(file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/MF-E11.5_rep1_ML8_071719/results/so_midface_E11.5.rds")

# Markers by cluster
out_file <- paste(out_folder, "/midface_E11.5.markers.txt", sep = "")
markers <- FindAllMarkers(so_midface_E11.5, 
                          min.pct = 0.1, 
                          logfc.threshold = 0.25, 
                          test.use = "wilcox")
markers_tib <- markers %>% as_tibble(rownames = "Symbol") %>% select(-Symbol) %>% select(gene, everything()) %>% mutate(avg_log2FC = round(avg_log2FC, 4))
write_tsv(markers_tib, file = out_file)

### Visualization
# UMAP: Visualize the S and G2M scores
p <- FeaturePlot(so_midface_E11.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_file <- paste(out_folder, "/midface_E11.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_file, width = 15, height = 10)
plot(p)
dev.off()

# UMAP plot (clusters)
p <- DimPlot(object = so_midface_E11.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_file <- paste(out_folder, "/midface_E11.5.UMAP.clusters.pdf", sep = "")
pdf(out_file, width = 15, height = 10)
plot(p)
dev.off()

# UMAP overlay of goi
goi <- c("Pbx1", "Pbx2", "Pbx3", "Pbx4",                                           # Pbx1-4
         "Cdh2", "Vim", "Fn1", "Snai1", "Snai2", "Twist1")                         # midface fate determinants (literature)
out_file <- paste(out_folder, "/midface_E11.5.UMAP.goi.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
for (gene in goi) {
  p <- FeaturePlot(so_midface_E11.5, gene)
  plot(p)
}
dev.off()

# Multipannel UMAP of goi
p <- FeaturePlot(so_midface_E11.5, features = goi,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value
out_file <- paste(out_folder, "/midface_E11.5.UMAP.goi.multipanel.pdf", sep = "")
pdf(out_file, width = 25, height = 20)
plot(p)
dev.off()

# Vioplot (by cluster) of goi
out_file <- paste(out_folder, "/midface_E11.5.vioplot.goi.pdf", sep = "")
pdf(out_file, width = 15, height = 10)
for (gene in goi) {
  p <- VlnPlot(so_midface_E11.5, 
               features = gene, 
               pt.size = 0, 
               group.by = "seurat_clusters") +
    stat_summary(fun = median, geom='point', size = 10, colour = "black", shape = 95) +
    theme(legend.position = 'none')
  plot(p)
}
dev.off()



############################# ADDITIONAL ANALYSIS ##############################

## Co-expression of Pbx1 and Draxin (Hand2 not expressed, therefore tested for Draxin)
# Define a list of gene sets (co-expressed genes)
gene_set_1 <- list(Coexpression = c("Pbx1",
                                    "Draxin"))

# Add module scores to the Seurat object
so_midface_E11.5 <- AddModuleScore(so_midface_E11.5, features = gene_set_1, name = "CoexpressionScore")

# Visualize the Coexpression score in UMAP
p <- FeaturePlot(so_midface_E11.5, features = "CoexpressionScore1", 
                 pt.size = 0.5) + 
  scale_color_gradient2(low = "red", mid = "white", high = "blue", 
                        midpoint = median(so_midface_E11.5$CoexpressionScore1)) +
  labs(title = "Co-expression in midface E11.5 (scRNA-seq)", 
       x = "UMAP 1", 
       y = "UMAP 2", 
       color = "Co-expression\n(Pbx1, Draxin)") +  # Custom legend title
  theme_minimal()
out_file <- paste(out_folder, "/midface_E11.5.UMAP.Pbx1+Draxin.co-expression.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
plot(p)
dev.off()


## Anti-/correlated expression of Pbx1 and Draxin
# Step 1: Calculate gene correlation
# Extract expression data for the two genes
Pbx1_expr <- FetchData(so_midface_E11.5, vars = "Pbx1")
Draxin_expr <- FetchData(so_midface_E11.5, vars = "Draxin")

# Calculate the correlation between the two genes
correlation <- cor(Pbx1_expr, Draxin_expr)
cat("Correlation between Pbx1 and Draxin: ", correlation, "\n")

## Testing whether the expression of Pbx1 and Draxin is mutually exclusive 
# 1. Define gene expression thresholds
# Extract expression data for Pbx1 and Draxin
Pbx1_expressed <- Pbx1_expr > 0
Draxin_expressed <- Draxin_expr > 0

# 2. Identify cells where only one gene is expressed
# Identify cells where Pbx1 is expressed but Draxin is not, and vice versa
mutually_exclusive_Pbx1 <- Pbx1_expressed & !Draxin_expressed
mutually_exclusive_Draxin <- Draxin_expressed & !Pbx1_expressed
both_expressed <- Pbx1_expressed & Draxin_expressed

# Count the number of mutually exclusive cells for each gene
mutually_exclusive_Pbx1_cells <- sum(mutually_exclusive_Pbx1)
mutually_exclusive_Draxin_cells <- sum(mutually_exclusive_Draxin)

# Count the number of cells where both genes are expressed
both_expressed_cells <- sum(both_expressed)

# Print the results
cat("Number of cells where Pbx1 is expressed but Draxin is not: ", mutually_exclusive_Pbx1_cells, "\n")
cat("Number of cells where Draxin is expressed but Pbx1 is not: ", mutually_exclusive_Draxin_cells, "\n")
cat("Number of cells where both Pbx1 and Draxin are expressed: ", both_expressed_cells, "\n")

# 3. Statistical test for mutual exclusivity
# Create a contingency table for the co-expression of Pbx1 and Draxin
contingency_table <- table(Pbx1_expressed, Draxin_expressed)

# Perform Fisher's Exact Test (for small numbers) or Chi-square test
# Fisher's exact test is recommended for small sample sizes (less than 5 in any cell)
fisher_test_result <- fisher.test(contingency_table)

# Print the p-value from the Fisher's test
cat("P-value for mutual exclusivity (Fisher's Exact Test): ", fisher_test_result$p.value, "\n")

# 4. Visualizing the results
# Create a new metadata column to label cells as mutually exclusive for each gene
so_midface_E11.5$MutualExclusive_Pbx1 <- mutually_exclusive_Pbx1
so_midface_E11.5$MutualExclusive_Draxin <- mutually_exclusive_Draxin
so_midface_E11.5$Both_Expressed <- both_expressed

# Visualize mutual exclusivity of Pbx1 and Draxin using UMAP
p <- FeaturePlot(so_midface_E11.5, 
                 features = c("MutualExclusive_Pbx1", "MutualExclusive_Draxin", "Both_Expressed"))
out_file <- paste(out_folder, "/midface_E11.5.UMAP.mutual.exclusive.Pbx1_Draxin.pdf", sep = "")
pdf(out_file, width = 10, height = 7)
plot(p)
dev.off()

# Visualize mutual exclusivity with a Venn diagram (with two overlapping circles)
library(grid)
library(futile.logger)
library(VennDiagram)

# Define the output file path for the Venn diagram
out_file <- paste(out_folder, "/midface_E11.5.VennDiagram.Pbx1+Draxin.pdf", sep = "")

# Open a PDF device to save the plot
pdf(out_file, width = 10, height = 5)

# Create the Venn diagram with two circles
venn.plot <- venn.diagram(
  x = list(
    "Pbx1" = which(Pbx1_expressed),
    "Draxin" = which(Draxin_expressed)
  ),
  category.names = c("Pbx1", "Draxin"),
  filename = NULL,  # We are using grid.draw() to plot, so no need for a file name here
  output = TRUE,
  lwd = 2,  # Line width of the circles
  fill = c("red", "blue"),  # Fill color for the circles
  alpha = c(0.5, 0.5),  # Transparency of the circles
  cex = 1.5,  # Text size for main title
  cat.cex = 1.5,  # Category text size
  cat.pos = 0,  # Category text position (0 is top-center)
  main = "Venn Diagram of Gene Expression, midface E11.5",  # Main title
  family = "sans",  # Use a generic sans-serif font (usually Helvetica or Arial)
  cat.fontface = 1,  # Regular font style for category names
  fontface = 1  # Regular font style for main title
)

# Plot the Venn diagram
grid.draw(venn.plot)

# Close the PDF device to save the plot
dev.off()


########################### OLD EXPLORATORY ANALYSIS (Pbx and Draxin expression)

# Retrieve dataset, pre-analyzed with standard work flow (SWF) from previous steps to visualize batch effects
so_midface_E11.5 <- readRDS(file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_midface/MF-E11.5_rep1_ML8_071719/results/so_midface_E11.5.rds")

# UMAP clusters
DimPlot(so_midface_E11.5, reduction = "umap",label = TRUE)

# Pbx1/2 and Draxin expression
plot <- FeaturePlot(so_midface_E11.5, 
                    features = c("Pbx1", "Pbx2", "Pbx3", "Draxin"),                    
                    cols = c('lightgray', 'blue'),
                    pt.size = 0.01)  # Adjust pt.size to your desired value

# Use patchwork to add a title to the whole plot grid
plot + plot_annotation(title = "scRNA-seq analysis of murine midface epithelium, E11.5")

# Expression of Pbx1/2, Draxin and epithelial markers
plot <- FeaturePlot(so_midface_E11.5, 
                    features = c("Pbx1", "Pbx2", "Pbx3", "Draxin", "Cdh1", "Tjp1"),                    
                    cols = c('lightgray', 'blue'),
                    pt.size = 0.01)  # Adjust pt.size to your desired value

# Use patchwork to add a title to the whole plot grid
plot + plot_annotation(title = "scRNA-seq of murine midface epithelium, E11.5: epithelial markers")


# Expression of Pbx1/2, Draxin and mesenchymal markers
plot <- FeaturePlot(so_midface_E11.5, 
                    features = c("Pbx1", "Pbx2", "Pbx3", "Draxin", "Cdh2", "Vim", "Fn1", "Snai1", "Snai2", "Twist1"),                    
                    cols = c('lightgray', 'blue'),
                    pt.size = 0.01)  # Adjust pt.size to your desired value

# Use patchwork to add a title to the whole plot grid
plot + plot_annotation(title = "scRNA-seq of murine midface epithelium, E11.5: mesenchymal markers")

# Define a list of gene sets (co-expressed genes)
gene_set_1 <- list(Coexpression = c("Pbx1", 
                                    "Draxin"))

# Add module scores to the Seurat object
data <- AddModuleScore(so_midface_E11.5, features = gene_set_1, name = "CoexpressionScore")

# Visualize the module score in UMAP
FeaturePlot(data, features = "CoexpressionScore1", 
            pt.size = 0.5) + 
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = median(data$CoexpressionScore1)) +
  labs(title = "Co-expression in midface E11.5 (scRNA-seq)", 
       x = "UMAP 1", 
       y = "UMAP 2", 
       color = "Co-expression\n(Pbx1, Draxin)") +  # Custom legend title
  theme_minimal()

# Violin plot showing expression of a given target in all clusters
VlnPlot(so_midface_E11.5, 
        features = "Draxin",            # Gene of interest
        group.by = "seurat_clusters", # Group by clusters (default is 'seurat_clusters' if you haven't customized the metadata)
        pt.size = 0.1,                # Adjust point size for individual cells, or set to 0 to hide points
        cols = NULL,                  # Optional: customize colors, NULL uses default palette
        y.max = NULL                  # Optional: limit the Y-axis range
)

