## scRNA-seq data analysis of hindlimb E10.5 with Seurat
# author: Vera Laub
# last edited: 2025-01-16

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
out_folder <- "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_E10.5_hindlimb/results/"

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
out_path <- paste(out_folder, "hindlimb_E10.5.qc.unfiltered.pdf", sep = "")
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
out_file <- paste(out_folder, "hindlimb_E10.5.data.qc.filtered.pdf", sep = "")
pdf(out_file, width = 4, height = 4)
plot(p)
dev.off()

# Exclude mitochondrial genes 
keep <- grep("^mt-", rownames(so_hindlimb_E10.5[["RNA"]]), invert = TRUE)
so_hindlimb_E10.5 <- CreateSeuratObject(counts = GetAssayData(so_hindlimb_E10.5[["RNA"]], layer = "counts")[keep, ], 
                                        project = "pancreas_E14.5", 
                                        meta.data = so_hindlimb_E10.5@meta.data)

# Check memory usage before and after merge (to monitor Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Normalization with SCTransform
n_features <- 3000
so_hindlimb_E10.5 <- SCTransform(so_hindlimb_E10.5,
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
pca_dim_sel <- 60                                      # Number of features selection by elbow method

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_hindlimb_E10.5, ndims = 100)
out_file <- paste(out_folder, "hindlimb_E10.5.qc.ellbowplot.pdf", sep = "")
pdf(out_file, width = 5, height = 5)
plot(p)
dev.off()

#UMAP
so_hindlimb_E10.5 <- RunUMAP(so_hindlimb_E10.5, dims = 1:pca_dim_sel)

# Clustering of cells (Leiden algorithm)                                        # Seurat uses graph-based approach to cluster cells
so_hindlimb_E10.5 <- FindNeighbors(so_hindlimb_E10.5, dims = 1:pca_dim_sel)     # "Cell clusters were identified by constructing a shared nearest neighbor graph followed by a modularity optimization-based clustering algorithm (Leiden algorithm) using the top 60 principal components as determined by PCA." (Losa et al., 2023)
so_hindlimb_E10.5 <- FindClusters(so_hindlimb_E10.5, resolution = 0.8)          # "Clustering was performed at multiple resolutions between 0.2 and 2, and optimal resolution was determined empirically based on the expression of known population markers (resolution = 0.8)." (Losa et al., 2023)

# Save Seurat object
saveRDS(so_hindlimb_E10.5, file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_E10.5_hindlimb/results/so_hindlimb_E10.5.rds")

# Read Seurat object (for post-hoc analysis)
so_hindlimb_E10.5 <- readRDS(file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_E10.5_hindlimb/results/so_hindlimb_E10.5.rds")

# Markers by cluster
out_file <- paste(out_folder, "/hindlimb_E10.5.markers.txt", sep = "")
markers <- FindAllMarkers(so_hindlimb_E10.5, 
                          min.pct = 0.1, 
                          logfc.threshold = 0.25, 
                          test.use = "wilcox")
markers_tib <- markers %>% as_tibble(rownames = "Symbol") %>% select(-Symbol) %>% select(gene, everything()) %>% mutate(avg_log2FC = round(avg_log2FC, 4))
write_tsv(markers_tib, file = out_file)

### Visualization
# UMAP: Visualize the S and G2M scores
p <- FeaturePlot(so_hindlimb_E10.5, 
                 features = c("S.Score", "G2M.Score"), 
                 reduction = "umap")
out_file <- paste(out_folder, "/hindlimb_E10.5.UMAP.CellCycleRegressed.pdf", sep = "")
pdf(out_file, width = 15, height = 10)
plot(p)
dev.off()

# UMAP plot (clusters)
p <- DimPlot(object = so_hindlimb_E10.5,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
out_file <- paste(out_folder, "/hindlimb_E10.5.UMAP.clusters.pdf", sep = "")
pdf(out_file, width = 15, height = 10)
plot(p)
dev.off()

# UMAP overlay of goi
goi <- c("Pbx1", "Pbx2", "Hand2",                                                  # Pbx1/2, Hand2
         "Prrx1", "Msx1", "Msx2", "Tfap2c", "Twist1", "Twist2",                    # hindlimb fate determinants (literature)
         "Hoxa1", "Hoxa2", "Hoxa3", "Hoxa4", "Hoxa5", "Hoxa6", "Hoxa7", "Hoxa9",   # Hoxa cluster
         "Hoxa10", "Hoxa11", "Hoxa13")
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

# Vioplot (by cluster) of goi
out_file <- paste(out_folder, "/hindlimb_E10.5.vioplot.goi.pdf", sep = "")
pdf(out_file, width = 15, height = 10)
for (gene in goi) {
  p <- VlnPlot(so_hindlimb_E10.5, 
               features = gene, 
               pt.size = 0, 
               group.by = "seurat_clusters") +
    stat_summary(fun = median, geom='point', size = 10, colour = "black", shape = 95) +
    theme(legend.position = 'none')
  plot(p)
}
dev.off()



############################# ADDITIONAL ANALYSIS ##############################

## Co-expression of Pbx1 and Hand2
# Define a list of gene sets (co-expressed genes)
gene_set_1 <- list(Coexpression = c("Pbx1",
                                  "Hand2"))

# Add module scores to the Seurat object
so_hindlimb_E10.5 <- AddModuleScore(so_hindlimb_E10.5, features = gene_set_1, name = "CoexpressionScore")

# Visualize the Coexpression score in UMAP
p <- FeaturePlot(so_hindlimb_E10.5, features = "CoexpressionScore1", 
                 pt.size = 0.5) + 
  scale_color_gradient2(low = "red", mid = "white", high = "blue", 
                        midpoint = median(so_hindlimb_E10.5$CoexpressionScore1)) +
  labs(title = "Co-expression in hindlimb E10.5 (scRNA-seq)", 
       x = "UMAP 1", 
       y = "UMAP 2", 
       color = "Co-expression\n(Pbx1, Hand2)") +  # Custom legend title
  theme_minimal()
out_file <- paste(out_folder, "/hindlimb_E10.5.UMAP.Pbx1+Hand2.co-expression.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
plot(p)
dev.off()


## Anti-/correlated expression of Pbx1 and Hand2
# Step 1: Calculate gene correlation
# Extract expression data for the two genes
Pbx1_expr <- FetchData(so_hindlimb_E10.5, vars = "Pbx1")
Hand2_expr <- FetchData(so_hindlimb_E10.5, vars = "Hand2")

# Calculate the correlation between the two genes
correlation <- cor(Pbx1_expr, Hand2_expr)
cat("Correlation between Pbx1 and Hand2: ", correlation, "\n")

## Testing whether the expression of Pbx1 and Hand2 is mutually exclusive 
# 1. Define gene expression thresholds
# Extract expression data for Pbx1 and Hand2
Pbx1_expressed <- Pbx1_expr > 0
Hand2_expressed <- Hand2_expr > 0

# 2. Identify cells where only one gene is expressed
# Identify cells where Pbx1 is expressed but Hand2 is not, and vice versa
mutually_exclusive_Pbx1 <- Pbx1_expressed & !Hand2_expressed
mutually_exclusive_Hand2 <- Hand2_expressed & !Pbx1_expressed
both_expressed <- Pbx1_expressed & Hand2_expressed

# Count the number of mutually exclusive cells for each gene
mutually_exclusive_Pbx1_cells <- sum(mutually_exclusive_Pbx1)
mutually_exclusive_Hand2_cells <- sum(mutually_exclusive_Hand2)

# Count the number of cells where both genes are expressed
both_expressed_cells <- sum(both_expressed)

# Print the results
cat("Number of cells where Pbx1 is expressed but Hand2 is not: ", mutually_exclusive_Pbx1_cells, "\n")
cat("Number of cells where Hand2 is expressed but Pbx1 is not: ", mutually_exclusive_Hand2_cells, "\n")
cat("Number of cells where both Pbx1 and Hand2 are expressed: ", both_expressed_cells, "\n")

# 3. Statistical test for mutual exclusivity
# Create a contingency table for the co-expression of Pbx1 and Hand2
contingency_table <- table(Pbx1_expressed, Hand2_expressed)

# Perform Fisher's Exact Test (for small numbers) or Chi-square test
# Fisher's exact test is recommended for small sample sizes (less than 5 in any cell)
fisher_test_result <- fisher.test(contingency_table)

# Print the p-value from the Fisher's test
cat("P-value for mutual exclusivity (Fisher's Exact Test): ", fisher_test_result$p.value, "\n")

# 4. Visualizing the results
# Create a new metadata column to label cells as mutually exclusive for each gene
so_hindlimb_E10.5$MutualExclusive_Pbx1 <- mutually_exclusive_Pbx1
so_hindlimb_E10.5$MutualExclusive_Hand2 <- mutually_exclusive_Hand2
so_hindlimb_E10.5$Both_Expressed <- both_expressed

# Visualize mutual exclusivity of Pbx1 and Hand2 using UMAP
p <- FeaturePlot(so_hindlimb_E10.5, 
                 features = c("MutualExclusive_Pbx1", "MutualExclusive_Hand2", "Both_Expressed"))
out_file <- paste(out_folder, "/hindlimb_E10.5.UMAP.mutual.exclusive.Pbx1_Hand2.pdf", sep = "")
pdf(out_file, width = 10, height = 7)
plot(p)
dev.off()

# Visualize mutual exclusivity with a Venn diagram (with two overlapping circles)
library(grid)
library(futile.logger)
library(VennDiagram)

# Define the output file path for the Venn diagram
out_file <- paste(out_folder, "/hindlimb_E10.5.VennDiagram.Pbx1+Hand2.pdf", sep = "")

# Open a PDF device to save the plot
pdf(out_file, width = 10, height = 5)

# Create the Venn diagram with two circles
venn.plot <- venn.diagram(
  x = list(
    "Pbx1" = which(Pbx1_expressed),
    "Hand2" = which(Hand2_expressed)
  ),
  category.names = c("Pbx1", "Hand2"),
  filename = NULL,  # We are using grid.draw() to plot, so no need for a file name here
  output = TRUE,
  lwd = 2,  # Line width of the circles
  fill = c("red", "blue"),  # Fill color for the circles
  alpha = c(0.5, 0.5),  # Transparency of the circles
  cex = 1.5,  # Text size for main title
  cat.cex = 1.5,  # Category text size
  cat.pos = 0,  # Category text position (0 is top-center)
  main = "Venn Diagram of Gene Expression, Hindlimb E10.5",  # Main title
  family = "sans",  # Use a generic sans-serif font (usually Helvetica or Arial)
  cat.fontface = 1,  # Regular font style for category names
  fontface = 1  # Regular font style for main title
)

# Plot the Venn diagram
grid.draw(venn.plot)

# Close the PDF device to save the plot
dev.off()


##### Co-expression of multiple factors
### 1. Co-expression of all genes
# Define a list of gene sets (co-expressed genes)
gene_set_all <- list(Coexpression = c("Pbx1", "Pbx2", "Hand1", "Hand2", "Tcf3", "Tcf4", "Tcf12"))

# Add module scores to the Seurat object for the entire set
so_hindlimb_E10.5 <- AddModuleScore(so_hindlimb_E10.5, 
                                    features = gene_set_all, 
                                    name = "CoexpressionScore")

# Visualize the Coexpression score in UMAP for the entire gene set
p <- FeaturePlot(so_hindlimb_E10.5, features = "CoexpressionScore1", 
                 pt.size = 0.5) + 
  scale_color_gradient2(low = "red", mid = "white", high = "blue", 
                        midpoint = median(so_hindlimb_E10.5$CoexpressionScore1)) +
  labs(title = "Co-expression in hindlimb E10.5 (scRNA-seq)", 
       x = "UMAP 1", 
       y = "UMAP 2", 
       color = "Co-expression\n(Pbx1, Pbx2, Hand1, Hand2, Tcf3, Tcf4, Tcf12)") +  # Updated legend title
  theme_minimal()

out_file <- paste(out_folder, "/hindlimb_E10.5.UMAP.Coexpression_Pbx1+2_Hand1+2_Tcf3+4+12.pdf", sep = "")
pdf(out_file, width = 7, height = 5)
plot(p)
dev.off()

### 2. Calculate Correlation Between All Genes
# Extract expression data for all genes of interest
genes_of_interest <- c("Pbx1", "Pbx2", "Hand1", "Hand2", "Tcf3", "Tcf4", "Tcf12")
expression_data <- FetchData(so_hindlimb_E10.5, vars = genes_of_interest)

# Calculate pairwise correlations
cor_matrix <- cor(expression_data)
print("Correlation matrix:")
print(cor_matrix)


### Testing Mutual Exclusivity for All Genes
# Define expression threshold (e.g., expressed if > 0)
expression_threshold <- 0

# Initialize a list to store mutual exclusivity results
mutually_exclusive_results <- list()

# Loop through all pairs of genes
for (i in 1:length(genes_of_interest)) {
  for (j in (i+1):length(genes_of_interest)) {
    gene1 <- genes_of_interest[i]
    gene2 <- genes_of_interest[j]
    
    # Extract expression data for the pair of genes
    gene1_expr <- FetchData(so_hindlimb_E10.5, vars = gene1)
    gene2_expr <- FetchData(so_hindlimb_E10.5, vars = gene2)
    
    # Define expression (0 or 1, based on threshold)
    gene1_expressed <- gene1_expr > expression_threshold
    gene2_expressed <- gene2_expr > expression_threshold
    
    # Calculate mutual exclusivity
    mutually_exclusive_gene1 <- gene1_expressed & !gene2_expressed
    mutually_exclusive_gene2 <- gene2_expressed & !gene1_expressed
    both_expressed <- gene1_expressed & gene2_expressed
    
    # Store results
    mutually_exclusive_results[[paste(gene1, gene2, sep = "_")]] <- list(
      mutually_exclusive_gene1 = sum(mutually_exclusive_gene1),
      mutually_exclusive_gene2 = sum(mutually_exclusive_gene2),
      both_expressed = sum(both_expressed)
    )
    
    # Create a contingency table for the pair
    contingency_table <- table(gene1_expressed, gene2_expressed)
    
    # Perform Fisher's exact test
    fisher_test_result <- fisher.test(contingency_table)
    
    # Print the p-value for mutual exclusivity
    cat("P-value for mutual exclusivity between", gene1, "and", gene2, ": ", fisher_test_result$p.value, "\n")
  }
}

# Print the results for mutual exclusivity
print(mutually_exclusive_results)


### Visualization of Mutual Exclusivity for All Pairs
# Loop through and plot mutual exclusivity for all gene pairs
for (i in 1:length(genes_of_interest)) {
  for (j in (i+1):length(genes_of_interest)) {
    gene1 <- genes_of_interest[i]
    gene2 <- genes_of_interest[j]
    
    # Create new metadata columns for mutual exclusivity
    so_hindlimb_E10.5[[paste("MutualExclusive_", gene1, gene2, sep = "")]] <- 
      (FetchData(so_hindlimb_E10.5, vars = gene1) > expression_threshold) & 
      !(FetchData(so_hindlimb_E10.5, vars = gene2) > expression_threshold)
    
    # Visualize mutual exclusivity using UMAP
    p <- FeaturePlot(so_hindlimb_E10.5, features = paste("MutualExclusive_", gene1, gene2, sep = ""))
    out_file <- paste(out_folder, "/hindlimb_E10.5.UMAP.mutual.exclusive_", gene1, "_", gene2, ".pdf", sep = "")
    pdf(out_file, width = 7, height = 5)
    plot(p)
    dev.off()
  }
}


### Venn Diagrams for Pairwise Gene Expression
# Create a Venn diagram for each gene pair
for (i in 1:length(genes_of_interest)) {
  for (j in (i+1):length(genes_of_interest)) {
    gene1 <- genes_of_interest[i]
    gene2 <- genes_of_interest[j]
    
    # Define the output file for the Venn diagram
    out_file <- paste(out_folder, "/hindlimb_E10.5.VennDiagram.Coexpression_", gene1, "_", gene2, ".pdf", sep = "")
    
    # Create the Venn diagram for the pair of genes
    venn.plot <- venn.diagram(
      x = list(
        gene1 = which(FetchData(so_hindlimb_E10.5, vars = gene1) > expression_threshold),
        gene2 = which(FetchData(so_hindlimb_E10.5, vars = gene2) > expression_threshold)
      ),
      category.names = c(gene1, gene2),
      filename = NULL,  # Plot it without saving to a file
      output = TRUE
    )
    
    # Save the plot to the specified file
    pdf(out_file, width = 10, height = 5)
    grid.draw(venn.plot)
    dev.off()
  }
}