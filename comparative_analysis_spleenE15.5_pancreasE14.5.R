########### COMPARATIVE ANALYSIS OF EMBRYONIC SPLEEN AND PANCREAS ##############

## Analysis of differentially regulated genes in Fgf9null E14.5 pancreas from bulk RNA-seq
# Load necessary libraries
library(readxl)
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)

# Step 1: Load bulk RNA-seq data (Fgf9 null E14.5 pancreas) from the Excel file
data_file <- "~/Documents/postdoc/collaboration/Maurizio/Fgf9null_datasets/bulk_RNA-seq/GSE210574_E13.5_E14.5_Control_Null_Counts.xls"
raw_data <- read_excel(data_file)

# Step 2: Prepare the data for pancreas E14.5
gene_data_pancreasE14.5 <- raw_data %>%
  select(gene_name, gene_id, starts_with("E14"))

# Step 3: Create a DESeq2-compatible format
# Extract count data for E14 samples
count_data_pancreasE14.5 <- gene_data_pancreasE14.5 %>%
  select(starts_with("E14")) %>%
  as.matrix()

# Sample information for E14 pancreas
col_data_pancreasE14.5 <- data.frame(
  condition = c("control", "control", "control", "Fgf9_null", "Fgf9_null"),
  row.names = colnames(count_data_pancreasE14.5)
)

# Step 4: DESeq2 analysis
dds_pancreasE14.5 <- DESeqDataSetFromMatrix(countData = count_data_pancreasE14.5, 
                                            colData = col_data_pancreasE14.5, 
                                            design = ~ condition)

# Run DESeq2
dds_pancreasE14.5 <- DESeq(dds_pancreasE14.5)

# Step 5: Extract results from DESeq2
res_pancreasE14.5 <- results(dds_pancreasE14.5, contrast = c("condition", "Fgf9_null", "control"))

# Convert DESeq2 results to a data frame
res_df_pancreasE14.5 <- as.data.frame(res_pancreasE14.5)
gene_info <- gene_data_pancreasE14.5 %>% select(gene_name, gene_id)
res_df_pancreasE14.5 <- cbind(gene_info, res_df_pancreasE14.5)

# Adjust p-values for multiple comparisons using FDR (False Discovery Rate)
res_df_pancreasE14.5$padj <- p.adjust(res_df_pancreasE14.5$pvalue, method = "BH")

# Step 6: Combine raw counts with DESeq2 results
# Extract raw counts for E14 samples (will be used for comparison)
raw_counts <- gene_data_pancreasE14.5 %>% select(starts_with("E14"))

# Combine DESeq2 results first, then raw counts
final_results <- cbind(res_df_pancreasE14.5, raw_counts)

# Step 7: Filter significantly dysregulated genes (optional)
significant_genes <- final_results %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

# Step 8: Output results
# Save the combined table (DESeq2 results first, then raw counts)
write.csv(final_results, 
          "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/bulkRNA-seq_pancreasE14.5_combined_results.csv", 
          row.names = FALSE)

# Save the filtered significant genes (optional)
write.csv(significant_genes, 
          "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/bulkRNA-seq_pancreasE14.5_significant_genes.csv", 
          row.names = FALSE)

# Filter all upregulated and downregulated genes
# Upregulated genes: log2FoldChange > 1
upregulated_genes <- final_results %>%
  filter(log2FoldChange > 1 & padj < 0.05)

# Downregulated genes: log2FoldChange < -1
downregulated_genes <- final_results %>%
  filter(log2FoldChange < -1 & padj < 0.05)

# Step 8: Output results
# Save the combined table (DESeq2 results first, then raw counts)
write.csv(final_results, 
          "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/bulkRNA-seq_pancreasE14.5_combined_results.csv", 
          row.names = FALSE)

# Save all upregulated genes
write.csv(upregulated_genes, 
          "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/bulkRNA-seq_pancreasE14.5_upregulated_genes.csv", 
          row.names = FALSE)

# Save all downregulated genes
write.csv(downregulated_genes, 
          "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/bulkRNA-seq_pancreasE14.5_downregulated_genes.csv", 
          row.names = FALSE)


# Step 9: Visualize results
# Plot a volcano plot to visualize upregulated vs downregulated genes
# Extract top 10 upregulated and downregulated genes based on log2FoldChange
top_upregulated <- res_df_pancreasE14.5 %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

top_downregulated <- res_df_pancreasE14.5 %>%
  filter(padj < 0.05 & log2FoldChange < 0) %>%
  arrange(log2FoldChange) %>%
  head(10)

# Combine top upregulated and downregulated genes for labeling
top_genes <- bind_rows(
  mutate(top_upregulated, direction = "Upregulated"),
  mutate(top_downregulated, direction = "Downregulated")
)

# Create the volcano plot with gene names displayed for the top 10 upregulated and downregulated genes
volcano_plot <- ggplot(res_df_pancreasE14.5, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha = 0.5) +
  scale_color_manual(
    values = c("black", "red"),
    labels = c("FALSE", "TRUE")  # Custom legend labels
  ) +
  labs(
    title = "Volcano Plot of Differential Expression (E14.5 pancreas bulk RNA-seq)",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value"
  ) +
  geom_text(
    data = top_genes, 
    aes(x = log2FoldChange, y = -log10(padj), label = gene_name),
    vjust = 1, hjust = 1, size = 3, color = "blue", fontface = "italic"
  ) +
  theme_minimal() +  # Use a minimal theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white"),  # White background for the plot area
    plot.background = element_rect(fill = "white"),   # White background for the whole plot
    panel.grid.major = element_line(color = "gray", size = 0.2),  # Optional: light grid lines
    panel.grid.minor = element_blank(),  # Optional: remove minor grid lines
    plot.title = element_text(hjust = 0.5),  # Center the plot title
    legend.position = c(0.8, 0.8),  # Move legend inside the plot at a specific location (top-right)
    legend.background = element_rect(fill = "white", color = "black"),  # Optional: add a border around the legend
    legend.title = element_text(face = "bold", size = 10)  # Title for the legend, optional styling
  ) +
  guides(
    color = guide_legend(title = "Differential expression")  # Add the legend title
  )

# Save the plot to the specified folder
# Define the file path
output_file_path <- "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/bulkRNA-seq_pancreasE14.5_differential_expression_volcano_plot.png"

# Save the plot as a PNG file
ggsave(output_file_path, plot = volcano_plot, width = 17, height = 10, dpi = 300)

# Print the output file path
cat("Volcano plot saved to:", output_file_path, "\n")


################################################################################
# How are differentially regulated genes in Fgf9 null E14.5 pancreas expressed in E15.5 spleen?

# Load E15.5 spleen scRNA-seq data
so_spleenE15.5_processed <- readRDS("/Users/veralaub/Documents/postdoc/collaboration/Maurizio/E15.5_spleen/results_spleenE15.5/so_spleenE15.5.rds")

##### UPREGULATED GENES in E14.5 Fgf9 null pancreas
# Extract the gene names of the upregulated genes as a vector of strings
upregulated_gene_names <- upregulated_genes$gene_name

# Print to check the result
print(upregulated_gene_names)

# Upregulated GOI: UMAP overlay
goi <- upregulated_gene_names
output_folder <- "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/"
outFile <- paste(output_folder, "/UMAP.goi_upregulated_Fgf9null_E14.5pancreas.in_scRNA-seq_spleenE15.5.pdf", sep = "")
pdf(outFile, width = 7, height = 5)

# Loop through each gene and check if it exists in the Seurat object
for (gene in goi) {
  if (gene %in% rownames(so_spleenE15.5_processed)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_spleenE15.5_processed, gene)
    plot(p)
  } else {
    # Print a message for missing genes (optional)
    message(paste("Gene not found in data: ", gene))
  }
}

dev.off()

# Upregulated GOI: UMAP Multipanel of top upregulated genes
# Step 1: Filter genes that have more than 50 reads in all control samples
control_samples <- c("E14_C1", "E14_C2", "E14_C3")  # Replace with actual control sample names
final_results_filtered <- final_results %>%
  filter(rowSums(select(., starts_with("E14_C")) > 50) == length(control_samples))

# Step 2: Select upregulated genes (log2FoldChange > 1) and calculate fold enrichment (2^log2FoldChange)
upregulated_genes <- final_results_filtered %>%
  filter(log2FoldChange > 1 & padj < 0.05)

# Calculate fold enrichment for upregulated genes
upregulated_genes$fold_enrichment <- 2^upregulated_genes$log2FoldChange

# Step 3: Sort by fold enrichment and select top 16 upregulated genes
top16_upregulated_genes <- upregulated_genes %>%
  arrange(desc(fold_enrichment)) %>%
  head(16) %>%
  pull(gene_name)  # Extract the gene names as a vector

# Print to check the selected genes
print(top16_upregulated_genes)

# Step 4: Multipanel UMAP of E15.5 spleen (top 16 upregulated genes in E14.5 pancreas)
goi <- top16_upregulated_genes  # Set the selected top 16 upregulated genes

# Create the UMAP plot
p <- FeaturePlot(so_spleenE15.5_processed, features = goi,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value

# Output file path for the multipanel UMAP plot
output_folder <- "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/"
outFile <- paste(output_folder, "/UMAP.goi_top16_upregulated_Fgf9null_E14.5pancreas.in_scRNA-seq_spleenE15.5.multipanel.pdf", sep = "")

# Save the plot as a PDF
pdf(outFile, width = 30, height = 20)
plot(p)
dev.off()

# Create the violin plot
p <- VlnPlot(so_spleenE15.5_processed, features = goi,
             pt.size = 0.1,  # Adjust point size if needed
             ncol = 4)  # Specify number of columns for the panel layout

# Output file path for the multipanel violin plot
output_folder <- "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/"
outFile <- paste(output_folder, "/ViolinPlot.goi_top16_upregulated_Fgf9null_E14.5pancreas.in_scRNA-seq_spleenE15.5.multipanel.pdf", sep = "")

# Save the plot as a PDF
pdf(outFile, width = 30, height = 20)
plot(p)
dev.off()


##### DOWNREGULATED GENES in E14.5 Fgf9 null pancreas
# Extract the gene names of the downregulated genes as a vector of strings
downregulated_gene_names <- downregulated_genes$gene_name

# Print to check the result
print(downregulated_gene_names)

# Downregulated GOI: UMAP overlay
goi <- downregulated_gene_names
output_folder <- "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/"
outFile <- paste(output_folder, "/UMAP.goi_downregulated_Fgf9null_E14.5pancreas.in_scRNA-seq_spleenE15.5.pdf", sep = "")
pdf(outFile, width = 7, height = 5)

# Loop through each gene and check if it exists in the Seurat object
for (gene in goi) {
  if (gene %in% rownames(so_spleenE15.5_processed)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_spleenE15.5_processed, gene)
    plot(p)
  } else {
    # Print a message for missing genes (optional)
    message(paste("Gene not found in data: ", gene))
  }
}

dev.off()

# Downregulated GOI: UMAP Multipanel of top downregulated genes
# Step 1: Filter genes that have more than 50 reads in all control samples
control_samples <- c("E14_C1", "E14_C2", "E14_C3")  # Replace with actual control sample names
final_results_filtered <- final_results %>%
  filter(rowSums(select(., starts_with("E14_C")) > 50) == length(control_samples))

# Step 2: Select downregulated genes (log2FoldChange < -1) and calculate fold enrichment (2^log2FoldChange)
downregulated_genes <- final_results_filtered %>%
  filter(log2FoldChange < -1 & padj < 0.05)

# Calculate fold enrichment for downregulated genes
downregulated_genes$fold_enrichment <- 2^downregulated_genes$log2FoldChange

# Step 3: Sort by fold enrichment and select top 16 downregulated genes
top16_downregulated_genes <- downregulated_genes %>%
  arrange(desc(fold_enrichment)) %>%
  head(16) %>%
  pull(gene_name)  # Extract the gene names as a vector

# Print to check the selected genes
print(top16_downregulated_genes)

# Step 4: Multipanel UMAP of E15.5 spleen (top 16 downregulated genes in E14.5 pancreas)
goi <- top16_downregulated_genes  # Set the selected top 16 downregulated genes

# Create the UMAP plot
p <- FeaturePlot(so_spleenE15.5_processed, features = goi,
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value

# Output file path for the multipanel UMAP plot
output_folder <- "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/"
outFile <- paste(output_folder, "/UMAP.goi_top16_downregulated_Fgf9null_E14.5pancreas.in_scRNA-seq_spleenE15.5.multipanel.pdf", sep = "")

# Save the plot as a PDF
pdf(outFile, width = 30, height = 20)
plot(p)
dev.off()

# Create the violin plot
p <- VlnPlot(so_spleenE15.5_processed, features = goi,
             pt.size = 0.1,  # Adjust point size if needed
             ncol = 4)  # Specify number of columns for the panel layout

# Output file path for the multipanel violin plot
output_folder <- "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/"
outFile <- paste(output_folder, "/ViolinPlot.goi_top16_downregulated_Fgf9null_E14.5pancreas.in_scRNA-seq_spleenE15.5.multipanel.pdf", sep = "")

# Save the plot as a PDF
pdf(outFile, width = 30, height = 20)
plot(p)
dev.off()


################################################################################
# How are marker genes from Tlx1+ cells in E15.5 spleen expressed in E14.5 pancreas scRNA-seq?

# Load data
so_Fgf9WT_processed <- readRDS("/Users/veralaub/Documents/postdoc/collaboration/Maurizio/Fgf9null_datasets/scRNA-seq/pancreas_Fgf9null/results_WT/so_pancreas_WT.rds")
output_folder <- "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/"

# Retrieve marker genes of Tlx1+ cells in E15.5 spleen from previous analysis
tlx1_spleen_markers <- read.csv("~/Documents/postdoc/collaboration/Maurizio/E15.5_spleen/results_spleenE15.5/markers_tlx1+cells.csv")
top100_tlx1_spleen_markers <- tlx1_spleen_markers %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
  arrange(desc(avg_log2FC)) %>%
  head(100)

print(top100_tlx1_spleen_markers)

goi <- top100_tlx1_spleen_markers[[1]]    # Store gene names of marker genes (Tlx1+ cells in E15.5 spleen) in a vector

# GOI: UMAP overlay (E14.5 pancreas scRNA-seq; Top100 Tlx1+ cell markers from E15.5 spleen, previous analysis)
goi <-  top100_tlx1_spleen_markers[[1]]
outFile <- paste(output_folder, "/UMAP.goi_top100_Tlx1+markers_E15.5spleen.in_scRNA-seq_E14.5pancreas.pdf", sep = "")
pdf(outFile, width = 7, height = 5)

# Loop through each gene and check if it exists in the Seurat object
for (gene in goi) {
  if (gene %in% rownames(so_Fgf9WT_processed)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_Fgf9WT_processed, gene)
    plot(p)
  } else {
    # Print a message for missing genes (optional)
    message(paste("Gene not found in data: ", gene))
  }
}
dev.off()

## Multipanel UMAP (E14.5 pancreas scRNA-seq; Top16 Tlx1+ cell markers from E15.5 spleen, previous analysis)
p <- FeaturePlot(so_Fgf9WT_processed, features = goi[1:19],
                 cols = c('lightgray', 'blue'),
                 pt.size = 0.01)   # Adjust pt.size to your desired value

# Output file path for the multipanel UMAP plot
output_folder <- "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/"
outFile <- paste(output_folder, "/UMAP.goi_top16_Tlx1+markers_E15.5spleen.in_scRNA-seq_E14.5pancreas.multipanel.pdf", sep = "")

# Save the plot as a PDF
pdf(outFile, width = 30, height = 20)
plot(p)
dev.off()

## Multipanel Violin Plot (E14.5 pancreas scRNA-seq; Top16 Tlx1+ cell markers from E15.5 spleen, previous analysis)
# Create the violin plot
p <- VlnPlot(so_Fgf9WT_processed, features = goi[1:19],
             pt.size = 0.1,  # Adjust point size if needed
             ncol = 4)  # Specify number of columns for the panel layout

# Output file path for the multipanel violin plot
output_folder <- "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/"
outFile <- paste(output_folder, "/ViolinPlot.goi_top16_Tlx1+markers_E15.5spleen.in_scRNA-seq_E14.5pancreas.multipanel.pdf", sep = "")

# Save the plot as a PDF
pdf(outFile, width = 30, height = 20)
plot(p)
dev.off()

## Heatmap of Top100 Tlx1+ cell markers from E15.5 spleen, with E14.5 pancreas scRNA-seq data
# 1. Extract the gene names of the top 100 markers from top100_tlx1_spleen_markers
top100_gene_names <- top100_tlx1_spleen_markers[[1]]  # Assuming this is the gene name column

# 2. Ensure the genes are present in the Seurat object
genes_in_seurat <- intersect(top100_gene_names, rownames(so_Fgf9WT_processed))
length(genes_in_seurat)  # Check how many of your top100 markers are present in the Seurat object

# 3. Extract expression data for the selected genes
expression_data <- FetchData(so_Fgf9WT_processed, vars = genes_in_seurat)
head(expression_data)   # Check if expression data is correctly extracted

# 4. Set the active assay to RNA if you want to plot raw expression
DefaultAssay(so_Fgf9WT_processed) <- "SCT"

# 5. Create the heatmap using raw expression data
p <- DoHeatmap(so_Fgf9WT_processed,
               features = genes_in_seurat, 
               group.by = "seurat_clusters",
               size = 3, 
               slot = "data")  # Use log-normalized data (LogNormalize in RNA assay)

output_folder <- "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/"
outFile <- paste(output_folder, "/Heatmap.goi_top100_Tlx1+markers_E15.5spleen.in_scRNA-seq_E14.5pancreas.pdf", sep = "")

# Save the plot as a PDF
pdf(outFile, width = 30, height = 20)
plot(p)
dev.off()

# UMAP plot of E14.5 pancreas scRNA-seq data (clusters)
p <- DimPlot(object = so_Fgf9WT_processed,
             reduction = 'umap',
             group.by = 'seurat_clusters',
             label = TRUE)
outFile <- paste(output_folder, "/UMAP.scRNA-seq_E14.5pancreas.clusters.pdf", sep = "")
pdf(outFile, width = 7, height = 5)
plot(p)
dev.off()


################################################################################
## Test for correlated expression of Tlx1 and Barx1 in E14.5 pancreas
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
outFile <- paste(output_folder, "/UMAP.scRNA-seq_E14.5pancreas.mutual_exclusive_Tlx1+Barx1.pdf", sep = "")
pdf(outFile, width = 10, height = 7)
plot(p)
dev.off()

# Visualize mutual exclusivity with a Venn diagram (with two overlapping circles)
library(VennDiagram)
library(grid)

# Define the output file path for the Venn diagram
outFile_venn_WT <- paste(output_folder, "/VennDiagram.scRNA-seq_E14.5pancreas.Tlx1_Barx1_expressed.pdf", sep = "")

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

################################################################################
## Is Nr2f2 (candidate from spleen) differentially regulated between WT and Fgf9 null pancreas?

# Retrieve dataset from previous analysis
so_pancreas_Fgf9mut <- readRDS("/Users/veralaub/Documents/postdoc/collaboration/Maurizio/Fgf9null_datasets/scRNA-seq/pancreas_Fgf9null/results/so_pancreas_Fgf9.rds")

# UMAP that differnetiates between datasets
# Plot UMAP, colored by 'orig.ident' (which stores the dataset origin)
p <- DimPlot(so_pancreas_Fgf9mut, 
             group.by = "orig.ident", 
             reduction = "umap")
outFile_pancreas_Fgf9 <- paste(output_folder, "/UMAP.Fgf9mut.orig.ident.pdf", sep = "")
pdf(outFile_pancreas_Fgf9, width = 7, height = 5)
plot(p)
dev.off()                                          

#GOI: UMAP overlay
goi <- c("Tlx1", "Barx1", "Nr2f2", "Fgf9")
outFile_pancreas_Fgf9 <- paste(output_folder, "/UMAP.Fgf9mut.goi.orig.ident.pdf", sep = "")
pdf(outFile_pancreas_Fgf9, width = 10, height = 5)
for (gene in goi) {
  p <- FeaturePlot(so_pancreas_Fgf9mut, 
                   features = gene,
                   split.by = "orig.ident")
  plot(p)
}
dev.off() 

################################################################################
## Integration of spleen E15.5 with pancreas E14.5 (both WT)

# Load required libraries
library(Seurat)
library(patchwork)
output_folder <- "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/integrated_analysis_spleen+pancreas/"

# Load raw data
d_spleenE15.5 <- Read10X(data.dir ="/Users/veralaub/Documents/postdoc/collaboration/Maurizio/E15.5_spleen/MR3_filtered_feature_bc_matrix/")
so_spleenE15.5 <- CreateSeuratObject(counts = d_spleenE15.5,
                                     project = "spleen_E15.5",
                                     min.cells = 3, min.features = 200)

d_pancreasE14.5 <- Read10X(data.dir ="/Users/veralaub/Documents/postdoc/collaboration/Maurizio/Fgf9null_datasets/scRNA-seq/pancreas_Fgf9null/E14.5_pancreas_Fgf9WT_GSM6434043/")
so_pancreasE14.5 <- CreateSeuratObject(counts = d_pancreasE14.5, 
                                                project = "pancreas_E14.5",
                                                min.cells = 3, min.features = 200)

# Merge E14.5 pancreas and E15.5 spleen datasets (stored in different layers)
so_spleenE15.5_pancreasE14.5 <- merge(so_spleenE15.5, y = so_pancreasE14.5,
                                      add.cell.ids = ls()[1:2],
                                      project = "spleenE15.5_pancreasE14.5_merged")
# Add mito fraction
so_spleenE15.5_pancreasE14.5[["percent.mt"]] <- PercentageFeatureSet(so_spleenE15.5_pancreasE14.5,
                                                                     pattern = "^mt-")

# QC stats before filtering
p <- RidgePlot(so_spleenE15.5_pancreasE14.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "data.qc.unfiltered.pdf", sep = "")
pdf(out_path, width = 7, height = 10)
plot(p)
dev.off()

# Define thresholds for filtering cells (can be adapted ad gusto)
# Only such cells that pass these criteria are kept for further analysis
# The current threshold filters for rather lowly expressed genes
percent.mt_max <- 5  # maximum percentage of mitochondrial genes (adjust as needed)
nFeature_RNA_min <- 1500  # minimum number of features per cell
nFeature_RNA_max <- 10000 # minimum number of features per cell
nCount_RNA_min <- 1000  # minimum number of RNA counts per cell
nCount_RNA_max <- 50000  # maximum number of RNA counts per cell

# Subset the Seurat object (filter based on thresholds above)
so_spleenE15.5_pancreasE14.5 <- subset(so_spleenE15.5_pancreasE14.5,
                                       subset = percent.mt <= percent.mt_max &
                                       nFeature_RNA >= nFeature_RNA_min &
                                       nFeature_RNA <= nFeature_RNA_max &
                                       nCount_RNA >= nCount_RNA_min &
                                       nCount_RNA <= nCount_RNA_max)

# QC stats after filtering
p <- RidgePlot(so_spleenE15.5_pancreasE14.5,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "data.qc.filtered.pdf", sep = "")
pdf(out_path, width = 7, height = 10)
plot(p)
dev.off()

# Parallelize process
library(future)
plan("multicore", workers = 10)  # Set the number of parallel workers

# Increase the maximum global size to 32 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 32 * 1024 * 1024 * 1024)

# Normalization with NormalizeData
so_spleenE15.5_pancreasE14.5 <- NormalizeData(so_spleenE15.5_pancreasE14.5)

# Identification of VariableFeatures
so_spleenE15.5_pancreasE14.5 <- FindVariableFeatures(so_spleenE15.5_pancreasE14.5)

# Scale Seurat object
so_spleenE15.5_pancreasE14.5 <- ScaleData(so_spleenE15.5_pancreasE14.5)

# Perform PCA analysis
so_spleenE15.5_pancreasE14.5 <- RunPCA(so_spleenE15.5_pancreasE14.5,
                                       verbose = FALSE, npcs = 20)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_spleenE15.5_pancreasE14.5)
out_path <- paste(output_folder, "data.qc.ellbowplot.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

# Find Neighbors and perform Clustering (Louvian)
so_spleenE15.5_pancreasE14.5 <- FindNeighbors(so_spleenE15.5_pancreasE14.5,
                                              dims = 1:20, 
                                              reduction = "pca")
so_spleenE15.5_pancreasE14.5 <- FindClusters(so_spleenE15.5_pancreasE14.5,
                                             resolution = 0.8, 
                                            # algorithm = 4,    # This would use Leiden algorithm, not feasible due to error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()"
                                             cluster.name = "unintegrated_clusters",
                                             verbose = TRUE)

# Run UMAP and save as Dimplot
so_spleenE15.5_pancreasE14.5 <- RunUMAP(so_spleenE15.5_pancreasE14.5, 
                                        dims = 1:10, 
                                        reduction = "pca", 
                                        reduction.name = "umap.unintegrated")


p <- DimPlot(so_spleenE15.5_pancreasE14.5, 
             reduction = "umap.unintegrated", 
             group.by = c("orig.ident", "seurat_clusters"))
out_path <- paste(output_folder, "data.qc.UMAP.unintergrated.pdf", sep = "")
pdf(out_path, width = 15, height = 5)
plot(p)
dev.off()

# Visualize goi as FeaturePlot (merged but unintegrated datasets)
goi <- c("Tlx1", "Barx1", "Nr2f2", "Fgf9")
out_path <- paste(output_folder, 
                  "/UMAP.spleenE15.5_pancreasE14.5_unintegrated.goi.orig.ident.pdf", sep = "")
pdf(out_path, width = 10, height = 5)
for (gene in goi) {
  p <- FeaturePlot(so_spleenE15.5_pancreasE14.5, 
                   features = gene,
                   split.by = "orig.ident")
  plot(p)
}
dev.off()  


######################## Data integration approach 1 ###########################
## Filtering of datasets individually, integration afterwards
# Add mito fraction
so_spleenE15.5[["percent.mt"]] <- PercentageFeatureSet(so_spleenE15.5,
                                                       pattern = "^mt-")
so_pancreasE14.5[["percent.mt"]] <- PercentageFeatureSet(so_pancreasE14.5,
                                                         pattern = "^mt-")

# QC stats before filtering
p <- RidgePlot(so_spleenE15.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "int_appr1.data.spleenE15.5.qc.unfiltered.pdf", sep = "")
pdf(out_path, width = 7, height = 10)
plot(p)
dev.off()

p <- RidgePlot(so_pancreasE14.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "int_appr1.data.pancreasE14.5.qc.unfiltered.pdf", sep = "")
pdf(out_path, width = 7, height = 10)
plot(p)
dev.off()

# Filter datasets
percent.mt_max <- 3  # maximum percentage of mitochondrial genes (adjust as needed)
nFeature_RNA_min <- 1000  # minimum number of features per cell
nFeature_RNA_max <- 10000 # minimum number of features per cell
nCount_RNA_min <- 1000  # minimum number of RNA counts per cell
nCount_RNA_max <- 50000  # maximum number of RNA counts per cell
so_spleenE15.5 <- subset(so_spleenE15.5,
                         subset = percent.mt <= percent.mt_max &
                         nFeature_RNA >= nFeature_RNA_min &
                         nFeature_RNA <= nFeature_RNA_max &
                         nCount_RNA >= nCount_RNA_min &
                         nCount_RNA <= nCount_RNA_max)

percent.mt_max <- 3  # maximum percentage of mitochondrial genes (adjust as needed)
nFeature_RNA_min <- 2000  # minimum number of features per cell
nFeature_RNA_max <- 10000 # minimum number of features per cell
nCount_RNA_min <- 1000  # minimum number of RNA counts per cell
nCount_RNA_max <- 50000  # maximum number of RNA counts per cell
so_pancreasE14.5 <- subset(so_pancreasE14.5,
                           subset = percent.mt <= percent.mt_max &
                           nFeature_RNA >= nFeature_RNA_min &
                           nFeature_RNA <= nFeature_RNA_max &
                           nCount_RNA >= nCount_RNA_min &
                           nCount_RNA <= nCount_RNA_max)

# QC stats after filtering
p <- RidgePlot(so_spleenE15.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "int_appr1.data.spleenE15.5.qc.filtered.pdf", sep = "")
pdf(out_path, width = 7, height = 10)
plot(p)
dev.off()

p <- RidgePlot(so_pancreasE14.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE)
out_path <- paste(output_folder, "int_appr1.data.pancreasE14.5.qc.filtered.pdf", sep = "")
pdf(out_path, width = 7, height = 10)
plot(p)
dev.off()

# Normalization of individual datasets with SCTransform (for variable features)
n_features <- 50
so_spleenE15.5 <- SCTransform(so_spleenE15.5,
                              verbose = TRUE,
                              variable.features.n = n_features)
so_pancreasE14.5 <- SCTransform(so_pancreasE14.5,
                                verbose = TRUE,
                                variable.features.n = n_features)

# Scale Seurat object
so_spleenE15.5 <- ScaleData(so_spleenE15.5)
so_pancreasE14.5 <- ScaleData(so_pancreasE14.5)

# Perform PCA analysis
so_spleenE15.5 <- RunPCA(so_spleenE15.5,
                         verbose = FALSE, npcs = 20)
so_pancreasE14.5 <- RunPCA(so_pancreasE14.5,
                           verbose = FALSE, npcs = 20)

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_spleenE15.5)
out_path <- paste(output_folder, "int_appr1.data.spleenE15.5.qc.ellbowplot.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

p <- ElbowPlot(so_pancreasE14.5)
out_path <- paste(output_folder, "int_appr1.data.pancreasE14.5.qc.ellbowplot.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

# Run FindVariableFeatures after SCTransform
so_spleenE15.5 <- FindVariableFeatures(so_spleenE15.5, 
                                       selection.method = "vst", 
                                       nfeatures = 3000)
so_pancreasE14.5 <- FindVariableFeatures(so_pancreasE14.5, 
                                         selection.method = "vst", 
                                         nfeatures = 3000)

# Change the default assay to "SCT"
# Set the default assay to "SCT" for both objects
DefaultAssay(so_spleenE15.5) <- "SCT"
DefaultAssay(so_pancreasE14.5) <- "SCT"

# Subset the Seurat objects to a smaller number of cells for testing purposes
# Used randomized subset of cells, since full datasets did not work (due to vector memory limit of 18.0 Gb)
so_spleenE15.5_sub <- so_spleenE15.5[, sample(1:ncol(so_spleenE15.5), 2000)]
so_pancreasE14.5_sub <- so_pancreasE14.5[, sample(1:ncol(so_pancreasE14.5), 2000)]

# Run the integration on the subset
anchors <- FindIntegrationAnchors(object.list = list(so_spleenE15.5_sub, so_pancreasE14.5_sub), 
                                  dims = 1:15, 
                                  anchor.features = 1000,
                                  verbose = TRUE)
so_spleenE15.5_pancreasE14.5_integrated <- IntegrateData(anchorset = anchors, 
                                                         dims = 1:15,
                                                         verbose = TRUE)

# Scale the integrated data
so_spleenE15.5_pancreasE14.5_integrated <- ScaleData(so_spleenE15.5_pancreasE14.5_integrated)

# Run PCA
so_spleenE15.5_pancreasE14.5_integrated <- RunPCA(so_spleenE15.5_pancreasE14.5_integrated, 
                                                  verbose = FALSE)

# Run UMAP
so_spleenE15.5_pancreasE14.5_integrated <- RunUMAP(so_spleenE15.5_pancreasE14.5_integrated, 
                                                   dims = 1:20)

# Find neighbors and perform clustering
so_spleenE15.5_pancreasE14.5_integrated <- FindNeighbors(so_spleenE15.5_pancreasE14.5_integrated, 
                                                         dims = 1:20)
so_spleenE15.5_pancreasE14.5_integrated <- FindClusters(so_spleenE15.5_pancreasE14.5_integrated, 
                                                        resolution = 0.8)

# Visualize the UMAP
p <- DimPlot(so_spleenE15.5_pancreasE14.5_integrated, 
             reduction = "umap",
             group.by = 'seurat_clusters',
             label = TRUE)
outFile <- paste(output_folder, "/int_appr1.UMAP.spleenE15.5_pancreasE14.5_integrated.clusters.pdf", sep = "")
pdf(outFile, width = 7, height = 5)
plot(p)
dev.off()

p <- DimPlot(so_spleenE15.5_pancreasE14.5_integrated, 
             reduction = "umap",
             group.by = "orig.ident",
             label = TRUE)
outFile <- paste(output_folder, "/int_appr1.UMAP.spleenE15.5_pancreasE14.5_integrated.orig.ident.pdf", sep = "")
pdf(outFile, width = 7, height = 5)
plot(p)
dev.off()

# Visualize as FeaturePlot
goi <- c("Tlx1", "Barx1", "Nr2f2", "Fgf9")
outFile <- paste(output_folder,
                 "/int_appr1.UMAP.spleenE15.5_pancreasE14.5_integrated.goi.orig.ident.pdf", sep = "")
pdf(outFile, width = 10, height = 5)
for (gene in goi) {
  p <- FeaturePlot(so_spleenE15.5_pancreasE14.5_integrated, 
                   features = gene,
                   split.by = "orig.ident")
  plot(p)
}
dev.off()           

# Save the Seurat object
outFile_spleenE15.5_pancreasE14.5_integrated <- paste(output_folder, "int_appr1.so_spleenE15.5_pancreasE14.5_integrated.rds", sep = "")
saveRDS(so_spleenE15.5_pancreasE14.5_integrated, file = outFile_spleenE15.5_pancreasE14.5_integrated)


######################## Data integration approach 2 ###########################
## Following Vignette on https://satijalab.org/seurat/articles/integration_introduction
# Using merged and pre-processed Seurat object with spleen E15.5 and pancreas E14.5 data

# Perform integration
so_spleenE15.5_pancreasE14.5 <- IntegrateLayers(object = so_spleenE15.5_pancreasE14.5, 
                                                method = CCAIntegration, 
                                                orig.reduction = "pca", 
                                                new.reduction = "integrated.cca",
                                                verbose = TRUE)

# Re-join layers after integration
so_spleenE15.5_pancreasE14.5[["RNA"]] <- JoinLayers(so_spleenE15.5_pancreasE14.5[["RNA"]])

# Elbow plot to determine the number of PCs
p <- ElbowPlot(so_spleenE15.5_pancreasE14.5)
out_path <- paste(output_folder, "int_appr2.data.spleenE15.5_pancreasE14.5_integrated.qc.ellbowplot.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()

# Find Neighbors and Clusters, based on PCs from Ellbowplot
so_spleenE15.5_pancreasE14.5 <- FindNeighbors(so_spleenE15.5_pancreasE14.5, 
                                              reduction = "integrated.cca", 
                                              dims = 1:20)
so_spleenE15.5_pancreasE14.5 <- FindClusters(so_spleenE15.5_pancreasE14.5, 
                                             resolution = 0.8)

so_spleenE15.5_pancreasE14.5 <- RunUMAP(so_spleenE15.5_pancreasE14.5, 
                                        dims = 1:20, 
                                        reduction = "integrated.cca")

# Visualization
p <- DimPlot(so_spleenE15.5_pancreasE14.5, 
             reduction = "umap", 
             group.by = c("orig.ident", "seurat_clusters"))
outFile <- paste(output_folder, "/int_appr2.UMAP.spleenE15.5_pancreasE14.5_integrated.orig.ident.pdf", sep = "")
pdf(outFile, width = 12, height = 5)
plot(p)
dev.off()


# Visualize as FeaturePlot
FeaturePlot(so_spleenE15.5_pancreasE14.5, 
            features = "Tlx1",
            reduction = "umap",
            split.by = "orig.ident")

goi <- c("Tlx1", "Barx1", "Nr2f2", "Fgf9", "Fgf10",
         "Mki67", "Cdk1", "Cdkn1c",
         "Pi15", "Neurl3", "Ackr3", "Tnni1", "Dpt", "Itih2", "Klhl41", "Cybrd1", # upregulated in bulk RNA-seq with Fgf9 mutation (pancreas E14.5)
         "Calcrl", "F2", "Fibin", "Myl9", "Atp1b4", "Col4a6", "Car2", "Sfrp2", 
         "Synpo2", "Col25a1", "Tpm2", "Mmp23", "Sema3e", "Cpz", "Chrm2", "Rarres2", 
         "Actg2", "Cxcl12", "Vwf", "Mgp", "Rerg", "Gm18194", "Ckm", "Clec11a", 
         "Atp2a1", "Tnnt3", "Vgll2", "Oit3", "Dcn", "Lum", "Tspan8", "Lyz2", 
         "Gli1", "Ndufa4l2", "Ank1", "Hmox1", "Ndrg4", "Tnnc1", "Gdf10", "Trac", 
         "Myh7", "Dmtn", "Ednrb", "Cnn1", "Bmper", "Apoa1", "Isl2", "Loxl1", 
         "Prss35", "Efemp1", "Slc4a1", "Abca8a", "Serpinb6b", "Edn1", "Gm48398", 
         "Colec11", "Ahnak2", "Nov", "Sox10", "Igfbp6", "Masp1", "Robo2", "Chodl", 
         "Tmem181b-ps", "Msln", "C3", "Dsc3", "Sh3rf2", "Onecut2", "Mpeg1", 
         "Anxa1", "Prkg1", "Acta2")
outFile <- paste(output_folder,
                 "/int_appr2.UMAP.spleenE15.5_pancreasE14.5_integrated.goi.orig.ident.pdf", 
                 sep = "")
pdf(outFile, width = 12, height = 5)
# Loop through each gene and check if it exists in the Seurat object
for (gene in goi) {
  if (gene %in% rownames(so_spleenE15.5_pancreasE14.5)) {
    # Plot only if the gene is found in the Seurat object
    p <- FeaturePlot(so_spleenE15.5_pancreasE14.5, 
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

# Save the Seurat object
outFile_spleenE15.5_pancreasE14.5_integrated <- paste(output_folder, "int_appr2.so_spleenE15.5_pancreasE14.5_integrated.rds", sep = "")
saveRDS(so_spleenE15.5_pancreasE14.5, file = outFile_spleenE15.5_pancreasE14.5_integrated)


############## Test for mutual exclusivity of candidate expression #############

## Is Tlx1 co-expressed or inversely correlated with Sox10 (TF that is overexpressed in Fgf9 null pancreas)?
# 1. Define gene expression thresholds
# Extract expression data for Tlx1 and Sox10
Tlx1_expr <- FetchData(so_spleenE15.5_pancreasE14.5, vars = "Tlx1")
Sox10_expr <- FetchData(so_spleenE15.5_pancreasE14.5, vars = "Sox10")

# Define a threshold for expression (e.g., greater than 0 means expressed)
Tlx1_expressed <- Tlx1_expr > 0
Sox10_expressed <- Sox10_expr > 0

# 2. Identify cells where only one gene is expressed
# Identify cells where Tlx1 is expressed but Sox10 is not, and vice versa
mutually_exclusive_Tlx1 <- Tlx1_expressed & !Sox10_expressed
mutually_exclusive_Sox10 <- Sox10_expressed & !Tlx1_expressed
both_expressed <- Tlx1_expressed & Sox10_expressed

# Count the number of mutually exclusive cells for each gene
mutually_exclusive_Tlx1_cells <- sum(mutually_exclusive_Tlx1)
mutually_exclusive_Sox10_cells <- sum(mutually_exclusive_Sox10)

# Count the number of cells where both genes are expressed
both_expressed_cells <- sum(both_expressed)

# Print the results
cat("Number of cells where Tlx1 is expressed but Sox10 is not: ", mutually_exclusive_Tlx1_cells, "\n")
cat("Number of cells where Sox10 is expressed but Tlx1 is not: ", mutually_exclusive_Sox10_cells, "\n")
cat("Number of cells where both Tlx1 and Sox10 are expressed: ", both_expressed_cells, "\n")

# 3. Statistical test for mutual exclusivity
# Create a contingency table for the co-expression of Tlx1 and Sox10
contingency_table <- table(Tlx1_expressed, Sox10_expressed)

# Perform Fisher's Exact Test (for small numbers) or Chi-square test
# Fisher's exact test is recommended for small sample sizes (less than 5 in any cell)
fisher_test_result <- fisher.test(contingency_table)

# Print the p-value from the Fisher's test
cat("P-value for mutual exclusivity (Fisher's Exact Test): ", fisher_test_result$p.value, "\n")

# 4. Visualizing the results
# Create a new metadata column to label cells as mutually exclusive for each gene
so_spleenE15.5_pancreasE14.5$MutualExclusive_Tlx1 <- mutually_exclusive_Tlx1
so_spleenE15.5_pancreasE14.5$MutualExclusive_Sox10 <- mutually_exclusive_Sox10
so_spleenE15.5_pancreasE14.5$Both_Expressed <- both_expressed

# Visualize mutual exclusivity of Tlx1 and Sox10 using UMAP
library(patchwork)

# Generate individual FeaturePlots for each feature
p1 <- FeaturePlot(so_spleenE15.5_pancreasE14.5, 
                  features = "MutualExclusive_Tlx1", 
                  reduction = "umap")
p2 <- FeaturePlot(so_spleenE15.5_pancreasE14.5, 
                  features = "MutualExclusive_Sox10", 
                  reduction = "umap")
p3 <- FeaturePlot(so_spleenE15.5_pancreasE14.5, 
                  features = "Both_Expressed", 
                  reduction = "umap")

# Combine the plots into one row
combined_plot <- p1 + p2 + p3 + plot_layout(ncol = 3)  # Arrange them in 1 row, 3 columns

# Define the output file path
outFile <- paste(output_folder, "/int_appr2.UMAP.Mutual_exclusive_Tlx1+Sox10.multipanel.pdf", sep = "")

# Save the combined plot as a PDF
pdf(outFile, width = 15, height = 5)  # Adjust the height/width for your preference
print(combined_plot)  # Print the combined plot to the PDF
dev.off()

# Visualize mutual exclusivity with a Venn diagram (with two overlapping circles)
library(VennDiagram)
library(grid)
library(futile.logger)

# Define the output file path for the Venn diagram
outFile_venn <- paste(output_folder, "/int_appr2.VennDiagram.MutualExclusiveExpression_Tlx1+Sox10.pdf", sep = "")

# Open a PDF device to save the plot
pdf(outFile_venn, width = 10, height = 5)

# Create the Venn diagram with two circles
venn.plot <- venn.diagram(
  x = list(
    "Tlx1 expressed" = which(Tlx1_expressed),
    "Sox10 expressed" = which(Sox10_expressed)
  ),
  category.names = c("Tlx1 expressed", "Sox10 expressed"),
  filename = NULL,  # We are using grid.draw() to plot, so no need for a file name here
  output = TRUE,
  lwd = 2,  # Line width of the circles
  fill = c("red", "blue"),  # Fill color for the circles
  alpha = c(0.5, 0.5),  # Transparency of the circles
  cex = 1.5,  # Text size for main title
  cat.cex = 1.5,  # Category text size
  cat.pos = 0,  # Category text position (0 is top-center)
  main = "Venn Diagram of Gene Expression, spleenE15.5 + pancreasE14.5",  # Main title
  family = "sans",  # Use a generic sans-serif font (usually Helvetica or Arial)
  cat.fontface = 1,  # Regular font style for category names
  fontface = 1  # Regular font style for main title
)

# Plot the Venn diagram
grid.draw(venn.plot)

# Close the PDF device to save the plot
dev.off()

### Stratification of mutual expression by orig.ident (pancreas and spleen data)
# Generate individual FeaturePlots for each feature, stratified by 'orig.ident'
p1 <- FeaturePlot(so_spleenE15.5_pancreasE14.5, 
                  features = "MutualExclusive_Tlx1", 
                  reduction = "umap",
                  split.by = "orig.ident")
p2 <- FeaturePlot(so_spleenE15.5_pancreasE14.5, 
                  features = "MutualExclusive_Sox10", 
                  reduction = "umap",
                  split.by = "orig.ident")
p3 <- FeaturePlot(so_spleenE15.5_pancreasE14.5, 
                  features = "Both_Expressed", 
                  reduction = "umap",
                  split.by = "orig.ident")

# Define the output file paths for each plot
outFile1 <- paste(output_folder, "/int_appr2.UMAP.MutualExclusive_Tlx1(noSox10).orig.ident.pdf", sep = "")
outFile2 <- paste(output_folder, "/int_appr2.UMAP.MutualExclusive_Sox10(noTlx1).orig.ident.pdf", sep = "")
outFile3 <- paste(output_folder, "/int_appr2.UMAP.Both_Expressed(Tlx1+Sox10).orig.ident.pdf", sep = "")

# Save each plot as a separate PDF
pdf(outFile1, width = 10, height = 5)
print(p1)  # Print the plot p1 to the PDF
dev.off()

pdf(outFile2, width = 10, height = 5)
print(p2)  # Print the plot p2 to the PDF
dev.off()

pdf(outFile3, width = 10, height = 5)
print(p3)  # Print the plot p3 to the PDF
dev.off()

# Generate individual ViolinPlots for each feature, stratified by 'seurat_clusters'
p1 <- VlnPlot(so_spleenE15.5_pancreasE14.5, 
              features = "MutualExclusive_Tlx1", 
              group.by = "seurat_clusters")  # Stratified by 'seurat_clusters'

p2 <- VlnPlot(so_spleenE15.5_pancreasE14.5, 
              features = "MutualExclusive_Sox10", 
              group.by = "seurat_clusters")  # Stratified by 'seurat_clusters'

p3 <- VlnPlot(so_spleenE15.5_pancreasE14.5, 
              features = "Both_Expressed", 
              group.by = "seurat_clusters")  # Stratified by 'seurat_clusters'

# Define the output file paths for each plot
outFile1 <- paste(output_folder, "/int_appr2.Violin.MutualExclusive_Tlx1(noSox10).seurat_clusters.pdf", sep = "")
outFile2 <- paste(output_folder, "/int_appr2.Violin.MutualExclusive_Sox10(noTlx1).seurat_clusters.pdf", sep = "")
outFile3 <- paste(output_folder, "/int_appr2.Violin.Both_Expressed(Tlx1+Sox10).seurat_clusters.pdf", sep = "")

# Save each plot as a separate PDF
pdf(outFile1, width = 10, height = 7)
print(p1)  # Print the plot p1 to the PDF
dev.off()

pdf(outFile2, width = 10, height = 7)
print(p2)  # Print the plot p2 to the PDF
dev.off()

pdf(outFile3, width = 10, height = 7)
print(p3)  # Print the plot p3 to the PDF
dev.off()


## Is Tlx1 co-expressed or inversely correlated with Klf4 (TF candidate from spleen project)?
# 1. Define gene expression thresholds
# Extract expression data for Tlx1 and Klf4
Tlx1_expr <- FetchData(so_spleenE15.5_pancreasE14.5, vars = "Tlx1")
Klf4_expr <- FetchData(so_spleenE15.5_pancreasE14.5, vars = "Klf4")

# Define a threshold for expression (e.g., greater than 0 means expressed)
Tlx1_expressed <- Tlx1_expr > 0
Klf4_expressed <- Klf4_expr > 0

# 2. Identify cells where only one gene is expressed
# Identify cells where Tlx1 is expressed but Klf4 is not, and vice versa
mutually_exclusive_Tlx1 <- Tlx1_expressed & !Klf4_expressed
mutually_exclusive_Klf4 <- Klf4_expressed & !Tlx1_expressed
both_expressed <- Tlx1_expressed & Klf4_expressed

# Count the number of mutually exclusive cells for each gene
mutually_exclusive_Tlx1_cells <- sum(mutually_exclusive_Tlx1)
mutually_exclusive_Klf4_cells <- sum(mutually_exclusive_Klf4)

# Count the number of cells where both genes are expressed
both_expressed_cells <- sum(both_expressed)

# Print the results
cat("Number of cells where Tlx1 is expressed but Klf4 is not: ", mutually_exclusive_Tlx1_cells, "\n")
cat("Number of cells where Klf4 is expressed but Tlx1 is not: ", mutually_exclusive_Klf4_cells, "\n")
cat("Number of cells where both Tlx1 and Klf4 are expressed: ", both_expressed_cells, "\n")

# 3. Statistical test for mutual exclusivity
# Create a contingency table for the co-expression of Tlx1 and Klf4
contingency_table <- table(Tlx1_expressed, Klf4_expressed)

# Perform Fisher's Exact Test (for small numbers) or Chi-square test
# Fisher's exact test is recommended for small sample sizes (less than 5 in any cell)
fisher_test_result <- fisher.test(contingency_table)

# Print the p-value from the Fisher's test
cat("P-value for mutual exclusivity (Fisher's Exact Test): ", fisher_test_result$p.value, "\n")

# 4. Visualizing the results
# Create a new metadata column to label cells as mutually exclusive for each gene
so_spleenE15.5_pancreasE14.5$MutualExclusive_Tlx1 <- mutually_exclusive_Tlx1
so_spleenE15.5_pancreasE14.5$MutualExclusive_Klf4 <- mutually_exclusive_Klf4
so_spleenE15.5_pancreasE14.5$Both_Expressed <- both_expressed

# Visualize mutual exclusivity of Tlx1 and Klf4 using UMAP
library(patchwork)

# Generate individual FeaturePlots for each feature
p1 <- FeaturePlot(so_spleenE15.5_pancreasE14.5, 
                  features = "MutualExclusive_Tlx1", 
                  reduction = "umap")
p2 <- FeaturePlot(so_spleenE15.5_pancreasE14.5, 
                  features = "MutualExclusive_Klf4", 
                  reduction = "umap")
p3 <- FeaturePlot(so_spleenE15.5_pancreasE14.5, 
                  features = "Both_Expressed", 
                  reduction = "umap")

# Combine the plots into one row
combined_plot <- p1 + p2 + p3 + plot_layout(ncol = 3)  # Arrange them in 1 row, 3 columns

# Define the output file path
outFile <- paste(output_folder, "/int_appr2.UMAP.Mutual_exclusive_Tlx1+Klf4.multipanel.pdf", sep = "")

# Save the combined plot as a PDF
pdf(outFile, width = 15, height = 5)  # Adjust the height/width for your preference
print(combined_plot)  # Print the combined plot to the PDF
dev.off()

# Visualize mutual exclusivity with a Venn diagram (with two overlapping circles)
library(VennDiagram)
library(grid)
library(futile.logger)

# Define the output file path for the Venn diagram
outFile_venn <- paste(output_folder, "/int_appr2.VennDiagram.MutualExclusiveExpression_Tlx1+Klf4.pdf", sep = "")

# Open a PDF device to save the plot
pdf(outFile_venn, width = 10, height = 5)

# Create the Venn diagram with two circles
venn.plot <- venn.diagram(
  x = list(
    "Tlx1 expressed" = which(Tlx1_expressed),
    "Klf4 expressed" = which(Klf4_expressed)
  ),
  category.names = c("Tlx1 expressed", "Klf4 expressed"),
  filename = NULL,  # We are using grid.draw() to plot, so no need for a file name here
  output = TRUE,
  lwd = 2,  # Line width of the circles
  fill = c("red", "blue"),  # Fill color for the circles
  alpha = c(0.5, 0.5),  # Transparency of the circles
  cex = 1.5,  # Text size for main title
  cat.cex = 1.5,  # Category text size
  cat.pos = 0,  # Category text position (0 is top-center)
  main = "Venn Diagram of Gene Expression, spleenE15.5 + pancreasE14.5",  # Main title
  family = "sans",  # Use a generic sans-serif font (usually Helvetica or Arial)
  cat.fontface = 1,  # Regular font style for category names
  fontface = 1  # Regular font style for main title
)

# Plot the Venn diagram
grid.draw(venn.plot)

# Close the PDF device to save the plot
dev.off()

### Stratification of mutual expression by orig.ident (pancreas and spleen data)
# Generate individual FeaturePlots for each feature, stratified by 'orig.ident'
p1 <- FeaturePlot(so_spleenE15.5_pancreasE14.5, 
                  features = "MutualExclusive_Tlx1", 
                  reduction = "umap",
                  split.by = "orig.ident")
p2 <- FeaturePlot(so_spleenE15.5_pancreasE14.5, 
                  features = "MutualExclusive_Klf4", 
                  reduction = "umap",
                  split.by = "orig.ident")
p3 <- FeaturePlot(so_spleenE15.5_pancreasE14.5, 
                  features = "Both_Expressed", 
                  reduction = "umap",
                  split.by = "orig.ident")

# Define the output file paths for each plot
outFile1 <- paste(output_folder, "/int_appr2.UMAP.MutualExclusive_Tlx1(noKlf4)_UMAP.orig.ident.pdf", sep = "")
outFile2 <- paste(output_folder, "/int_appr2.UMAP.MutualExclusive_Klf4(noTlx1)_UMAP.orig.ident.pdf", sep = "")
outFile3 <- paste(output_folder, "/int_appr2.UMAP.Both_Expressed(Tlx1+Klf4).orig.ident.pdf", sep = "")

# Save each plot as a separate PDF
pdf(outFile1, width = 10, height = 5)
print(p1)  # Print the plot p1 to the PDF
dev.off()

pdf(outFile2, width = 10, height = 5)
print(p2)  # Print the plot p2 to the PDF
dev.off()

pdf(outFile3, width = 10, height = 5)
print(p3)  # Print the plot p3 to the PDF
dev.off()

# Generate individual ViolinPlots for each feature, stratified by 'seurat_clusters'
p1 <- VlnPlot(so_spleenE15.5_pancreasE14.5, 
              features = "MutualExclusive_Tlx1", 
              group.by = "seurat_clusters")  # Stratified by 'seurat_clusters'

p2 <- VlnPlot(so_spleenE15.5_pancreasE14.5, 
              features = "MutualExclusive_Klf4", 
              group.by = "seurat_clusters")  # Stratified by 'seurat_clusters'

p3 <- VlnPlot(so_spleenE15.5_pancreasE14.5, 
              features = "Both_Expressed", 
              group.by = "seurat_clusters")  # Stratified by 'seurat_clusters'

# Define the output file paths for each plot
outFile1 <- paste(output_folder, "/int_appr2.Violin.MutualExclusive_Tlx1(noKlf4).seurat_clusters.pdf", sep = "")
outFile2 <- paste(output_folder, "/int_appr2.Violin.MutualExclusive_Klf4(noTlx1).seurat_clusters.pdf", sep = "")
outFile3 <- paste(output_folder, "/int_appr2.Violin.Both_Expressed(Tlx1+Klf4).seurat_clusters.pdf", sep = "")

# Save each plot as a separate PDF
pdf(outFile1, width = 10, height = 7)
print(p1)  # Print the plot p1 to the PDF
dev.off()

pdf(outFile2, width = 10, height = 7)
print(p2)  # Print the plot p2 to the PDF
dev.off()

pdf(outFile3, width = 10, height = 7)
print(p3)  # Print the plot p3 to the PDF
dev.off()


######################## Data integration approach 3 ###########################
# Using Nuos E15.5 spleen Seurat object and pre-processed Seurat object of pancreas E14.5 data

### PREPARATION OF DATASETS
# Load required libraries
library(Seurat)
library(patchwork)

# Increase the maximum global size to 32 GB (2 * 1024^3 bytes)
options(future.globals.maxSize = 32 * 1024 * 1024 * 1024)

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/integrated_analysis_spleen+pancreas/"

# Load raw data
so_spleenE15.5 <- readRDS("/Users/veralaub/Documents/postdoc/collaboration/Maurizio/E15.5_spleen/E15.5_MR4_subsetNoWeird.rds")

so_pancreasE14.5 <- readRDS("/Users/veralaub/Documents/postdoc/collaboration/Maurizio/Fgf9null_datasets/scRNA-seq/pancreas_Fgf9null/results_WT/so_pancreas_WT.rds")

# Check memory usage before and after merge (due to Error "vector memory limit of 18.0 Gb reached, see mem.maxVSize()")
pryr::mem_used()

# Merge E14.5 pancreas and E15.5 spleen datasets (stored in different layers)
# Rename the datasets directly in the add.cell.ids argument
so_spleenE15.5_pancreasE14.5 <- merge(so_spleenE15.5, y = so_pancreasE14.5,
                                      add.cell.ids = c("spleenE15.5", "pancreasE14.5"),
                                      project = "spleenE15.5_pancreasE14.5_merged")

# Change the 'orig.ident' metadata to reflect the new names
# Adjust the 'gsub' pattern based on actual 'orig.ident' values
so_spleenE15.5_pancreasE14.5$orig.ident <- gsub("^MR4", "spleen_E15.5", so_spleenE15.5_pancreasE14.5$orig.ident)
so_spleenE15.5_pancreasE14.5$orig.ident <- gsub("^pancreas_E14.5", "pancreas_E14.5", so_spleenE15.5_pancreasE14.5$orig.ident)

# QC stats before filtering
p <- RidgePlot(so_spleenE15.5_pancreasE14.5, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 1, log = TRUE,
               group.by = "orig.ident")
out_path <- paste(output_folder, "int_appr3.data.spleenE15.5_pancreasE14.5.qcRidgepPlot.pdf", sep = "")
pdf(out_path, width = 7, height = 10)
plot(p)
dev.off()

# Normalization with SCTransform
n_features <- 100
so_spleenE15.5_pancreasE14.5 <- SCTransform(so_spleenE15.5_pancreasE14.5,
                                            verbose = TRUE,
                                            variable.features.n = n_features)

# Scale Seurat object
so_spleenE15.5_pancreasE14.5 <- ScaleData(so_spleenE15.5_pancreasE14.5)

# Perform PCA analysis
so_spleenE15.5_pancreasE14.5 <- RunPCA(so_spleenE15.5_pancreasE14.5,
                                       verbose = FALSE, npcs = 30)

# Run UMAP
so_spleenE15.5_pancreasE14.5 <- RunUMAP(so_spleenE15.5_pancreasE14.5, 
                                        dims = 1:30, 
                                        reduction = "pca")

# Visualize datasets as UMAP before Integration
p <- DimPlot(so_spleenE15.5_pancreasE14.5, 
             reduction = "umap", 
             group.by = c("orig.ident", "seurat_clusters"))
outFile <- paste(output_folder, "/int_appr3.UMAP.spleenE15.5_pancreasE14.5_NotIntegrated.pdf", sep = "")
pdf(outFile, width = 12, height = 5)
plot(p)
dev.off()

# Number of features selection by elbow method (you can use elbow plot to decide on the number of PCs)
p <- ElbowPlot(so_spleenE15.5_pancreasE14.5,
               ndims = 30)
out_path <- paste(output_folder, "int_appr3.data.qc.ellbowplot.pdf", sep = "")
pdf(out_path, width = 5, height = 5)
plot(p)
dev.off()


### ACTUAL INTEGRATION PART
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
outFile <- paste(output_folder, "so_spleenE15.5_int_appr3.rds", sep = "")
saveRDS(so_spleenE15.5, file = outFile)
outFile <- paste(output_folder, "so_pancreasE14.5_int_appr3.rds", sep = "")
saveRDS(so_pancreasE14.5, file = outFile)

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

# Remove the RNA assay from the anchors object to save memory
anchors@object.list[[1]]@assays[["RNA"]] <- NULL
anchors@object.list[[2]]@assays[["RNA"]] <- NULL

# Step 5: Integrate the datasets using the found anchors
so_spleenE15.5_pancreasE14.5_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

# Step 6: Perform scaling and PCA on the integrated data
so_spleenE15.5_pancreasE14.5_integrated <- ScaleData(so_spleenE15.5_pancreasE14.5_integrated)
so_spleenE15.5_pancreasE14.5_integrated <- RunPCA(so_spleenE15.5_pancreasE14.5_integrated, verbose = FALSE)

# Perform UMAP on integrated data
so_spleenE15.5_pancreasE14.5_integrated <- RunUMAP(so_spleenE15.5_pancreasE14.5_integrated, 
                                        dims = 1:13, 
                                        reduction = "pca", 
                                        reduction.name = "umap.integrated")

# Change the 'orig.ident' metadata of E15.5 spleen (to match name)
so_spleenE15.5_pancreasE14.5_integrated$orig.ident <- gsub("^MR4", "spleen_E15.5", so_spleenE15.5_pancreasE14.5_integrated$orig.ident)

# Visualize datasets as UMAP after Integration
p <- DimPlot(so_spleenE15.5_pancreasE14.5_integrated, 
             reduction = "umap.integrated", 
             group.by = c("orig.ident", "seurat_clusters"))
outFile <- paste(output_folder, "/int_appr3.UMAP.spleenE15.5_pancreasE14.5_integrated.pdf", sep = "")
pdf(outFile, width = 12, height = 5)
plot(p)
dev.off()

# Save the Seurat object
outFile <- paste(output_folder, "int_appr3.so_spleenE15.5_pancreasE14.5_integrated.rds", sep = "")
saveRDS(so_spleenE15.5_pancreasE14.5_integrated, file = outFile)

### Visualization
# Change the default assay to "SCT"
DefaultAssay(so_spleenE15.5_pancreasE14.5_integrated) <- "SCT"

# Visualize as FeaturePlot
FeaturePlot(so_spleenE15.5_pancreasE14.5_integrated, 
            features = "Tlx1",
            reduction = "umap.integrated",
            split.by = "orig.ident")

goi <- c("Tlx1", "Barx1", "Nr2f2", "Klf4", "Fgf9", "Fgf10",
         "Mki67", "Cdk1", "Cdkn1c", 
         "Pi15", "Neurl3", "Ackr3", "Tnni1", "Dpt", "Itih2", "Klhl41", "Cybrd1", # upregulated in bulk RNA-seq with Fgf9 mutation (pancreas E14.5)
         "Calcrl", "F2", "Fibin", "Myl9", "Atp1b4", "Col4a6", "Car2", "Sfrp2", 
         "Synpo2", "Col25a1", "Tpm2", "Mmp23", "Sema3e", "Cpz", "Chrm2", "Rarres2", 
         "Actg2", "Cxcl12", "Vwf", "Mgp", "Rerg", "Gm18194", "Ckm", "Clec11a", 
         "Atp2a1", "Tnnt3", "Vgll2", "Oit3", "Dcn", "Lum", "Tspan8", "Lyz2", 
         "Gli1", "Ndufa4l2", "Ank1", "Hmox1", "Ndrg4", "Tnnc1", "Gdf10", "Trac", 
         "Myh7", "Dmtn", "Ednrb", "Cnn1", "Bmper", "Apoa1", "Isl2", "Loxl1", 
         "Prss35", "Efemp1", "Slc4a1", "Abca8a", "Serpinb6b", "Edn1", "Gm48398", 
         "Colec11", "Ahnak2", "Nov", "Sox10", "Igfbp6", "Masp1", "Robo2", "Chodl", 
         "Tmem181b-ps", "Msln", "C3", "Dsc3", "Sh3rf2", "Onecut2", "Mpeg1", 
         "Anxa1", "Prkg1", "Acta2")
outFile <- paste(output_folder,
                 "/int_appr3.UMAP.spleenE15.5_pancreasE14.5_integrated.goi.orig.ident.pdf", 
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

# Visualize marker genes as FeaturePlot
goi <- c("Tlx1", "Barx1", "Pdx1", "Hoxa2", "Cpa1", "Epcam", "Upk3b", "Nkx2-5", 
         "Lum", "Col6a2", "Cdkn1c", "Cdk1")
outFile <- paste(output_folder,
                 "/int_appr3.UMAP.spleenE15.5_pancreasE14.5_integrated.goi-markers.orig.ident.pdf", 
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

############## Test for mutual exclusivity of candidate expression #############

## Is Tlx1 co-expressed or inversely correlated with Sox10 (TF that is overexpressed in Fgf9 null pancreas)?
# 1. Define gene expression thresholds
# Extract expression data for Tlx1 and Sox10
Tlx1_expr <- FetchData(so_spleenE15.5_pancreasE14.5_integrated, vars = "Tlx1")
Sox10_expr <- FetchData(so_spleenE15.5_pancreasE14.5_integrated, vars = "Sox10")

# Define a threshold for expression (e.g., greater than 0 means expressed)
Tlx1_expressed <- Tlx1_expr > 0
Sox10_expressed <- Sox10_expr > 0

# 2. Identify cells where only one gene is expressed
# Identify cells where Tlx1 is expressed but Sox10 is not, and vice versa
mutually_exclusive_Tlx1 <- Tlx1_expressed & !Sox10_expressed
mutually_exclusive_Sox10 <- Sox10_expressed & !Tlx1_expressed
both_expressed <- Tlx1_expressed & Sox10_expressed

# Count the number of mutually exclusive cells for each gene
mutually_exclusive_Tlx1_cells <- sum(mutually_exclusive_Tlx1)
mutually_exclusive_Sox10_cells <- sum(mutually_exclusive_Sox10)

# Count the number of cells where both genes are expressed
both_expressed_cells <- sum(both_expressed)

# Print the results
cat("Number of cells where Tlx1 is expressed but Sox10 is not: ", mutually_exclusive_Tlx1_cells, "\n")
cat("Number of cells where Sox10 is expressed but Tlx1 is not: ", mutually_exclusive_Sox10_cells, "\n")
cat("Number of cells where both Tlx1 and Sox10 are expressed: ", both_expressed_cells, "\n")

# 3. Statistical test for mutual exclusivity
# Create a contingency table for the co-expression of Tlx1 and Sox10
contingency_table <- table(Tlx1_expressed, Sox10_expressed)

# Perform Fisher's Exact Test (for small numbers) or Chi-square test
# Fisher's exact test is recommended for small sample sizes (less than 5 in any cell)
fisher_test_result <- fisher.test(contingency_table)

# Print the p-value from the Fisher's test
cat("P-value for mutual exclusivity (Fisher's Exact Test): ", fisher_test_result$p.value, "\n")

# 4. Visualizing the results
# Create a new metadata column to label cells as mutually exclusive for each gene
so_spleenE15.5_pancreasE14.5_integrated$MutualExclusive_Tlx1 <- mutually_exclusive_Tlx1
so_spleenE15.5_pancreasE14.5_integrated$MutualExclusive_Sox10 <- mutually_exclusive_Sox10
so_spleenE15.5_pancreasE14.5_integrated$Both_Expressed <- both_expressed

# Visualize mutual exclusivity of Tlx1 and Sox10 using UMAP
library(patchwork)

# Generate individual FeaturePlots for each feature
p1 <- FeaturePlot(so_spleenE15.5_pancreasE14.5_integrated, 
                  features = "MutualExclusive_Tlx1", 
                  reduction = "umap.integrated",
                  split.by = "orig.ident")
p2 <- FeaturePlot(so_spleenE15.5_pancreasE14.5_integrated, 
                  features = "MutualExclusive_Sox10", 
                  reduction = "umap.integrated",
                  split.by = "orig.ident")
p3 <- FeaturePlot(so_spleenE15.5_pancreasE14.5_integrated, 
                  features = "Both_Expressed", 
                  reduction = "umap.integrated",
                  split.by = "orig.ident")

# Save each plot separately
# Output for Tlx1 expression
outFile_Tlx1 <- paste(output_folder, "/int_appr3.UMAP.MutualExclusive_Tlx1(noSox10).pdf", sep = "")
pdf(outFile_Tlx1, width = 10, height = 5)  # Adjust the height/width for your preference
print(p1)  # Print the first plot
dev.off()

# Output for Sox10 expression
outFile_Sox10 <- paste(output_folder, "/int_appr3.UMAP.MutualExclusive_Sox10(noTlx1).pdf", sep = "")
pdf(outFile_Sox10, width = 10, height = 5)  # Adjust the height/width for your preference
print(p2)  # Print the second plot
dev.off()

# Output for Both Tlx1 and Sox10 expression
outFile_Both <- paste(output_folder, "/int_appr3.UMAP.Both_Expressed(Tlx1+Sox10).pdf", sep = "")
pdf(outFile_Both, width = 10, height = 5)  # Adjust the height/width for your preference
print(p3)  # Print the third plot
dev.off()

# Visualize mutual exclusivity with a Venn diagram (with two overlapping circles)
library(grid)
library(futile.logger)
library(VennDiagram)

# Define the output file path for the Venn diagram
outFile_venn <- paste(output_folder, "/int_appr3.VennDiagram.MutualExclusiveExpression_Tlx1+Sox10.pdf", sep = "")

# Open a PDF device to save the plot
pdf(outFile_venn, width = 10, height = 5)

# Create the Venn diagram with two circles
venn.plot <- venn.diagram(
  x = list(
    "Tlx1 expressed" = which(Tlx1_expressed),
    "Sox10 expressed" = which(Sox10_expressed)
  ),
  category.names = c("Tlx1", "Sox10"),
  filename = NULL,  # Output as an object
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  
  # Center the category names and adjust positioning
  cat.pos = c(0, 180),  # Place the categories on opposite sides (180 degrees apart)
  cat.dist = c(0.05, 0.05)  # Adjust distance of the category labels from the circles
  
)

# Draw the Venn diagram
grid.draw(venn.plot)

# Close the PDF device
dev.off()

## Is Tlx1 co-expressed or inversely correlated with Klf4 (TF that is overexpressed in Fgf9 null pancreas)?
# 1. Define gene expression thresholds
# Extract expression data for Tlx1 and Klf4
Tlx1_expr <- FetchData(so_spleenE15.5_pancreasE14.5_integrated, vars = "Tlx1")
Klf4_expr <- FetchData(so_spleenE15.5_pancreasE14.5_integrated, vars = "Klf4")

# Define a threshold for expression (e.g., greater than 0 means expressed)
Tlx1_expressed <- Tlx1_expr > 0
Klf4_expressed <- Klf4_expr > 0

# 2. Identify cells where only one gene is expressed
# Identify cells where Tlx1 is expressed but Klf4 is not, and vice versa
mutually_exclusive_Tlx1 <- Tlx1_expressed & !Klf4_expressed
mutually_exclusive_Klf4 <- Klf4_expressed & !Tlx1_expressed
both_expressed <- Tlx1_expressed & Klf4_expressed

# Count the number of mutually exclusive cells for each gene
mutually_exclusive_Tlx1_cells <- sum(mutually_exclusive_Tlx1)
mutually_exclusive_Klf4_cells <- sum(mutually_exclusive_Klf4)
both_expressed_cells <- sum(both_expressed)

# Print the results
cat("Number of cells where Tlx1 is expressed but Klf4 is not: ", mutually_exclusive_Tlx1_cells, "\n")
cat("Number of cells where Klf4 is expressed but Tlx1 is not: ", mutually_exclusive_Klf4_cells, "\n")
cat("Number of cells where both Tlx1 and Klf4 are expressed: ", both_expressed_cells, "\n")

# 3. Statistical test for mutual exclusivity
# Create a contingency table for the co-expression of Tlx1 and Klf4
contingency_table <- table(Tlx1_expressed, Klf4_expressed)

# Perform Fisher's Exact Test (for small numbers) or Chi-square test
# Fisher's exact test is recommended for small sample sizes (less than 5 in any cell)
fisher_test_result <- fisher.test(contingency_table)

# Print the p-value from the Fisher's test
cat("P-value for mutual exclusivity (Fisher's Exact Test): ", fisher_test_result$p.value, "\n")

# 4. Visualizing the results
# Create a new metadata column to label cells as mutually exclusive for each gene
so_spleenE15.5_pancreasE14.5_integrated$MutualExclusive_Tlx1 <- mutually_exclusive_Tlx1
so_spleenE15.5_pancreasE14.5_integrated$MutualExclusive_Klf4 <- mutually_exclusive_Klf4
so_spleenE15.5_pancreasE14.5_integrated$Both_Expressed <- both_expressed

# Visualize mutual exclusivity of Tlx1 and Klf4 using UMAP
library(patchwork)

# Generate individual FeaturePlots for each feature
p1 <- FeaturePlot(so_spleenE15.5_pancreasE14.5_integrated, 
                  features = "MutualExclusive_Tlx1", 
                  reduction = "umap.integrated",
                  split.by = "orig.ident")
p2 <- FeaturePlot(so_spleenE15.5_pancreasE14.5_integrated, 
                  features = "MutualExclusive_Klf4", 
                  reduction = "umap.integrated",
                  split.by = "orig.ident")
p3 <- FeaturePlot(so_spleenE15.5_pancreasE14.5_integrated, 
                  features = "Both_Expressed", 
                  reduction = "umap.integrated",
                  split.by = "orig.ident")

# Save each plot separately
# Output for Tlx1 expression
outFile_Tlx1 <- paste(output_folder, "/int_appr3.UMAP.MutualExclusive_Tlx1(noKlf4).pdf", sep = "")
pdf(outFile_Tlx1, width = 10, height = 5)  # Adjust the height/width for your preference
print(p1)  # Print the first plot
dev.off()

# Output for Klf4 expression
outFile_Klf4 <- paste(output_folder, "/int_appr3.UMAP.MutualExclusive_Klf4(noTlx1).pdf", sep = "")
pdf(outFile_Klf4, width = 10, height = 5)  # Adjust the height/width for your preference
print(p2)  # Print the second plot
dev.off()

# Output for Both Tlx1 and Klf4 expression
outFile_Both <- paste(output_folder, "/int_appr3.UMAP.Both_Expressed(Tlx1+Klf4).pdf", sep = "")
pdf(outFile_Both, width = 10, height = 5)  # Adjust the height/width for your preference
print(p3)  # Print the third plot
dev.off()

# Visualize mutual exclusivity with a Venn diagram (with two overlapping circles)
library(grid)
library(futile.logger)
library(VennDiagram)

# Define the output file path for the Venn diagram
outFile_venn <- paste(output_folder, "/int_appr3.VennDiagram.MutualExclusiveExpression_Tlx1+Klf4.pdf", sep = "")

# Open a PDF device to save the plot
pdf(outFile_venn, width = 10, height = 5)

# Create the Venn diagram with two circles
venn.plot <- venn.diagram(
  x = list(
    "Tlx1 expressed" = which(Tlx1_expressed),
    "Klf4 expressed" = which(Klf4_expressed)
  ),
  category.names = c("Tlx1", "Klf4"),
  filename = NULL,  # Output as an object
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  
  # Center the category names and adjust positioning
  cat.pos = c(0, 180),  # Place the categories on opposite sides (180 degrees apart)
  cat.dist = c(0.05, 0.05)  # Adjust distance of the category labels from the circles
  
)

# Draw the Venn diagram
grid.draw(venn.plot)

# Close the PDF device
dev.off()

## Is Tlx1 co-expressed or inversely correlated with Barx1
# 1. Define gene expression thresholds
# Extract expression data for Tlx1 and Barx1
Tlx1_expr <- FetchData(so_spleenE15.5_pancreasE14.5_integrated, vars = "Tlx1")
Barx1_expr <- FetchData(so_spleenE15.5_pancreasE14.5_integrated, vars = "Barx1")

# Define a threshold for expression (e.g., greater than 0 means expressed)
Tlx1_expressed <- Tlx1_expr > 0
Barx1_expressed <- Barx1_expr > 0

# 2. Identify cells where only one gene is expressed
# Identify cells where Tlx1 is expressed but Barx1 is not, and vice versa
mutually_exclusive_Tlx1 <- Tlx1_expressed & !Barx1_expressed
mutually_exclusive_Barx1 <- Barx1_expressed & !Tlx1_expressed
both_expressed <- Tlx1_expressed & Barx1_expressed

# Count the number of mutually exclusive cells for each gene
mutually_exclusive_Tlx1_cells <- sum(mutually_exclusive_Tlx1)
mutually_exclusive_Barx1_cells <- sum(mutually_exclusive_Barx1)
both_expressed_cells <- sum(both_expressed)

# Print the results
cat("Number of cells where Tlx1 is expressed but Barx1 is not: ", mutually_exclusive_Tlx1_cells, "\n")
cat("Number of cells where Barx1 is expressed but Tlx1 is not: ", mutually_exclusive_Barx1_cells, "\n")
cat("Number of cells where both Tlx1 and Barx1 are expressed: ", both_expressed_cells, "\n")

# 3. Statistical test for mutual exclusivity
# Create a contingency table for the co-expression of Tlx1 and Barx1
contingency_table <- table(Tlx1_expressed, Barx1_expressed)

# Perform Fisher's Exact Test (for small numbers) or Chi-square test
# Fisher's exact test is recommended for small sample sizes (less than 5 in any cell)
fisher_test_result <- fisher.test(contingency_table)

# Print the p-value from the Fisher's test
cat("P-value for mutual exclusivity (Fisher's Exact Test): ", fisher_test_result$p.value, "\n")

# 4. Visualizing the results
# Create a new metadata column to label cells as mutually exclusive for each gene
so_spleenE15.5_pancreasE14.5_integrated$MutualExclusive_Tlx1 <- mutually_exclusive_Tlx1
so_spleenE15.5_pancreasE14.5_integrated$MutualExclusive_Barx1 <- mutually_exclusive_Barx1
so_spleenE15.5_pancreasE14.5_integrated$Both_Expressed <- both_expressed

# Visualize mutual exclusivity of Tlx1 and Barx1 using UMAP
library(patchwork)

# Generate individual FeaturePlots for each feature
p1 <- FeaturePlot(so_spleenE15.5_pancreasE14.5_integrated, 
                  features = "MutualExclusive_Tlx1", 
                  reduction = "umap.integrated",
                  split.by = "orig.ident")
p2 <- FeaturePlot(so_spleenE15.5_pancreasE14.5_integrated, 
                  features = "MutualExclusive_Barx1", 
                  reduction = "umap.integrated",
                  split.by = "orig.ident")
p3 <- FeaturePlot(so_spleenE15.5_pancreasE14.5_integrated, 
                  features = "Both_Expressed", 
                  reduction = "umap.integrated",
                  split.by = "orig.ident")

# Save each plot separately
# Output for Tlx1 expression
outFile_Tlx1 <- paste(output_folder, "/int_appr3.UMAP.MutualExclusive_Tlx1(noBarx1).pdf", sep = "")
pdf(outFile_Tlx1, width = 10, height = 5)  # Adjust the height/width for your preference
print(p1)  # Print the first plot
dev.off()

# Output for Barx1 expression
outFile_Barx1 <- paste(output_folder, "/int_appr3.UMAP.MutualExclusive_Barx1(noTlx1).pdf", sep = "")
pdf(outFile_Barx1, width = 10, height = 5)  # Adjust the height/width for your preference
print(p2)  # Print the second plot
dev.off()

# Output for Both Tlx1 and Barx1 expression
outFile_Both <- paste(output_folder, "/int_appr3.UMAP.Both_Expressed(Tlx1+Barx1).pdf", sep = "")
pdf(outFile_Both, width = 10, height = 5)  # Adjust the height/width for your preference
print(p3)  # Print the third plot
dev.off()

# Visualize mutual exclusivity with a Venn diagram (with two overlapping circles)
library(grid)
library(futile.logger)
library(VennDiagram)

# Define the output file path for the Venn diagram
outFile_venn <- paste(output_folder, "/int_appr3.VennDiagram.MutualExclusiveExpression_Tlx1+Barx1.pdf", sep = "")

# Open a PDF device to save the plot
pdf(outFile_venn, width = 10, height = 5)

# Create the Venn diagram with two circles
venn.plot <- venn.diagram(
  x = list(
    "Tlx1 expressed" = which(Tlx1_expressed),
    "Barx1 expressed" = which(Barx1_expressed)
  ),
  category.names = c("Tlx1", "Barx1"),
  filename = NULL,  # Output as an object
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  
  # Center the category names and adjust positioning
  cat.pos = c(0, 180),  # Place the categories on opposite sides (180 degrees apart)
  cat.dist = c(0.05, 0.05)  # Adjust distance of the category labels from the circles
  
)

# Draw the Venn diagram
grid.draw(venn.plot)

# Close the PDF device
dev.off()


##### Comparative analysis of additional genes 
## Nr2f2, Tlx1, Barx1

# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)

# Load Seurat object
so_spleenE15.5_pancreasE14.5_integrated <- readRDS("/Users/veralaub/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/integrated_analysis_spleen+pancreas/int_appr3.so_spleenE15.5_pancreasE14.5_integrated.rds")

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/integrated_analysis_spleen+pancreas/"

## Is Tlx1 co-expressed or inversely correlated with Nr2f2
# 1. Define gene expression thresholds
# Extract expression data for Tlx1 and Nr2f2
Tlx1_expr <- FetchData(so_spleenE15.5_pancreasE14.5_integrated, vars = "Tlx1")
Nr2f2_expr <- FetchData(so_spleenE15.5_pancreasE14.5_integrated, vars = "Nr2f2")

# Define a threshold for expression (e.g., greater than 0 means expressed)
Tlx1_expressed <- Tlx1_expr > 0
Nr2f2_expressed <- Nr2f2_expr > 0

# 2. Identify cells where only one gene is expressed
# Identify cells where Tlx1 is expressed but Nr2f2 is not, and vice versa
mutually_exclusive_Tlx1 <- Tlx1_expressed & !Nr2f2_expressed
mutually_exclusive_Nr2f2 <- Nr2f2_expressed & !Tlx1_expressed
both_expressed <- Tlx1_expressed & Nr2f2_expressed

# Count the number of mutually exclusive cells for each gene
mutually_exclusive_Tlx1_cells <- sum(mutually_exclusive_Tlx1)
mutually_exclusive_Nr2f2_cells <- sum(mutually_exclusive_Nr2f2)
both_expressed_cells <- sum(both_expressed)

# Print the results
cat("Number of cells where Tlx1 is expressed but Nr2f2 is not: ", mutually_exclusive_Tlx1_cells, "\n")
cat("Number of cells where Nr2f2 is expressed but Tlx1 is not: ", mutually_exclusive_Nr2f2_cells, "\n")
cat("Number of cells where both Tlx1 and Nr2f2 are expressed: ", both_expressed_cells, "\n")

# 3. Statistical test for mutual exclusivity
# Create a contingency table for the co-expression of Tlx1 and Nr2f2
contingency_table <- table(Tlx1_expressed, Nr2f2_expressed)

# Perform Fisher's Exact Test (for small numbers) or Chi-square test
# Fisher's exact test is recommended for small sample sizes (less than 5 in any cell)
fisher_test_result <- fisher.test(contingency_table)

# Print the p-value from the Fisher's test
cat("P-value for mutual exclusivity (Fisher's Exact Test): ", fisher_test_result$p.value, "\n")

# 4. Visualizing the results
# Create a new metadata column to label cells as mutually exclusive for each gene
so_spleenE15.5_pancreasE14.5_integrated$MutualExclusive_Tlx1 <- mutually_exclusive_Tlx1
so_spleenE15.5_pancreasE14.5_integrated$MutualExclusive_Nr2f2 <- mutually_exclusive_Nr2f2
so_spleenE15.5_pancreasE14.5_integrated$Both_Expressed <- both_expressed

# Generate individual FeaturePlots for each feature
p1 <- FeaturePlot(so_spleenE15.5_pancreasE14.5_integrated, 
                  features = "MutualExclusive_Tlx1", 
                  reduction = "umap.integrated",
                  split.by = "orig.ident")
p2 <- FeaturePlot(so_spleenE15.5_pancreasE14.5_integrated, 
                  features = "MutualExclusive_Nr2f2", 
                  reduction = "umap.integrated",
                  split.by = "orig.ident")
p3 <- FeaturePlot(so_spleenE15.5_pancreasE14.5_integrated, 
                  features = "Both_Expressed", 
                  reduction = "umap.integrated",
                  split.by = "orig.ident")

# Save each plot separately
# Output for Tlx1 expression
outFile_Tlx1 <- paste(output_folder, "/int_appr3.UMAP.MutualExclusive_Tlx1(noNr2f2).pdf", sep = "")
pdf(outFile_Tlx1, width = 10, height = 5)  # Adjust the height/width for your preference
print(p1)  # Print the first plot
dev.off()

# Output for Nr2f2 expression
outFile_Nr2f2 <- paste(output_folder, "/int_appr3.UMAP.MutualExclusive_Nr2f2(noTlx1).pdf", sep = "")
pdf(outFile_Nr2f2, width = 10, height = 5)  # Adjust the height/width for your preference
print(p2)  # Print the second plot
dev.off()

# Output for Both Tlx1 and Nr2f2 expression
outFile_Both <- paste(output_folder, "/int_appr3.UMAP.Both_Expressed(Tlx1+Nr2f2).pdf", sep = "")
pdf(outFile_Both, width = 10, height = 5)  # Adjust the height/width for your preference
print(p3)  # Print the third plot
dev.off()

# Visualize mutual exclusivity with a Venn diagram (with two overlapping circles)
library(grid)
library(futile.logger)
library(VennDiagram)

# Define the output file path for the Venn diagram
outFile_venn <- paste(output_folder, "/int_appr3.VennDiagram.MutualExclusiveExpression_Tlx1+Nr2f2.pdf", sep = "")

# Open a PDF device to save the plot
pdf(outFile_venn, width = 10, height = 5)

# Create the Venn diagram with two circles
venn.plot <- venn.diagram(
  x = list(
    "Tlx1 expressed" = which(Tlx1_expressed),
    "Nr2f2 expressed" = which(Nr2f2_expressed)
  ),
  category.names = c("Tlx1", "Nr2f2"),
  filename = NULL,  # Output as an object
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  
  # Center the category names and adjust positioning
  cat.pos = c(0, 180),  # Place the categories on opposite sides (180 degrees apart)
  cat.dist = c(0.05, 0.05)  # Adjust distance of the category labels from the circles
)

# Draw the Venn diagram
grid.draw(venn.plot)

# Close the PDF device
dev.off()
