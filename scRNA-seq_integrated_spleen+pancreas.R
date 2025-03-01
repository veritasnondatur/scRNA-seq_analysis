##### Comparative analysis of integrated pancreas_E14.5 and spleen_E15.5 datasets

# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)

# Load Seurat object
so_spleenE15.5_pancreasE14.5_integrated_old <- readRDS("/Users/veralaub/Documents/postdoc/collaboration/Maurizio/comparative_analysis_spleen+pancreas/integrated_analysis_spleen+pancreas/int_appr3.so_spleenE15.5_pancreasE14.5_integrated.rds")

# Define output folder (for results)
output_folder <- "~/Documents/postdoc/collaboration/Maurizio/WIP_scRNA-seq_integrated_spleen+pancreas/results/"

# Definition of GOIs
goi <- c("Tlx1", "Barx1", "Klf4", "Nr2f2", "Nkx2-5",
         "Fgf9", "Fgfr1", "Fgfr2", "Fgfr3", 
         "Cdk1", "Cdkn1c", "Mki67",
         "Epcam", "Upk3b", "Lum")
outFile <- paste(output_folder,
                 "/int_appr3.UMAP.spleenE15.5_pancreasE14.5_integrated.goi.orig.ident.pdf", 
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

