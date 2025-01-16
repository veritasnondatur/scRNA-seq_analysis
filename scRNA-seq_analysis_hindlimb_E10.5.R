## scRNA-seq data analysis of E10.5 hindlimb with Seurat
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

# Read barcodes, features (genes) and matrix file: https://www.youtube.com/watch?v=43Z13DS_emQ&list=PLOLdjuxsfI4N1SdaQQYXGoa5Z93hPxWVY
data <- Read10X(data.dir = '/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_E10.5_hindlimb/data_zip')

# Preview of loaded data
data

# Create Seurat objects: https://www.youtube.com/watch?v=cMT92ZExyAQ&list=PLOLdjuxsfI4N1SdaQQYXGoa5Z93hPxWVY&index=2
data <- CreateSeuratObject(counts = data, 
                           project = "hindlimb_E10.5", 
                           min.cells = 3, min.features = 1000)  # "Cells were filtered to ensure inclusion of only those showing a number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes" (Losa et al., 2023)
class(data)

# View Seurat objects
data
colnames(data)
rownames(data)
#view(data)                                                                        # or click D in Data list in upper right window
#view(data@meta.data)

# Save Seurat object
saveRDS(data, file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_E10.5_hindlimb/analysis/SeuratObject_scRNA-seq_E10.5_hindlimb.RDS")

# Read Seurat Objects
data <- readRDS("/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_E10.5_hindlimb/analysis/SeuratObject_scRNA-seq_E10.5_hindlimb.RDS")


######## Standard pre-processing workflow
# Quality control and selecting cells for further analysis
# Data normalization (NormalizeData)
# Identification of high variability features [feature selection] (FindVariableFeatures)
# Data scaling (ScaleData)

### QC metrics: "nFeature_RNA", "nCount_RNA", "percent.mt"

# Low quality cells or empty droplets often have very few genes (low nFeature_RNA & nCount_RNA)
# Cell doublets or multiples have high values of nFeature_RNA & nCount_RNA
# low-quality/dying cells often have high percentage of mitochondrial genes (percent.mt)

# View QC metrics
view(data@meta.data)
range(data$nFeature_RNA)
range(data$nCount_RNA)

# Store mitochondrial percentage in object meta.data
data <- PercentageFeatureSet(data, 
                             pattern = "^MT-", 
                             col.name = "percent.mt")

#view(data@meta.data)
range(data$percent.mt)


### Selecting cells for further analysis

# Visualize QC metrics as a violin plot
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA",                        # comment: Marta etc. seem to have uploaded some sort of prefiltered dataset, as there are no mitochondrial genes present here
                                           "percent.mt"), ncol = 3)

# Use subset function to retrieve a certain subset of cells
data <- subset(data,                                                            # implement filters according to Losa et al., 2023: Cells with...
               subset = nFeature_RNA >3000 & nFeature_RNA <25000                # "number of total expressed transcripts between 3000 and 25,000, corresponding to at least 1000 expressed genes"
             # & nCount_RNA < 20000                            
               & percent.mt < 10)                                               # "the mass of transcripts derived from the mitochondrial chromosomes representing less than 10% of the total."

# Save Seurat object
saveRDS(data, file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_E10.5_hindlimb/analysis/SeuratObject_scRNA-seq_E10.5_hindlimb_subset.RDS")


### Data normalization (NormalizeData)

# Retrieve merged dataset from previous steps
data <- readRDS(file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_E10.5_hindlimb/analysis/SeuratObject_scRNA-seq_E10.5_hindlimb_subset.RDS")

# Inspect counts for dataset
data@assays[["RNA"]]@layers[["counts"]]@x
range(data@assays[["RNA"]]@layers[["counts"]]@x)

# Normalizing the data
data <- SCTransform(data,                                 # "Data were normalized using scTransform, using the best 5000 features." (Losa et al., 2023)
                    method = "glmGamPoi",                 # This is the default method
                    variable.features.n = 5000, 
                    verbose = TRUE)

                      
# Inspect normalized counts for control (D) and nicotine treatment (N)
data@assays[["SCT"]]@data
range(data@assays[["SCT"]]@data)


### Identification of highly variable features (feature selection)
data <- FindVariableFeatures(data,              # "Variable expressed genes across the single cells under cutoff with an average expression of more than 0.0125 and less than 3, and dispersion of more than 0.5, were detected for down-stream analysis" (Guo et al., 2019)
                             selection.method = "vst",          # You can also use "mean.var.plot" or other methods
                             nfeatures = 5000,                  # Specify the number of features to return
                             clip.min = 0.0125,                 # Clip min value to 0.0125 for average expression
                             clip.max = 3,                      # Clip max value to 3 for average expression
                             dispersion.threshold = 0.2         # Set the minimum dispersion threshold
)


# Create the variable feature plot
data@assays[["SCT"]]@var.features
VariableFeaturePlot(data)

# Data scaling, 2000 variable features are used for scaling                     # ScaleData function shifts expression of each gene, such that mean expression across cells for the gene is zero and then scales the expression of each gene so that the variance across all cells is one
data <- ScaleData(data)                                                         # 5000 identified variable features
data@assays[["SCT"]]@scale.data

# Alternativelym, one can perform Data scaling with all genes (takes longer time)
# all.genes <- rownames(data)
# data <- ScaleData(data, features = all.genes)

# Save Seurat object
saveRDS(data, file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_E10.5_hindlimb/analysis/SeuratObject_scRNA-seq_E10.5_hindlimb_norm_scale.RDS")


### Perform PCA on the scaled data (linear dimension reduction)
data <- readRDS(file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_E10.5_hindlimb/analysis/SeuratObject_scRNA-seq_E10.5_hindlimb_norm_scale.RDS")

# Regress out unwanted variables (mitochondrial percentage and number of UMIs)
data <- ScaleData(data,
                  vars.to.regress = c("percent.mt", "nCount_RNA"))

# Perform PCA
data <- RunPCA(data,                                                            # "Cell clusters were identified by constructing a shared nearest neighbor graph followed by a modularity optimization-based clustering algo- rithm (Leiden algorithm106) using the top 60 principal components as determined by PCA." (Losa et al., 2023)
               npcs = 60,                                                       # Number of principal components to compute
               features = VariableFeatures(object = data),                      # Use variable features for PCA
               ndims.print = 1:5,                                               # Print details for the first 5 PCs
               nfeatures.print = 30                                             # Print details for the top 30 features
               )

# Examine and visualize PCA results a few different ways
# DimPlot(), VizDimReduction() and DimHeatmap
DimPlot(data, reduction = "pca", dims = c(1,10))
DimPlot(data, reduction = "pca", dims = c(1,60))
DimPlot(data, reduction = "pca", dims = c(1,100))


### Determine the 'dimensionalily' of the dataset

# With JackStraw Function, runs much slower than ElbowPlot Function
#data <- JackStraw(data, num.replicate = 100)
#data <- ScoreJackStraw(data, dims = 1:20)
#JackStrawPlot(data, dims = 1:20)

# With ElbowPlot Function, more straight forward to use                         # you should include all PCs up to the ellobw shift (?) fo cell clustering, but you always can try to include higher PCs to see how the cell cluster looks like later
ElbowPlot(data)                                                                 # function uses ndims = 20 as default, since we use npca = 20 above increaseing this above the default is not feasible in this case
ElbowPlot(data, ndims = 10, reduction = "pca")                 

# Clustering of cells                                                           # Seurat uses graph-based approach to cluster cells
data <- FindNeighbors(data, dims = 1:60)                                        # "Cell clusters were identified by constructing a shared nearest neighbor graph followed by a modularity optimization-based clustering algorithm (Leiden algorithm) using the top 60 principal components as determined by PCA." (Losa et al., 2023)
data <- FindClusters(data, resolution = 0.8)                                    # "Clustering was performed at multiple resolutions between 0.2 and 2, and optimal resolution was determined empirically based on the expression of known population markers (resolution = 0.8)." (Losa et al., 2023)

# Run non-linear dimensionality reduction (UMAP/t-SNE)
data <- RunUMAP(data, dims = 1:60)
DimPlot(data, reduction = "umap", label = TRUE, repel = TRUE)

#data <- RunTSNE(object = data, dims = 1:20)     # "We ran t-SNE with the same number of PCs and default parameters to visualize the clustering results." (Guo et al., 2019)
#DimPlot(object = data, reduction = "tsne",label = TRUE, repel = TRUE)

View(data)

# Save Seurat object
saveRDS(data, file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_E10.5_hindlimb/analysis/SeuratObject_scRNA-seq_E10.5_hindlimb_SWF.RDS")



### Visualization as UMAP

# Retrieve dataset, pre-analyzed with standard work flow (SWF) from previous steps to visualize batch effects
data_SWF <- readRDS(file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_E10.5_hindlimb/analysis/SeuratObject_scRNA-seq_E10.5_hindlimb_SWF.RDS")
data_SWF

# UMAP clusters
DimPlot(data_SWF, reduction = "umap",label = TRUE)

# PBX1/2, HAND2 and individual genes
FeaturePlot(data_SWF, features = c("Pbx1", 
                                   "Pbx2", 
                                   "Pbx3",
                                   "Hand2"),                    
            cols = c('lightgray', 'blue'),
            pt.size = 0.01)  # Adjust pt.size to your desired value


# PBX1/2, HAND2 and hindlimb fate determinants (literature)
FeaturePlot(data_SWF, features = c("Pbx1", 
                                   "Pbx2", 
                                   "Hand2", 
                                   "Prrx1", 
                                   "Msx1",
                                   "Msx2",
                                   "Tfap2c",
                                   "Twist1",
                                   "Twist2"),                    
            cols = c('lightgray', 'blue'),
            pt.size = 0.01)  # Adjust pt.size to your desired value

# Define a list of gene sets (co-expressed genes)
gene_set_1 <- list(Coexpression = c("Pbx1", 
                                  "Pbx2", 
                                  "Hand2", 
                                  "Prrx1", 
                                  "Msx1",
                                  "Msx2",
                                  "Tfap2c",
                                  "Twist1",
                                  "Twist2"))

# Add module scores to the Seurat object
data <- AddModuleScore(data_SWF, features = gene_set_1, name = "CoexpressionScore")

# Visualize the module score in UMAP
FeaturePlot(data, features = "CoexpressionScore1", 
            pt.size = 0.5) + 
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = median(data$CoexpressionScore1)) +
  labs(title = "Co-expression in E10.5 hindlimb scRNA-seq", 
       x = "UMAP 1", 
       y = "UMAP 2", 
       color = "Co-expression\n(Pbx1, Pbx2, Hand2,\nPrrx1, Msx1, Msx2,\nTfap2c, Twist1, Twist2)") +  # Custom legend title
  theme_minimal()


# PBX1/2, HAND2 and fate determinants (literature); possibly apical-ectodermal ridge cluster?
FeaturePlot(data_SWF, features = c("Pbx1", 
                                   "Pbx2", 
                                   "Hand2",
                                   "Sp6",
                                   "Sp8",
                                   "Krt8",
                                   "Vwa2",
                                   "Fgf8",
                                   "Bmp2",
                                   "Dlx2",
                                   "Dlx5",
                                   "Sox10"),                    
            cols = c('lightgray', 'blue'),
            pt.size = 0.01)  # Adjust pt.size to your desired value

# Define a list of gene sets (co-expressed genes)
gene_set_2 <- list(Coexpression = c("Pbx1", 
                                    "Pbx2",
                                    "Sp6",
                                    "Sp8",
                                    "Krt8",
                                    "Vwa2",
                                    "Fgf8",
                                    "Bmp2",
                                    "Dlx2",
                                    "Dlx5"))

# Add module scores to the Seurat object
data <- AddModuleScore(data_SWF, features = gene_set_2, name = "CoexpressionScore")

# Visualize the module score in UMAP
FeaturePlot(data, features = "CoexpressionScore1", 
            pt.size = 0.5) + 
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = median(data$CoexpressionScore1)) +
  labs(title = "Co-expression in E10.5 hindlimb scRNA-seq", 
       x = "UMAP 1", 
       y = "UMAP 2", 
       color = "Co-expression\n(Pbx1, Pbx2,\nSp6, Sp8, Krt8, \nVwa2, Fgf8, Bmp2, \nDlx2, Dlx5)") +  # Custom legend title
  theme_minimal()

# Violin plot showing expression of a given target in all clusters
VlnPlot(data_SWF, 
        features = c("Pbx1", "Pbx2", "Pbx3"),            # Gene of interest
        group.by = "seurat_clusters", # Group by clusters (default is 'seurat_clusters' if you haven't customized the metadata)
        pt.size = 0.1,                # Adjust point size for individual cells, or set to 0 to hide points
        cols = NULL,                  # Optional: customize colors, NULL uses default palette
        y.max = NULL                  # Optional: limit the Y-axis range
)

################################################################################
### scRNA-seq data integration
# Serves to remove batch effects between different data collections or experimental manipulations.
DimPlot(data_SWF, reduction = "tsne",label = TRUE, group.by = 'orig.ident')

# Load Seurat object without pre-processing (SWF) for data integration
data <- readRDS(file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_nicotine-treated_EB/analysis/data_subset.RDS")
view(data@meta.data)
view(data)


### DETOUR (based on ChatGPT)
# When integrating datasets with different numbers of features (genes) and cells in Seurat, you can follow these steps to prepare and merge your datasets effectively:
# 1. Identify Common Features: First, identify the common genes between the two datasets. This will ensure that you are integrating only the overlapping features:

# Read Seurat Objects
dataset1 <- readRDS("/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_nicotine-treated_EB/analysis/D.RDS")
dataset2 <- readRDS("/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_nicotine-treated_EB/analysis/N.RDS")

# Get the features from both datasets
features1 <- rownames(dataset1)  # Replace with your actual dataset
features2 <- rownames(dataset2)  # Replace with your actual dataset

# Find common features
common_features <- intersect(features1, features2)

# 2. Subset Each Dataset to Common Features
# Once you have the common features, subset each dataset to keep only these features:
dataset1_subset <- dataset1[common_features, ]
dataset2_subset <- dataset2[common_features, ]
  
# 3. Normalize Each Dataset Independently
# Normalize each dataset independently to remove any systematic differences:
library(Seurat)

# Normalize dataset 1
dataset1_subset <- NormalizeData(dataset1_subset)
dataset1_subset <- FindVariableFeatures(dataset1_subset, selection.method = 'vst', nfeatures = 2000)

# Normalize dataset 2
dataset2_subset <- NormalizeData(dataset2_subset)
dataset2_subset <- FindVariableFeatures(dataset2_subset, selection.method = 'vst', nfeatures = 2000)

# 4. Integrate the Datasets
# Now that both datasets are normalized and contain the same features, you can proceed with integration. Use the IntegrateData function:
# Create a list of datasets
data.list <- list(dataset1_subset, dataset2_subset)

#### DETOUR OVER

# Normalize and identify variable features for each dataset independently
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = 'vst', nfeatures = 2000)
  })

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = data.list)

# Find Integration Anchors (will be used to correct the technical difference between the two datasets)
data.anchors <- FindIntegrationAnchors(object.list = data.list,
                                                      anchor.features = features)

# Perform integration to create an 'integrated' data assay
data.integrated <- IntegrateData(anchorset = data.anchors)

# Save Seurat object
saveRDS(data.integrated, "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_nicotine-treated_EB/analysis/merged_ctrl_integrated.RDS")

### Perform an integrated analysis
DefaultAssay(data.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
data.integrated <- ScaleData(data.integrated, verbose = FALSE)
data.integrated <- RunPCA(data.integrated, npcs = 50, verbose = FALSE)
data.integrated <- FindNeighbors(data.integrated ,
                                                 reduction = "pca", dims = 1:20)
data.integrated <- FindClusters(data.integrated, resolution = 0.8)
data.integrated <- RunTSNE(object = data.integrated, 
                                           reduction = "pca",
                                           dims = 1:20)     # "We ran t-SNE with the same number of PCs and default parameters to visualize the clustering results." (Guo et al., 2019)

# Visualization
DimPlot(data.integrated, reduction = "tsne",label = TRUE, group.by = 'orig.ident')

# Compare
plot1 <- DimPlot(data_SWF, reduction = "tsne",label = TRUE, group.by = 'orig.ident')
plot2 <- DimPlot(data.integrated, reduction = "tsne",label = TRUE, group.by = 'orig.ident')
plot1+plot2

plot1_features <- FeaturePlot(data_SWF, features = c("ZFHX3", 
                                                                     "PBX1"),                    
                              cols = c('lightgray', 'blue'))

plot2_features <- FeaturePlot(data.integrated, features = c("ZFHX3", 
                                                                            "PBX1"),                    
                              cols = c('lightgray', 'blue'))

plot3 <- DimPlot(data.integrated, reduction = "tsne",label = TRUE, group.by = 'orig.ident')
plot1_features+plot2_features+plot3

# Save Seurat object
saveRDS(data.integrated, "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_nicotine-treated_EB/analysis/merged_ctrl_integrated.RDS")

############################# FURTHER ANALYSIS/ DETOURS ##############################


######################################################################################
### Clustering and UMAP with ChatGPT
data <- readRDS(file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_nicotine-treated_EB/analysis/data_normalized_scaled.RDS")

# Step 1: Regress out unwanted variables (mitochondrial percentage and number of UMIs)
data <- ScaleData(data, 
                                  vars.to.regress = c("percent.mt", "nCount_RNA"))

# Step 2: Scale the data (this is done in ScaleData)
# Step 3: Run PCA
data <- RunPCA(data,
                               npcs = 50,
                               features = VariableFeatures(object = data),
                               ndims.print = 1:5,
                               nfeatures.print = 30)

# Step 4: Cluster the cells using the top 20 PCs
data <- FindNeighbors(data, dims = 1:20)
data <- FindClusters(data, resolution = 0.8)

# Optional: Run UMAP or t-SNE for visualization
data <- RunUMAP(data, dims = 1:20)

# Print the clustering results
print(head(data@meta.data))


### Create a UMAP Plot
# Step 1: Ensure you have run the necessary steps before this
# (e.g., ScaleData, RunPCA, FindNeighbors, FindClusters)

# Step 2: Run UMAP using the top 20 principal components
data <- RunUMAP(data, dims = 1:20)

# Step 3: Create a UMAP plot, colored by cluster identities
umap_plot <- DimPlot(data, reduction = "umap", group.by = "ident") +
  ggtitle("UMAP of scRNA-seq Data") +
  theme_minimal()

# Display the plot
print(umap_plot)

#########################################################################################################
### Retrieve expression levels of a specific gene for all cells in Seurat object without data integration
library(ggplot2)

# Load normalized dataset
data <- readRDS(file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_nicotine-treated_EB/analysis/data_normalized_scaled.RDS")

# Specify the gene of interest
gene_of_interest <- "ZFHX3"  # Replace with your gene name
rownames(data)

# Check available genes
available_genes <- rownames(data)
print(gene_of_interest %in% available_genes)

# Extract expression levels for the gene of interest
gene_expression_levels <- FetchData(data, vars = gene_of_interest)

# Check the fetched expression levels
print(head(gene_expression_levels))

# Assuming you have a metadata variable for identity
meta_data <- data@meta.data

# Create the data frame for plotting
data_for_plot <- data.frame(
  expression = gene_expression_levels[[gene_of_interest]],  # Ensure this retrieves the correct column
  identity = meta_data$orig.ident  # Change this to your grouping variable if different
)

# Filter out cells with zero expression
data_for_plot <- data_for_plot[data_for_plot$expression > 0, ]

# Check the structure of data_for_plot
str(data_for_plot)

# Convert identity to a factor
data_for_plot$identity <- as.factor(data_for_plot$identity)

# Split the data by group
group1 <- data_for_plot$expression[data_for_plot$identity == "D"]  # Replace with your actual group names
group2 <- data_for_plot$expression[data_for_plot$identity == "N"]

# Perform the t-test
t_test_result <- t.test(group1, group2)

# Extract the p-value
p_value <- t_test_result$p.value

# Create the boxplot
boxplot <- ggplot(data_for_plot, aes(x = identity, y = expression)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = paste("Expression of", gene_of_interest),
       y = "Expression Level [AU], normalized",
       x = "Treatment Group") +
  scale_x_discrete(labels = c("D" = "Control", "N" = "Nicotine")) +  # Replace D and N with Control and Nicotine
  # Increase label size and place inside the plot
  geom_text(aes(x = 1.5, y = median(expression, na.rm = TRUE),  # Centered between the groups
                label = paste("p-value =", format(p_value, digits = 2))), 
            size = 4,  # Adjust size to match axis labels
            vjust = -2,  # Center vertically
            hjust = 0.5)  # Center horizontally

# Display the plot
print(boxplot)

#########################################################################################################
### Retrieve expression levels of a specific gene for all cells in Seurat object using integrated Seurat object
library(ggplot2)

# Load normalized dataset
data.integrated <- readRDS(file = "/Users/veralaub/Documents/postdoc/bioinformatics/data/scRNA-seq/scRNA-seq_nicotine-treated_EB/analysis/merged_ctrl_integrated.RDS")

# Specify the gene of interest
gene_of_interest <- "RELN"  # Replace with your gene name
rownames(data.integrated)

# Check available genes
available_genes <- rownames(data.integrated)
print(gene_of_interest %in% available_genes)

# Extract expression levels for the gene of interest
gene_expression_levels <- FetchData(data.integrated, vars = gene_of_interest)

# Check the fetched expression levels
print(head(gene_expression_levels))

# Assuming you have a metadata variable for identity
meta_data <- data.integrated@meta.data

# Create the data frame for plotting
data_for_plot <- data.frame(
  expression = gene_expression_levels[[gene_of_interest]],  # Ensure this retrieves the correct column
  identity = meta_data$orig.ident  # Change this to your grouping variable if different
)

# Filter out cells with zero expression
#data_for_plot <- data_for_plot[data_for_plot$expression > 0, ]

# Check the structure of data_for_plot
str(data_for_plot)

# Convert identity to a factor
data_for_plot$identity <- as.factor(data_for_plot$identity)

# Split the data by group
group1 <- data_for_plot$expression[data_for_plot$identity == "D"]  # Replace with your actual group names
group2 <- data_for_plot$expression[data_for_plot$identity == "N"]

# Perform the t-test
t_test_result <- t.test(group1, group2)

# Extract the p-value
p_value <- t_test_result$p.value

# Create the boxplot
boxplot <- ggplot(data_for_plot, aes(x = identity, y = expression)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = paste("Expression of", gene_of_interest),
       y = "Expression Level [AU], normalized",
       x = "Treatment Group") +
  scale_x_discrete(labels = c("D" = "Control", "N" = "Nicotine")) +  # Replace D and N with Control and Nicotine
  # Increase label size and place inside the plot
  geom_text(aes(x = 1.5, y = median(expression, na.rm = TRUE),  # Centered between the groups
                label = paste("p-value =", format(p_value, digits = 2))), 
            size = 4,  # Adjust size to match axis labels
            vjust = -2,  # Center vertically
            hjust = 0.5)  # Center horizontally

# Display the plot
print(boxplot)

#########################################################################################################
## Setup of Seurat

# Enter commands in R (or R studio, if installed)
install.packages('Seurat')

# Load packages developed by other labs that can substantially enhance speed 
# and performance
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 
                                       'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))

# We also recommend installing these additional packages, which are used in our 
# vignettes, and enhance the functionality of Seurat:
  
# Signac: analysis of single-cell chromatin data
# SeuratData: automatically load datasets pre-packaged as Seurat objects
# Azimuth: local annotation of scRNA-seq and scATAC-seq queries across multiple organs and tissues
# SeuratWrappers: enables use of additional integration and differential expression methods

# Install the remotes package
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
install.packages('Signac')
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
remotes::install_github("satijalab/azimuth", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)
