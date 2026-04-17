### Analysis of snRNA-seq human data, CS12-CS20 retrieved from public domain (Khouri-Farah et al., 2026)
# Author: Vera Laub
# Last edited: 04/17/2026

# ================================
# Load libraries
# ================================
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(edgeR)
library(limma)

# ================================
# File paths
# ================================
input_file <- "/Users/veralaub/Documents/postdoc/collaboration/Pauline/Gpr50_project/data_human/human_face_no-neuro_clustering_cellrangerARC-raw_emptyDrops_singlets_finalannot_27Mar.rds"

output_dir <- "/Users/veralaub/Documents/postdoc/collaboration/Pauline/Gpr50_project/results_human_face_Khouri-Farah2026/"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ================================
# Load Seurat object
# ================================
seurat_obj <- readRDS(input_file)

# Quick overview
print(seurat_obj)
dim(seurat_obj)
head(seurat_obj@meta.data)

# ================================
# Check metadata fields
# ================================
colnames(seurat_obj@meta.data)

# Likely relevant:
# - seurat_clusters
# - stage
# - orig.ident

# =====================================================================
# Set active parameters in Seurat objects to those likely used in paper
# =====================================================================
# Check current active clustering
Idents(seurat_obj)
table(seurat_obj$seurat_clusters)

# Check which resolutions are used and how many clusters are in each
sapply(
  grep("RNA_snn_res", colnames(seurat_obj@meta.data), value = TRUE),
  function(res) length(unique(seurat_obj[[res]][,1]))
)

# Check where 8 clusters reported in paper are stored
table(seurat_obj$celltype)
table(seurat_obj$subtype_reduced)

# ================================
# UMAP of entire dataset (8 cell types)
# ================================
pdf(file.path(output_dir, "UMAP_entireDataset_celltypes.pdf"), width = 8, height = 6)

p1 <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE   # avoids overlapping labels
) + ggtitle("UMAP - Cell types (paper)")

print(p1)

dev.off()


# ================================
# UMAP split by developmental stage
# ================================
pdf(file.path(output_dir, "UMAP_entireDataset_celltypes_split_by_stage.pdf"), width = 12, height = 8)

p2 <- DimPlot(
  seurat_obj,
  reduction = "umap",
  split.by = "stage",
  group.by = "celltype",
  label = FALSE,
  ncol = 3
) + ggtitle("UMAP split by stage (cell types)")

print(p2)

dev.off()


# ================================
# Feature expression (UMAP)
# ================================
genes_of_interest <- c("GPR50", "DGKK", "HAND2", "DLX6")

pdf(file.path(output_dir, "FeaturePlots_entireDataset_genes_of_interest.pdf"), width = 10, height = 8)

p3 <- FeaturePlot(
  seurat_obj,
  features = genes_of_interest,
  reduction = "umap",
  ncol = 2,
  order = TRUE,
  cols = c("lightgrey", "red"),
  min.cutoff = "q05",
  max.cutoff = "q95"
)

print(p3)

dev.off()


# ================================
# Violin plots (by celltype instead of clusters)
# ================================
pdf(file.path(output_dir, "ViolinPlots_entireDataset_celltypes.pdf"), width = 10, height = 8)

p4 <- VlnPlot(
  seurat_obj,
  features = genes_of_interest,
  group.by = "celltype",
  pt.size = 0
)

print(p4)

dev.off()


# ================================
# Violin plots split by stage
# ================================
pdf(file.path(output_dir, "ViolinPlots_entireDataset_by_stage.pdf"), width = 12, height = 8)

p5 <- VlnPlot(
  seurat_obj,
  features = genes_of_interest,
  group.by = "stage",
  pt.size = 0
)

print(p5)

dev.off()

# ================================================
# PSEUDOBULK DIFFERENTIAL EXPRESSION
# ================================================
# ================================================
# STEP 1 — subset clusters
# ================================================
clusters_of_interest <- c(
  "mesenchyme",
  "ectoderm",
  "SOX2+_SOX10-"
)

subset_obj <- subset(
  seurat_obj,
  subset = celltype %in% clusters_of_interest
)

# ================================================
# STEP 2 — clean labels (CRITICAL)
# ================================================
subset_obj$celltype_clean <- as.character(subset_obj$celltype)

subset_obj$celltype_clean <- recode(
  subset_obj$celltype_clean,
  "SOX2+_SOX10-" = "SOX2pos_SOX10neg"
)

subset_obj$celltype_clean <- make.names(subset_obj$celltype_clean)

# ================================================
# STEP 3 — pseudobulk aggregation
# ================================================
pb <- AggregateExpression(
  subset_obj,
  group.by = c("celltype_clean", "orig.ident"),
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)

counts <- pb$RNA

# ================================================
# STEP 4 — CLEAN COLUMN NAMES (CRITICAL FIX)
# ================================================
colnames(counts) <- gsub("-", "_", colnames(counts))
colnames(counts) <- gsub("\\+", "pos", colnames(counts))

# ================================================
# STEP 5 — build metadata safely
# ================================================
split_names <- strsplit(colnames(counts), "_")

meta <- data.frame(
  celltype = sapply(split_names, `[`, 1),
  orig.ident = sapply(split_names, `[`, 2)
)

meta$celltype <- factor(meta$celltype)

# ================================================
# SAFETY CHECKS
# ================================================
cat("\nCelltype counts:\n")
print(table(meta$celltype))

if (length(levels(meta$celltype)) < 2) {
  stop("ERROR: fewer than 2 cell types found.")
}

if (!all(colnames(counts) == paste(meta$celltype, meta$orig.ident, sep = "_"))) {
  warning("Column/meta mismatch detected — continuing but check carefully.")
}

# ================================================
# STEP 6 — filter low genes
# ================================================
keep_genes <- rowSums(counts > 0) >= 2
counts <- counts[keep_genes, ]

# ================================================
# STEP 7 — design matrix
# ================================================
design <- model.matrix(~0 + meta$celltype)
colnames(design) <- levels(meta$celltype)

cat("\nDesign matrix columns:\n")
print(colnames(design))

# ================================================
# STEP 8 — edgeR normalization
# ================================================
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)

# ================================================
# STEP 9 — voom transformation
# ================================================
v <- voom(dge, design, plot = FALSE)

# ================================================
# STEP 10 — linear model
# ================================================
fit <- lmFit(v, design)

# ================================================
# STEP 11 — contrasts (FINAL SAFE VERSION)
# ================================================
contrasts <- makeContrasts(
  mes_vs_ect = mesenchyme - ectoderm,
  mes_vs_sox = mesenchyme - SOX2pos,
  ect_vs_sox = ectoderm - SOX2pos,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# ================================================
# STEP 12 — results
# ================================================
res_mes_ect <- topTable(fit2, coef = "mes_vs_ect", number = Inf)
res_mes_sox <- topTable(fit2, coef = "mes_vs_sox", number = Inf)
res_ect_sox <- topTable(fit2, coef = "ect_vs_sox", number = Inf)

# ================================================
# STEP 13 — save results
# ================================================
write.csv(res_mes_ect, file.path(output_dir, "PB_mes_vs_ect.csv"))
write.csv(res_mes_sox, file.path(output_dir, "PB_mes_vs_sox.csv"))
write.csv(res_ect_sox, file.path(output_dir, "PB_ect_vs_sox.csv"))

# ================================================
# STEP 14 — sanity check
# ================================================
grep("GPR50", rownames(res_mes_ect), value = TRUE)


# ================================================
# GPR50 HIGH vs OTHER — PSEUDOBULK PIPELINE
# (memory-safe + fully consistent)
# ================================================
# ================================================
# STEP 1 — DEFINE GPR50 HIGH CELLS
# ================================================
gpr50_vals <- FetchData(seurat_obj, vars = "GPR50")[,1]

cutoff <- quantile(gpr50_vals, 0.95, na.rm = TRUE)

seurat_obj$GPR50_status <- ifelse(
  gpr50_vals > cutoff,
  "GPR50_high",
  "GPR50_other"
)

table(seurat_obj$GPR50_status)

# Set identity
Idents(seurat_obj) <- "GPR50_status"

# ================================================
# STEP 2 — PSEUDOBULK AGGREGATION
# ================================================
pb <- AggregateExpression(
  seurat_obj,
  group.by = c("GPR50_status", "orig.ident"),
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)

counts <- pb$RNA

# ================================================
# STEP 3 — CLEAN COLUMN NAMES (SAFE)
# ================================================
colnames(counts) <- gsub("-", "_", colnames(counts))
colnames(counts) <- gsub("\\+", "pos", colnames(counts))

# ================================================
# STEP 4 — BUILD METADATA (NO STRING SPLITTING)
# ================================================
meta <- data.frame(
  sample = colnames(counts),
  group = ifelse(
    grepl("^GPR50_high", colnames(counts)),
    "GPR50_high",
    "GPR50_other"
  )
)

meta$group <- factor(meta$group)

cat("\nGroup table:\n")
print(table(meta$group))

if (length(levels(meta$group)) < 2) {
  stop("ERROR: fewer than 2 groups detected in pseudobulk.")
}

# ================================================
# STEP 5 — DESIGN MATRIX
# ================================================
design <- model.matrix(~0 + meta$group)
colnames(design) <- levels(meta$group)

cat("\nDesign matrix:\n")
print(colnames(design))

# ================================================
# STEP 6 — FILTER LOWLY EXPRESSED GENES
# ================================================
keep <- rowSums(counts > 0) >= 2
counts <- counts[keep, ]

# ================================================
# STEP 7 — edgeR NORMALIZATION
# ================================================
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)

# ================================================
# STEP 8 — VOOM TRANSFORM
# ================================================
v <- voom(dge, design, plot = FALSE)

# ================================================
# STEP 9 — LINEAR MODEL
# ================================================
fit <- lmFit(v, design)

# ================================================
# STEP 10 — CONTRASTS
# ================================================
contrasts <- makeContrasts(
  GPR50_high_vs_other = GPR50_high - GPR50_other,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# ================================================
# STEP 11 — RESULTS
# ================================================
res <- topTable(fit2, coef = "GPR50_high_vs_other", number = Inf)

# ================================================
# STEP 12 — SAVE OUTPUT
# ================================================
write.csv(
  res,
  file.path(output_dir, "GPR50_high_vs_other_pseudobulk_FINAL.csv")
)

# ================================================
# STEP 13 — SANITY CHECK
# ================================================
head(res)


# ================================================
# GPR50 HIGH vs OTHER — STRICT PSEUDOBULK PIPELINE
# ================================================

library(Seurat)
library(edgeR)
library(limma)
library(dplyr)

# ================================================
# STEP 1 — DEFINE GPR50 HIGH CELLS (QUANTILE)
# ================================================
gpr50_vals <- FetchData(seurat_obj, "GPR50")[,1]

cutoff <- quantile(gpr50_vals, 0.95, na.rm = TRUE)

seurat_obj$GPR50_status <- ifelse(
  gpr50_vals > cutoff,
  "GPR50_high",
  "GPR50_other"
)

table(seurat_obj$GPR50_status)

Idents(seurat_obj) <- "GPR50_status"

# ================================================
# STEP 2 — PSEUDOBULK AGGREGATION
# ================================================
pb <- AggregateExpression(
  seurat_obj,
  group.by = c("GPR50_status", "orig.ident"),
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)

counts <- pb$RNA

# ================================================
# STEP 3 — CLEAN COLUMN NAMES
# ================================================
colnames(counts) <- gsub("-", "_", colnames(counts))
colnames(counts) <- gsub("\\+", "pos", colnames(counts))

# ================================================
# STEP 4 — STRONG FILTERING (IMPORTANT)
# ================================================

# 1. remove very lowly expressed genes
keep1 <- rowSums(counts >= 5) >= 3   # expressed in at least 3 pseudobulk samples

# 2. remove extremely sparse genes
keep2 <- rowSums(counts > 0) >= 10

# 3. keep only genes passing both
keep <- keep1 & keep2

counts <- counts[keep, ]

cat("Genes retained after filtering:", nrow(counts), "\n")

# ================================================
# STEP 5 — META (ROBUST GROUP ASSIGNMENT)
# ================================================
meta <- data.frame(
  sample = colnames(counts),
  group = ifelse(
    grepl("^GPR50_high", colnames(counts)),
    "GPR50_high",
    "GPR50_other"
  )
)

meta$group <- factor(meta$group)

table(meta$group)

if (length(levels(meta$group)) < 2) {
  stop("ERROR: fewer than 2 groups after filtering")
}

# ================================================
# STEP 6 — DESIGN MATRIX
# ================================================
design <- model.matrix(~0 + meta$group)
colnames(design) <- levels(meta$group)

design

# ================================================
# STEP 7 — edgeR NORMALIZATION
# ================================================
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)

# OPTIONAL but recommended for snRNA-seq:
keep_expr <- filterByExpr(dge, design)
dge <- dge[keep_expr, , keep.lib.sizes = FALSE]

# ================================================
# STEP 8 — VOOM TRANSFORMATION
# ================================================
v <- voom(dge, design, plot = FALSE)

# ================================================
# STEP 9 — LINEAR MODEL
# ================================================
fit <- lmFit(v, design)

# ================================================
# STEP 10 — CONTRASTS
# ================================================
contrasts <- makeContrasts(
  GPR50_high_vs_other = GPR50_high - GPR50_other,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# ================================================
# STEP 11 — RESULTS (STRICT FILTER)
# ================================================
res <- topTable(fit2, number = Inf)

# stronger filtering for interpretation
res_filtered <- res %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 0.5)

# ================================================
# STEP 12 — SAVE
# ================================================
output_dir <- "/Users/veralaub/Documents/postdoc/collaboration/Pauline/Gpr50_project/results_human_face_Khouri-Farah2026"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  res_filtered,
  file.path(output_dir, "GPR50_high_vs_other_STRICT_pseudobulk.csv")
)

# ================================================
# STEP 13 — QUICK CHECK
# ================================================
cat("\nSignificant genes:", nrow(res_filtered), "\n")
head(res_filtered)


# =====================================================
# GPR50 HIGH vs OTHER — UP/DOWN REGULATED VISUALIZATION
# =====================================================
DefaultAssay(seurat_obj) <- "RNA"

# ================================================
# STEP 1 — ENSURE GPR50 LABEL EXISTS
# ================================================
gpr50_vals <- FetchData(seurat_obj, "GPR50")[,1]

cutoff <- quantile(gpr50_vals, 0.95, na.rm = TRUE)

seurat_obj$GPR50_status <- ifelse(
  gpr50_vals > cutoff,
  "GPR50_high",
  "GPR50_other"
)

table(seurat_obj$GPR50_status)

# ================================================
# STEP 2 — LOAD PSEUDOBULK RESULTS
# ================================================
markers <- read.csv(
  file.path(output_dir, "GPR50_high_vs_other_STRICT_pseudobulk.csv"),
  row.names = 1
)

markers$gene <- rownames(markers)

# ================================================
# STEP 3 — SPLIT UP / DOWN REGULATED GENES
# ================================================
up_genes <- markers %>%
  filter(logFC > 0, adj.P.Val < 0.05) %>%
  arrange(desc(logFC)) %>%
  pull(gene) %>%
  head(50)

down_genes <- markers %>%
  filter(logFC < 0, adj.P.Val < 0.05) %>%
  arrange(logFC) %>%
  pull(gene) %>%
  head(50)

# ================================================
# STEP 4 — KEEP ONLY GENES PRESENT IN OBJECT
# ================================================
up_genes <- intersect(up_genes, rownames(seurat_obj))
down_genes <- intersect(down_genes, rownames(seurat_obj))

cat("Up genes used:", length(up_genes), "\n")
cat("Down genes used:", length(down_genes), "\n")

# ================================================
# STEP 5 — SCALE DATA (IMPORTANT FOR VISUALIZATION)
# ================================================
all_genes <- unique(c(up_genes, down_genes))

seurat_obj <- ScaleData(
  seurat_obj,
  features = all_genes,
  verbose = FALSE
)

# ================================================
# STEP 6 — UMAP PLOTS (UP REGULATED)
# ================================================
pdf(file.path(output_dir, "GPR50_UP_genes_UMAP.pdf"),
    width = 7, height = 6)

for (g in up_genes) {
  
  if (!g %in% rownames(seurat_obj)) next
  
  p <- FeaturePlot(
    seurat_obj,
    features = g,
    reduction = "umap",
    order = TRUE,
    pt.size = 0.1
  ) + ggtitle(g)
  
  print(p)
}

dev.off()

# ================================================
# STEP 7 — UMAP PLOTS (DOWN REGULATED)
# ================================================
pdf(file.path(output_dir, "GPR50_DOWN_genes_UMAP.pdf"),
    width = 7, height = 6)

for (g in down_genes) {
  
  if (!g %in% rownames(seurat_obj)) next
  
  p <- FeaturePlot(
    seurat_obj,
    features = g,
    reduction = "umap",
    order = TRUE,
    pt.size = 0.1
  ) + ggtitle(g)
  
  print(p)
}

dev.off()

# ================================================
# STEP 8 — VIOLIN PLOTS (UP REGULATED)
# ================================================
pdf(file.path(output_dir, "GPR50_UP_genes_Violin.pdf"),
    width = 7, height = 5)

for (g in up_genes) {
  
  if (!g %in% rownames(seurat_obj)) next
  
  p <- VlnPlot(
    seurat_obj,
    features = g,
    group.by = "GPR50_status",
    pt.size = 0
  ) + ggtitle(g)
  
  print(p)
}

dev.off()

# ================================================
# STEP 9 — VIOLIN PLOTS (DOWN REGULATED)
# ================================================
pdf(file.path(output_dir, "GPR50_DOWN_genes_Violin.pdf"),
    width = 7, height = 5)

for (g in down_genes) {
  
  if (!g %in% rownames(seurat_obj)) next
  
  p <- VlnPlot(
    seurat_obj,
    features = g,
    group.by = "GPR50_status",
    pt.size = 0
  ) + ggtitle(g)
  
  print(p)
}

dev.off()

# ================================
# Save session info (reproducibility)
# ================================
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "session_info.txt"))
