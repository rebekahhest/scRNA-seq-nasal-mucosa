## Libraries ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)

## Load data ----
seurat_obj <- readRDS("seurat_ass4.rds")

## QC ----
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
VlnPlot(seurat_obj, features = "nFeature_RNA", raster = FALSE, pt.size = 0)
VlnPlot(seurat_obj, features = "nCount_RNA", raster = FALSE, pt.size = 0)
VlnPlot(seurat_obj, features = "percent.mt", raster = FALSE, pt.size = 0)

seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & percent.mt < 10)

# Set system memory to 32GB
Sys.setenv(R_MAX_VSIZE = "32Gb")

## Normalize ----
seurat_obj <- NormalizeData(seurat_obj)
# seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", method = "glmGamPoi") # too computationally expensive

## Feature selection
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scale
seurat_obj <- ScaleData(seurat_obj)

## PCA ----
seurat_obj <- RunPCA(seurat_obj)

# Elbow plot
ElbowPlot(seurat_obj, ndims = 50)

## Clustering ----
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5) # well-separated clusters
# seurat_obj <- FindClusters(seurat_obj, resolution = 0.8) # over fragmented

# UMAP with clusters displayed
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

## Batch Effect Correction ----
# Check for batch effects
DimPlot(seurat_obj, group.by = "mouse_id") # by mouse ID
DimPlot(seurat_obj, group.by = "time") # by time
DimPlot(seurat_obj, group.by = "disease__ontology_label") # by disease
# no distinct clustering observed so no batch effect correction applied

## Manual Annotation ----
# Find cluster markers
markers <- FindAllMarkers(seurat_obj,
                          features = VariableFeatures(seurat_obj),
                          only.pos = TRUE,
                          logfc.threshold = 0.25)

# Retrieve top genes per cluster by average log2FC
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Feature plot of cluster 4 (consistent cell type)
FeaturePlot(seurat_obj, features = c("Iglc2",
                                     "Fcmr",
                                     "Iglc1",
                                     "Ighd",
                                     "Ms4a1"))

# Feature plot of cluster 33 (mixed cell types)
FeaturePlot(seurat_obj, features = c("Sult1e1",
                                     "Tac4",
                                     "Hist1h2ap",
                                     "Col9a1",
                                     "Svopl"))

# Feature plot across clusters of different cell types
FeaturePlot(seurat_obj, features = c("Krt15",
                                     "Ms4a7",
                                     "Ly6g",
                                     "Ms4a1",
                                     "Ptprb"))

# Review current labels
levels(seurat_obj)

# Identify cluster labels
new_labels <- c(
  "Olfactory receptors",    # 0
  "Olfactory receptors",    # 1
  "Macrophages",            # 2
  "Basal epithelial",       # 3
  "B lymphocytes",          # 4
  "Neurons",                # 5
  "Macrophages",            # 6
  "Endothelial",            # 7
  "T/NK",                   # 8
  "Neutrophils/myeloid",    # 9
  "Epithelial",             # 10
  "Olfactory neurons",      # 11
  "Fibroblasts",            # 12
  "Ciliated epithelium",    # 13
  "Dendritic",              # 14
  "Secretory epithelial",   # 15
  "Secretory epithelial",   # 16
  "Neutrophils",            # 17
  "Progenitor epithelial",  # 18
  "Smooth muscle/pericytes",# 19
  "Osteoblasts",            # 20
  "Secretory epithelial",   # 21
  "Neutrophils",            # 22
  "Secretory epithelial",   # 23
  "Secretory epithelial",   # 24
  "Tuft",                   # 25
  "Glial",                  # 26
  "Proliferating cells",    # 27
  "Proliferating cells",    # 28
  "Proliferating cells",    # 29
  "Proliferating cells",    # 30
  "Proliferating cells",    # 31
  "Cancer",                 # 32
  "Secretory epithelial"    # 33
)

# Assign cluster names
names(new_labels) <- levels(seurat_obj)

# Apply labels
seurat_obj <- RenameIdents(seurat_obj, new_labels)

# Plot UMAP with labels
DimPlot(seurat_obj, label = TRUE, repel = TRUE) + NoLegend()

## Differential Expression Analysis ----
# Review for NA/blank data
# Time
unique(seurat_obj$time) # no empty or NA

# Condition
unique(seurat_obj$disease__ontology_label) # no empty or NA

# Tissue
unique(seurat_obj$organ_custom) # no empty or NA

# Mouse ID
unique(seurat_obj$mouse_id) # empty mouse id ""
# Convert empty mouse id to NA
seurat_obj$mouse_id[seurat_obj$mouse_id == ""] <- NA
# Check NA count
sum(is.na(seurat_obj$mouse_id)) # 54055
# contaminated cDNA --> remove samples for statistical analyses

# Subset only samples with existing Mouse ID
seurat_de <- subset(seurat_obj, subset = !is.na(mouse_id))

# Subset only Macrophage clusters
seurat_mac <- subset(seurat_de, idents = "Macrophages")

seurat_mac$mouse <- gsub("(m[0-9])_.*", "\\1", seurat_mac$mouse_id)       # m1, m2, m3
seurat_mac$tissue <- gsub("m[0-9]_(.*)_.*", "\\1", seurat_mac$mouse_id)   # ET, RT, Sinus

# Pseudobulk by mouse ID x timepoint
pseudo_seurat_time <- AggregateExpression(
  seurat_mac,
  return.seurat = TRUE,
  group.by = c("mouse", "time")  
)

# Set identity
Idents(pseudo_seurat_time) <- "time"

# Define function for Differential Expression using DESeq2
run_de <- function(obj, t1, t2) {
  
  df <- FindMarkers(
    obj,
    ident.1 = t1,
    ident.2 = t2,
    test.use = "DESeq2"
  )
  
  df$significant <- ifelse(
    df$p_val_adj < 0.05 & abs(df$avg_log2FC) > 0.25,
    ifelse(df$avg_log2FC > 0, "Up", "Down"),
    "Not Sig"
  )
  
  df <- df[!is.na(df$p_val_adj) & !is.na(df$avg_log2FC), ]
}

# Define function to extract count of significant differentially expressed genes
sig_genes <- function(df) {
  df %>% 
  filter(significant %in% c("Up", "Down")) %>% 
  nrow()
}

# Call functions for each time point pairwise comparison
# dea_naive_d02 <- run_de(pseudo_seurat_time, "Naive", "D02")
# dea_naive_d05 <- run_de(pseudo_seurat_time, "Naive", "D05")
# dea_naive_d08 <- run_de(pseudo_seurat_time, "Naive", "D08")
# dea_naive_d14 <- run_de(pseudo_seurat_time, "Naive", "D14")
dea_d02_d05 <- run_de(pseudo_seurat_time, "D02", "D05")
dea_d02_d08 <- run_de(pseudo_seurat_time, "D02", "D08")
dea_d02_d14 <- run_de(pseudo_seurat_time, "D02", "D14")
# dea_d05_d08 <- run_de(pseudo_seurat_time, "D05", "D08")
# dea_d05_d14 <- run_de(pseudo_seurat_time, "D05", "D14")
# dea_d08_d14 <- run_de(pseudo_seurat_time, "D08", "D14")

sig_d02_d05 <- sig_genes(dea_d02_d05) # 594
sig_d02_d08 <- sig_genes(dea_d02_d08) # 667
sig_d02_d14 <- sig_genes(dea_d02_d14) # 267

# Plot
p1 <- ggplot(dea_d02_d05, aes(x=avg_log2FC, y=-log10(p_val_adj), colour = significant)) +
  geom_point(alpha=0.5, size=1) +
  scale_color_manual(values=c("Down"="blue", "Not Sig"="gray", "Up"="red")) +
  labs(x="Log2 Fold Change", y="-Log10 Adjusted p-value", title = "A) D02 vs D05") +
  theme_minimal()

p2 <- ggplot(dea_d02_d08, aes(x=avg_log2FC, y=-log10(p_val_adj), colour = significant)) +
  geom_point(alpha=0.5, size=1) +
  scale_color_manual(values=c("Down"="blue", "Not Sig"="gray", "Up"="red")) +
  labs(x="Log2 Fold Change", y="-Log10 Adjusted p-value", title = "B) D02 vs D08") +
  theme_minimal()

p3 <- ggplot(dea_d02_d14, aes(x=avg_log2FC, y=-log10(p_val_adj), colour = significant)) +
  geom_point(alpha=0.5, size=1) +
  scale_color_manual(values=c("Down"="blue", "Not Sig"="gray", "Up"="red")) +
  labs(x="Log2 Fold Change", y="-Log10 Adjusted p-value", title = "C) D02 vs D14") +
  theme_minimal()

p1 + p2 + p3

# Pseudobulk by mouse ID x tissue
pseudo_seurat_tissue <- AggregateExpression(
  seurat_mac,
  return.seurat = TRUE,
  group.by = c("mouse", "organ_custom")  
)

# Set identity
Idents(pseudo_seurat_tissue) <- "organ_custom"

# Call function for each tissue type pairwise comparison
dea_OM_RM <- run_de(pseudo_seurat_tissue, "OM", "RM")
dea_OM_LNG <- run_de(pseudo_seurat_tissue, "OM", "LNG")
dea_RM_LNG <- run_de(pseudo_seurat_tissue, "RM", "LNG")

sig_OM_RM <- sig_genes(dea_OM_RM) # 832
sig_OM_LNG <- sig_genes(dea_OM_LNG) # 773
sig_RM_LNG <- sig_genes(dea_RM_LNG) # 1114

# Volcano plots
p1 <- ggplot(dea_OM_RM, aes(x=avg_log2FC, y=-log10(p_val_adj), colour = significant)) +
  geom_point(alpha=0.5, size=1) +
  scale_color_manual(values=c("Down"="blue", "Not Sig"="gray", "Up"="red")) +
  labs(x="Log2 Fold Change", y="-Log10 Adjusted p-value", title = "A) OM vs RM") +
  theme_minimal()

p2 <- ggplot(dea_OM_LNG, aes(x=avg_log2FC, y=-log10(p_val_adj), colour = significant)) +
  geom_point(alpha=0.5, size=1) +
  scale_color_manual(values=c("Down"="blue", "Not Sig"="gray", "Up"="red")) +
  labs(x="Log2 Fold Change", y="-Log10 Adjusted p-value", title = "B) OM vs LNG") +
  theme_minimal()

p3 <- ggplot(dea_RM_LNG, aes(x=avg_log2FC, y=-log10(p_val_adj), colour = significant)) +
  geom_point(alpha=0.5, size=1) +
  scale_color_manual(values=c("Down"="blue", "Not Sig"="gray", "Up"="red")) +
  labs(x="Log2 Fold Change", y="-Log10 Adjusted p-value", title = "C) RM vs LNG") +
  theme_minimal()

p1 + p2 + p3

## Functional Enrichment Analysis ----
# Define function to run enrichment analysis (ORA: GO BP)
run_go_enrichment <- function(dea_df) {
  
  # Add gene column from rownames
  dea_df$gene <- rownames(dea_df)
  
  # Extract significant genes
  sig_genes <- dea_df %>%
    filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>%
    pull(gene) %>%
    na.omit() %>%
    unique()
  
  # Define background universe
  all_genes <- dea_df %>%
    pull(gene) %>%
    na.omit() %>%
    unique()
  
  # Run GO enrichment
  ego_bp <- clusterProfiler::enrichGO(
    gene          = sig_genes,
    universe      = all_genes,
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  return(ego_bp)
}

# Call function for each time point pairwise comparison
# res_naive_d02 <- run_go_enrichment(dea_naive_d02)
# res_naive_d05 <- run_go_enrichment(dea_naive_d05)
# res_naive_d08 <- run_go_enrichment(dea_naive_d08)
# res_naive_d14 <- run_go_enrichment(dea_naive_d14)
res_d02_d05 <- run_go_enrichment(dea_d02_d05)
res_d02_d08 <- run_go_enrichment(dea_d02_d08)
res_d02_d14 <- run_go_enrichment(dea_d02_d14)
# res_d05_d08 <- run_go_enrichment(dea_d05_d08)
# res_d05_d14 <- run_go_enrichment(dea_d05_d14)
# res_d08_d14 <- run_go_enrichment(dea_d08_d14)

# Extract number of significant terms
nrow(as.data.frame(res_d02_d05)) # 617
nrow(as.data.frame(res_d02_d08)) # 995
nrow(as.data.frame(res_d02_d14)) # 14

# Plot
p1 <- dotplot(res_d02_d05, showCategory = 10, title = "A) D02 vs D05")
p2 <- dotplot(res_d02_d08, showCategory = 10, title = "B) D02 vs D08")
p3 <- dotplot(res_d02_d14, showCategory = 10, title = "C) D02 vs D14")

p1 + p2 + p3

# Call function for each tissue type pairwise comparison
res_OM_RM <- run_go_enrichment(dea_OM_RM)
res_OM_LNG <- run_go_enrichment(dea_OM_LNG)
res_RM_LNG <- run_go_enrichment(dea_RM_LNG)

# Extract number of significant terms
nrow(as.data.frame(res_OM_RM)) # 987
nrow(as.data.frame(res_OM_LNG)) # 965
nrow(as.data.frame(res_RM_LNG)) # 1231

# Plot
p1 <- dotplot(res_OM_RM, showCategory = 10, title = "A) OM vs RM")
p2 <- dotplot(res_OM_LNG, showCategory = 10, title = "B) OM vs LNG")
p3 <- dotplot(res_RM_LNG, showCategory = 10, title = "C) RM vs LNG")

p1 + p2 + p3