# scRNA-seq-nasal-mucosa
Single-cell RNA-seq analysis of mouse nasal mucosa to characterize cellular heterogeneity and host response to viral infection

**Author:** Rebekah Hest  
**Course:** BINF6110 - Genomic Methods for Bioinformatics

## Assignment 4

This repository contains Assignment 4, a scRNA-seq analysis of mouse nasal mucosa using Seurat, including clustering, manual cell type annotation, differential expression, and pathway enrichment to investigate cellular responses to viral infection

### Introduction

Single-cell RNA sequencing (scRNA-seq) has emerged as a powerful tool for studying gene expression in complex biological systems, offering key advantages over traditional bulk RNA sequencing (RNA-seq). While inexpensive, bulk RNA-seq methods average gene expression across many cells at the tissue level which can mask cell diversity and obscure cell-specific information (Mao et al., 2023). In contrast, scRNA-seq can isolate and capture transcriptional profiles at the individual cell-level, enabling characterization of heterogeneous cells and identification of distinct cell populations (Tzec-Interián, González-Padilla and Góngora-Castillo, 2025). This high-resolution technology can be used for the detection and classification of rare but functionally important cell types with applications in disease diagnosis and treatment (Jovic et al., 2022).
<br><br>
This approach is particularly valuable for investigating virus-host interactions and informing the development of therapeutics and vaccines (Chang et al., 2024). Respiratory viral infections, including Influenza A Virus (IAV), typically enter the body through the mouth or nose, where the nasal epithelium serves as a critical immunological barrier (Denney and Ho, 2018). The nasal epithelium is a complex and heterogeneous tissue composed of multiple major and rare cells and the nasal mucosa is highly dynamic, exhibiting changes in response to environmental and infectious stimuli (Cybulski et al., 2026). Such heterogeneity and plasticity make it essential to study gene expression at single-cell resolution in order to fully capture the range of cellular responses to infection and understand how these responses contribute to disease progression and immune defense.
	However, the advantages of scRNA-seq come with increased analytical complexity. Compared to bulk RNA-seq, scRNA-seq data are high-dimensional, sparse, and subject to technical variability such as dropout events and batch effects (Kharchenko, 2021). As a result, analysis requires a multi-step pipeline including quality control, normalization, dimensionality reduction, clustering, and cell type annotation before meaningful biological comparisons can be made (Chen, Ning and Shi, 2019).
  <br><br>
Within the scRNA-seq framework, several methodological choices exist at each stage of analysis. To reduce technical noise of scRNA-seq data, normalization is an essential step that will benefit downstream analyses. Traditional log-normalization methods, such as those implemented in Seurat, scale gene expression values by sequencing depth and apply a logarithmic transformation, offering a simple and computationally efficient approach (Stuart et al., 2019). In contrast, more advanced methods such as SCTransform use a model-based framework based on regularized negative binomial regression to account for technical noise and variance across genes, often improving clustering and differential expression analyses (Hafemeister and Satija, 2019). Clustering is a key step in scRNA-seq analysis, as it enables the identification of distinct cell populations based on transcriptional similarity. Partition-based methods such as k-means require the number of clusters to be specified a priori, which can be limiting in complex and heterogeneous datasets and can be sensitive to outliers resulting in failures to detect rare cell types (Zhang et al., 2023). In contrast, graph-based approaches (e.g., Seurat’s FindClusters) use cell–cell similarity networks to identify clusters in a data-driven manner, making them better suited for capturing the structure of single-cell data. Finally, differential expression analysis in scRNA-seq can be performed either at the single-cell level or by aggregating counts across cells within clusters (pseudobulk). Single-cell methods, such as MAST, model expression at the level of individual cells but can be sensitive to technical noise and inflated false positive rates. In contrast, pseudobulk approaches using tools such as DESeq2 aggregate counts by sample and cell type, providing more robust statistical inference and improved control of false positives, making them a preferred approach for many scRNA-seq studies (Kalantari-Dehaghi, Ghohabi-Esfahani and Emadi-Baygi, 2025).
<br><br>
In this study, a scRNA-seq dataset from Kazer et al. (2025) of murine nasal mucosa will be analyzed to investigate the spatial and temporal variation in host response across distinct tissue types. To address these objectives, a Seurat-based analysis pipeline including clustering and UMAP visualization, cell type annotation, differential expression analysis, and pathway enrichment analysis will be applied. These methods were selected to balance biological interpretability with statistical rigor, enabling the identification of key cell populations and transcriptional pathways involved in the immune response to viral infection.

### Methods

#### Data Acquisition and Quality Control
Single-cell RNA sequencing data were obtained from Kazer et al. (2025) as a pre-processed Seurat (v.5.4.0) object (`readRDS`) containing gene expression counts and associated metadata. The dataset included cells from murine nasal tissues collected across multiple time points following influenza A virus infection.
<br><br>
Quality control was performed to remove low-quality cells based on visualization of gene count, UMI count, and mitochondrial gene expression distributions using violin plots (`VlnPlot`). The percentage of mitochondrial gene expression was calculated using genes with the prefix “mt-” (`PercentageFeatureSet`), and cells with greater than 10% mitochondrial content (`percent.mt < 10`) were removed. Cells with fewer than 500 detected genes (`nFeature_RNA > 500`) were excluded to remove low-quality cells or empty droplets. 

#### Normalization and Feature Selection
Gene expression data were normalized using the log-normalization method implemented in Suerat (`NormalizeData`). Although SCTransform provides a model-based approach that accounts for technical variation, it was not used due to computational constraints associated with the size of the dataset.
<br><br>
Highly variable genes were identified (`FindVariableFeatures`) using the variance-stabilizing transformation (vst) method, where the top 2000 features were selected to capture genes that drive cell variability and reduce noise. The data were then scaled (`ScaleData`) to equalize gene expression values for downstream analyses.

#### Dimensionality Reduction and Clustering
Principal component analysis (PCA) was performed to reduce dimensionality and highlight similarity patterns in the data (`runPCA`). An elbow plot was generated to determine the number of principal components used for downstream analysis, with the first 15 components selected as they captured the majority of variance prior to plateauing (`ElbowPlot`).
<br><br>
Graph-based clustering was performed using the Seurat `FindNeighbors` and `FindClusters` functions based on the selected principal components. Multiple clustering resolutions were evaluated, and a resolution of 0.5 was selected as it produced well-separated clusters without the over-fragmentation observed at higher resolutions (e.g., 0.8). Uniform Manifold Approximation and Projection (UMAP) was used to efficiently visualize high-dimensional single-cell data into interpretable clusters (`RunUMAP`, `DimPlot`).

#### Batch Effect Correction
Potential batch effects were assessed by visualizing UMAPs coloured by sample metadata, including time points, mouse identities, and disease condition.

#### Cell Type Annotation and Visualization
Cell type annotation was performed using a manual, marker-based approach. Cluster-specific marker genes were identified using Seurat’s `FindAllMarkers` function with a log fold-change threshold of 0.25. For each cluster, the top five marker genes ranked by average log2 fold-change were selected to represent dominant expression signatures.
<br><br>
Clusters were classified into cell types according to the known cell type(s) of the dominant marker genes through a literature search with consideration of the biological context of nasal mucosa tissue. In cases where clusters exhibited mixed or ambiguous marker profiles, cell type identities were consolidated based on the predominant expression pattern across multiple markers.
<br><br>
Feature plots were generated throughout the annotation process to visualize the spatial distribution and specificity of marker gene expression across clusters, enabling validation of candidate cell type assignments (`FeaturePlot`). Representative feature plots were generated for both well-defined clusters and clusters with mixed marker signatures, as well as for key marker genes across major cell types.
Following annotation, cluster identities were relabeled using Seurat’s `RenameIdents` function, and annotated cell types were visualized using UMAP with cluster labels displayed.

#### Differential Expression Analysis
The Seurat object was inspected for NA or missing data and samples missing Mouse identifiers (IDs) were removed to ensure independent biological replicates and to exclude potential doublets or ambiguous cell assignments. 
<br><br>
To account for biological replication, gene expression counts were aggregated across cells for each mouse–timepoint combination (n = 10) and mouse–tissue combination (n = 3) using Seurat’s `AggregateExpression()` function. 
This pseudobulk approach enabled differential expression testing at the level of biological replicates, rather than individual cells. Differential expression analysis was performed using Seurat’s `FindMarkers` function with the `DESeq2` method to account for between-sample variability by modelling counts using a negative binomial distribution and reduce false positives (Love et al., 2014).
<br><br>
Volcano plots were generated for the time point comparisons relative to day 2 (peak viral infection), excluding Naïve (URT infection) (n=3), and for each tissue type comparison (n=3). 

#### Functional Enrichment Analysis
Over-representation analysis (ORA) was conducted using the `enrichGO` function from clusterProfiler (v.4.16.0) using the `org.Mm.eg.db` Mouse annotation database (v.21.0). Genes were filtered based on an adjusted p-value < 0.05 and an absolute log2 fold-change > 0.25, consistent with commonly used thresholds in single-cell RNA-seq Seurat workflows to capture biologically meaningful expression changes while filtering out technical noise. Gene Ontology (GO) Biological Process (BP) terms were evaluated using gene symbols as identifiers, with Benjamin-Hochberg multiple testing correction applied. A background universe consisting of a unique list of the tested genes was used to control for selection bias and results were visualized using a dot plot (`dotplot`).

### Results

#### Cell Clustering Reveals Distinct Populations
Quality control filtering was performed by assessing feature counts and mitochondrial gene expression. Comparable distribution of these metrics was observed across time points, and thresholds were applied accordingly to remove low-quality cells and potential outliers (n=149064; Figure 1).
x	x
Principal component analysis (PCA) was performed to reduce dimensionality, and the first 15 principal components were selected for downstream clustering based on the elbow plot (Figure 2). UMAP visualization revealed clear separation of cells into distinct clusters, indicating substantial cellular heterogeneity within the dataset (Figure 3). Clustering at a resolution of 0.5 produced well-defined and biologically interpretable groups within excessive fragmentation. UMAPs grouped by timepoints, disease condition and mouse identities showed strong mixing of cells by these variables (Figure 4). This suggests that clustering was driven by biological variation rather than technical artifacts and therefore, no batch effect correction was applied. 

#### Manual Annotation Identifies Major Cell Types
Clusters were annotated based on the top five expressed genes by average log2 fold change. Feature plots demonstrated that several clusters exhibited consistent and cell type-specific expression patterns. For example, cluster 4 showed strong and localized expression of B cell markers including Iglc2, Fcmr, Iglc1, Ighd, and Ms4a1 supporting its annotation as B lymphocytes (Figure 5). In contrast, some clusters displayed more heterogenous expression. Cluster 33 contained expression of genes associated with multiple cell types, including glandular cells (Sult1e1), immune cells (Tac4), chondrocytes (Col9a1), epithelial (Svopl) and possible associations with proliferating cells (Hist1h2ap), suggesting a mixed or less well-defined cell type (Figure 6). Feature plots of marker genes across multiple clusters further confirmed the separation of cell type populations, including basal epithelium (Krt15), macrophages (Ms4a7), neutrophils (Ly6g), B cells (Ms4a1), and endothelial (Ptprb) (Figure 7). Given the focus on immune response to viral infection, macrophages were selected for downstream analysis. Genes associated with macrophages in cluster 2, including Cd209f, Cd5l, and Ms4a7 showed clustering, confirming the cell type of this population (Figure 8). 

### Discussion

## References
