# scRNA-seq-nasal-mucosa
Single-cell RNA-seq analysis of mouse nasal mucosa to characterize cellular heterogeneity and host response to viral infection

**Author:** Rebekah Hest  
**Course:** BINF6110 - Genomic Methods for Bioinformatics

## Assignment 4

This repository contains Assignment 4, a scRNA-seq analysis of mouse nasal mucosa using Seurat, including clustering, manual cell type annotation, differential expression, and pathway enrichment to investigate cellular responses to viral infection

### Introduction

Single-cell RNA sequencing (scRNA-seq) has emerged as a powerful tool for studying gene expression in complex biological systems, offering key advantages over traditional bulk RNA sequencing (RNA-seq). While inexpensive, bulk RNA-seq methods average gene expression across many cells at the tissue level which can mask cell diversity and obscure cell-specific information (Mao et al., 2023). In contrast, scRNA-seq can isolate and capture transcriptional profiles at the individual cell-level, enabling characterization of heterogeneous cells and identification of distinct cell populations (Tzec-Interián, González-Padilla and Góngora-Castillo, 2025). This high-resolution technology can be used for the detection and classification of rare but functionally important cell types with applications in disease diagnosis and treatment (Jovic et al., 2022).
<br><br>
This approach is particularly valuable for investigating virus-host interactions and informing the development of therapeutics and vaccines (Chang et al., 2024). Respiratory viral infections, including SARS-CoV-2, typically enter the body through the mouth or nose, where the nasal epithelium serves as a critical immunological barrier (Chen and Wang, 2023; Ziegler et al., 2021). The nasal epithelium is a complex and heterogeneous tissue composed of multiple major and rare cells and the nasal mucosa is highly dynamic, exhibiting changes in response to environmental and infectious stimuli (Cybulski et al., 2026). Such heterogeneity and plasticity make it essential to study gene expression at single-cell resolution in order to fully capture the range of cellular responses to infection and understand how these responses contribute to disease progression and immune defense.
	However, the advantages of scRNA-seq come with increased analytical complexity. Compared to bulk RNA-seq, scRNA-seq data are high-dimensional, sparse, and subject to technical variability such as dropout events and batch effects (Kharchenko, 2021). As a result, analysis requires a multi-step pipeline including quality control, normalization, dimensionality reduction, clustering, and cell type annotation before meaningful biological comparisons can be made (Chen, Ning and Shi, 2019).
  <br><br>
Within the scRNA-seq framework, several methodological choices exist at each stage of analysis. To reduce technical noise of scRNA-seq data, normalization is an essential step that will benefit downstream analyses. Traditional log-normalization methods, such as those implemented in Seurat, scale gene expression values by sequencing depth and apply a logarithmic transformation, offering a simple and computationally efficient approach (Stuart et al., 2019). In contrast, more advanced methods such as SCTransform use a model-based framework based on regularized negative binomial regression to account for technical noise and variance across genes, often improving clustering and differential expression analyses (Hafemeister and Satija, 2019). Clustering is a key step in scRNA-seq analysis, as it enables the identification of distinct cell populations based on transcriptional similarity. Partition-based methods such as k-means require the number of clusters to be specified a priori, which can be limiting in complex and heterogeneous datasets and can be sensitive to outliers resulting in failures to detect rare cell types (Zhang et al., 2023). In contrast, graph-based approaches (e.g., Seurat’s FindClusters) use cell–cell similarity networks to identify clusters in a data-driven manner, making them better suited for capturing the structure of single-cell data. Finally, differential expression analysis in scRNA-seq can be performed either at the single-cell level or by aggregating counts across cells within clusters (pseudobulk). Single-cell methods, such as MAST, model expression at the level of individual cells but can be sensitive to technical noise and inflated false positive rates. In contrast, pseudobulk approaches using tools such as DESeq2 aggregate counts by sample and cell type, providing more robust statistical inference and improved control of false positives, making them a preferred approach for many scRNA-seq studies (Kalantari-Dehaghi, Ghohabi-Esfahani and Emadi-Baygi, 2025).
<br><br>
In this study, a scRNA-seq dataset from Kazer et al. (2025) of murine nasal mucosa will be analyzed to investigate the spatial and temporal varion in host response across distinct tissue types. To address these objectives, a Seurat-based analysis pipeline including clustering and UMAP visualization, cell type annotation, differential expression analysis, and pathway enrichment analysis will be applied. These methods were selected to balance biological interpretability with statistical rigor, enabling the identification of key cell populations and transcriptional pathways involved in the immune response to viral infection.

### Methods

#### Data Acquisition and Preprocessing

#### Normalization and Feature Selection

#### Dimensionality Reduction and Clustering

#### Batch Effect Correction

#### Cell Type Annotation and Visualization

#### Differential Expression Analysis

#### Functional Enrichment Analysis

### Results

### Discussion

## References
