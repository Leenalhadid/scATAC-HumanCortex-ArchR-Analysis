# 🧬 Single-Cell ATAC-seq Analysis of Human Brain Cortex Cells

[![R](https://img.shields.io/badge/R--based-ArchR-blueviolet)](https://www.archrproject.com/)  
 
---

## 🧾 Project Overview

This repository presents the full pipeline and results of a single-cell ATAC-seq analysis of human brain cortex cells. The data were originally published by [Trevino et al. (2021)](https://www.cell.com/cell/fulltext/S0092-8674(21)00942-9) and analyzed using the ArchR framework.

### 🔍 Objectives
- Preprocessing and quality control
- Dimensionality reduction and batch correction
- Clustering and peak calling
- Identification of marker peaks and genes
- TF motif activity analysis
- Integration with scRNA-seq data
- Peak-to-gene linkage analysis
- Differential accessibility and motif enrichment

---

## 📦 Dataset Description

The dataset includes three replicates of human brain cortex cells:
- `w21_dc1r3_r1`
- `w21_dc2r2_r1`
- `w21_dc2r2_r2`

**Download**: [scbi_p2.zip](https://icbb-share.s3.eu-central-1.amazonaws.com/single-cell-bioinformatics/scbi_p2.zip)

---



### 🧪 Environment (R setup)

Install required packages in R:

\`\`\`r
install.packages("BiocManager")
BiocManager::install("ArchR")
BiocManager::install("chromVARmotifs")
\`\`\`

Optional dependencies for plots:

\`\`\`r
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("patchwork")
\`\`\`


