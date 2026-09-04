# Spatial Transcriptomics Toolbox

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<img width="2940" height="1670" alt="screenshot" src="https://github.com/user-attachments/assets/5f286714-b30b-4568-b130-9da77aab7a60" />

**Interactive spatial transcriptomics dashboard for 10x Visium data.**  
Provides cluster/gene expression views, UMAP switching, gene search, and violin plots. Includes optional cell‑type deconvolution using NNLS with a reference scRNA‑seq dataset (GSE127465). Built with **Scanpy**, **Squidpy**, **Plotly**, and **Pandas**.

---

## 📌 About

This toolbox offers two main components:

1. **Interactive Dashboard** – a self‑contained HTML report with spatial and UMAP projections, gene expression maps, violin plots, and an integrated gene search.
2. **Cell‑type Deconvolution (NNLS)** – estimates cell type proportions per spot using non‑negative least squares, based on a reference single‑cell RNA‑seq dataset from GSE127465.

The project was developed for **Human Lung Cancer (FFPE)** Visium data but can be adapted to other spatial transcriptomics datasets.

---

## ✨ Features

### Dashboard
- **Spatial and UMAP projections** – switch between tissue image and UMAP embedding.
- **Cluster visualisation** – colour‑coded Leiden clusters with interactive hover.
- **Gene expression maps** – display log2‑CPM for any gene with a Viridis color scale.
- **Violin plots** – show gene expression distribution across clusters.
- **Gene search** – type-ahead dropdown for quick gene selection.
- **Opacity control** – adjust point transparency.
- **Three viewing modes**:
  - **Clusters** – view Leiden clusters on spatial/UMAP coordinates.
  - **Feature** – visualise gene expression for any gene.
  - **Patch** – view the raw tissue image without overlaid points.

### Deconvolution
- **Reference‑based cell‑type deconvolution** using NNLS (non‑negative least squares).
- Cleans reference data by removing patient‑specific cell types.
- Outputs predicted cell type proportions and dominant cell type per spot.
- Results are automatically integrated into the dashboard.

---

## 🧬 Methodological Approach

The toolbox follows a standard spatial transcriptomics analysis pipeline:

1. **Data loading** – reads 10x Visium `.h5` files and spatial metadata (images, coordinates).
2. **Preprocessing** – filters cells and genes, normalises (log‑CPM), scales, and selects highly variable genes.
3. **Dimensionality reduction** – PCA, UMAP, and Leiden clustering.
4. **Optional deconvolution** – uses NNLS with a reference scRNA‑seq dataset to predict cell type proportions.
5. **Interactive visualisation** – generates a self‑contained HTML dashboard with spatial and UMAP views, gene expression, and violin plots.

---

## 🚀 Installation

### 1. Clone the repository

```bash
git clone https://github.com/kzhezhel/transcriptomics-toolbox.git
cd transcriptomics-toolbox
