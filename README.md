# transcriptomics-toolbox
Interactive spatial transcriptomics dashboard for Visium data with cluster/gene expression views, UMAP switching, gene finder and violin plot function.

# Spatial Transcriptomics Toolbox
<img width="2940" height="1677" alt="image" src="https://github.com/user-attachments/assets/6b9115fc-02af-46d8-860b-9fb971fa4363" />

## About

A self-contained HTML report with spatial and UMAP views, gene search and gene expression maps along with violin plots. Using deconvolve_reference.py HTML dashboard can apply resulting cell types to specific data. Deconvolution is based on a representative file available at the NCBI GEO (GSE127465). Notice: the project was developed to work with Visium data and used specifically Human Lung Cancer (FFPE), but may be adapted to other data as well. 

Built with **Scanpy**, **Squidpy**, **Plotly**, and **Pandas**.

---

## Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/kzhezhel/transcriptomics-toolbox.git
cd transcriptomics-toolbox
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
