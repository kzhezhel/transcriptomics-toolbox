import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.io import mmread
import scanpy as sc
import os
import warnings
import argparse
import plotly.graph_objects as go

warnings.filterwarnings('ignore', category=RuntimeWarning)

parser = argparse.ArgumentParser(
    description="Deconvolve Visium spatial transcriptomics data using NNLS with a reference scRNA-seq dataset."
)
parser.add_argument(
    "--visium",
    type=str,
    default="lung_cancer_processed.h5ad",
    help="Path to the processed Visium AnnData file (default: lung_cancer_processed.h5ad)."
)
parser.add_argument(
    "--ref_dir",
    type=str,
    default="spatial/reference",
    help="Path to the directory containing reference files (default: spatial/reference)."
)
parser.add_argument(
    "--output",
    type=str,
    default="lung_cancer_deconvolved_nnls.h5ad",
    help="Output file name for deconvolved AnnData (default: lung_cancer_deconvolved_nnls.h5ad)."
)
parser.add_argument(
    "--no_plot",
    action="store_true",
    help="Skip generating the spatial plot."
)
args = parser.parse_args()

# Existing files test

if not os.path.exists(args.visium):
    raise FileNotFoundError(f"Visium file not found: {args.visium}")

ref_files = [
    "GSE127465_human_counts_normalized_54773x41861.mtx.gz",
    "GSE127465_gene_names_human_41861.tsv.gz",
    "GSE127465_human_cell_metadata_54773x25.tsv.gz"
]
for f in ref_files:
    if not os.path.exists(os.path.join(args.ref_dir, f)):
        raise FileNotFoundError(f"Reference file not found: {os.path.join(args.ref_dir, f)}")

print("=" * 60)
print("NNLS DECONVOLUTION PIPELINE")
print("=" * 60)


# Data uploading

print("\n[1/6] Loading data...")
adata_vis = sc.read(args.visium)
print(f"Visium: {adata_vis.n_obs} spots, {adata_vis.n_vars} genes")


# Uploading reference data

print("\n[2/6] Loading reference data...")
mat = mmread(os.path.join(args.ref_dir, "GSE127465_human_counts_normalized_54773x41861.mtx.gz")).tocsr()
genes = pd.read_csv(
    os.path.join(args.ref_dir, "GSE127465_gene_names_human_41861.tsv.gz"),
    header=None, compression='gzip'
)[0].values
meta = pd.read_csv(
    os.path.join(args.ref_dir, "GSE127465_human_cell_metadata_54773x25.tsv.gz"),
    sep='\t', compression='gzip'
)

if mat.shape[0] == len(genes):
    mat = mat.T

adata_ref = sc.AnnData(X=mat)
adata_ref.var_names = genes
adata_ref.obs = meta
print(f"Reference: {adata_ref.n_obs} cells, {adata_ref.n_vars} genes")

# 3. Filtering reference

print("\n[3/6] Preprocessing reference...")
sc.pp.log1p(adata_ref)

if 'Major cell type' not in adata_ref.obs.columns:
    raise ValueError("Column 'Major cell type' not found in reference metadata")

def clean_cell_type(major):
    if 'Patient' in major and 'specific' in major:
        return None
    return major

adata_ref.obs['clean_cell_type'] = adata_ref.obs['Major cell type'].apply(clean_cell_type)
initial_count = adata_ref.n_obs
adata_ref = adata_ref[adata_ref.obs['clean_cell_type'].notna()].copy()
print(f"Removed {initial_count - adata_ref.n_obs} cells with Patient-specific types")
print(f"Remaining: {adata_ref.n_obs} cells")

# Common genes
print("\n[4/6] Finding common genes...")
common_genes = np.intersect1d(adata_vis.var_names, adata_ref.var_names)
print(f"Common genes: {len(common_genes)}")
adata_vis = adata_vis[:, common_genes].copy()
adata_ref = adata_ref[:, common_genes].copy()

# Calculating cell types

print("\n[5/6] Computing cell type signatures...")
cell_types = adata_ref.obs['clean_cell_type'].unique()
cell_types = [ct for ct in cell_types if ct != "ND"]
print(f"Cell types (excluding ND): {cell_types}")

signatures = {}
for ct in cell_types:
    cells = adata_ref.obs_names[adata_ref.obs['clean_cell_type'] == ct]
    if len(cells) == 0:
        continue
    expr = adata_ref[cells, :].X.mean(axis=0)
    signatures[ct] = np.asarray(expr).flatten()

S = np.column_stack([signatures[ct] for ct in signatures.keys()])
type_names = list(signatures.keys())
print(f"Found {len(type_names)} cell types: {type_names}")

# Deconvolution for each spot

print("\n[6/6] Deconvolving spots...")
n_spots = adata_vis.n_obs
proportions = np.zeros((n_spots, len(type_names)))
X_spots = adata_vis.X

for i in range(n_spots):
    x = X_spots[i].toarray().flatten() if hasattr(X_spots, 'toarray') else X_spots[i]
    coef, _ = nnls(S, x)
    s = coef.sum()
    if s > 0:
        coef = coef / s
    proportions[i, :] = coef

pred_df = pd.DataFrame(proportions, index=adata_vis.obs_names, columns=type_names)

for ct in type_names:
    adata_vis.obs[f'pred_{ct}'] = pred_df[ct].values
adata_vis.obs['pred_cell_type'] = pred_df.idxmax(axis=1)


# Save results

print(f"\nSaving results to {args.output}...")
adata_vis.write(args.output)
print(f"✅ Deconvolution complete. Results saved to {args.output}")

# Statistics and visualisation

print("\n" + "=" * 60)
print("STATISTICS")
print("=" * 60)

print("\nPredicted cell type distribution:")
print(adata_vis.obs['pred_cell_type'].value_counts())

if not args.no_plot:
    print("\nGenerating spatial plot...")
    coords = adata_vis.obsm['spatial']
    fig = go.Figure()
    fig.add_trace(go.Scattergl(
        x=coords[:, 0],
        y=coords[:, 1],
        mode='markers',
        marker=dict(
            size=4,
            color=adata_vis.obs['pred_cell_type'].astype('category').cat.codes,
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title="Cell type")
        ),
        text=adata_vis.obs['pred_cell_type'],
        hoverinfo='text'
    ))
    fig.update_layout(
        title="Predicted Cell Types (NNLS Deconvolution)",
        width=800,
        height=600
    )
    fig.write_html("nnls_deconvolution_full.html")
    print("✅ Spatial plot saved to nnls_deconvolution_full.html")

print("\n✅ Done!")
