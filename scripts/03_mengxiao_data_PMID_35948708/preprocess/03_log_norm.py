#!/usr/bin/env python3

# --- Aim
# # Use scanpy's sc.pp.normalize_total & sc.pp.log1p to normalize the data


# --- Define environment
# env: preprocessing

# --- Load packages
import argparse
import scanpy as sc
import pandas as pd
import anndata as ad

import numpy as np
import os


sc.logging.print_header()
sc.settings.verbosity = 3

# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & save as .h5ad file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--raw_data_dir", help="Path to raw data")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--min_count", help="Minimal number of counts per spot", type=int)
parser.add_argument("--min_gene", help="Minimal number of genes per spot", type=int)
parser.add_argument("--outs_path", help="Path to spaceranger 'outs'")
parser.add_argument("--use_data", help="Which data to use (log_norm or SCT)")


args = vars(parser.parse_args())


# --- Arguments
raw_data_dir = args["raw_data_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]

min_count=args["min_count"]
min_gene=args["min_gene"]


# --- Load the raw anndata object

input_adata_dir = os.path.join(project_dir, h5ad_dir, 'raw')
adata = ad.read(os.path.join(input_adata_dir, f'{sample_id}_raw.h5ad'))


# --- Do some very conservative filtering: only spots with 0 counts

sc.pp.filter_cells(adata, min_counts=1)

adata.layers['filtered_counts'] = adata.X


# --- Normalization

adata.layers['norm'] = sc.pp.normalize_total(adata, layer='filtered_counts', target_sum=1e4, inplace=False)['X']

# log transform
adata.layers['log_norm'] = sc.pp.log1p(adata.layers['norm'], copy=True)

# Scale gene expression to unit variance (clip values that exceed 10 standard deviations)
adata.layers['scaled'] = sc.pp.scale(adata.layers['log_norm'], max_value=10, copy=True)

# Detect highly variable genes as well
sc.pp.highly_variable_genes(adata, layer='log_norm', flavor='seurat', n_top_genes=2000)


# --- PCA, UMAP + clustering

seed = 42

# create a new object
pca_obj = ad.AnnData(adata.layers['scaled'])
pca_obj.var = adata.var
sc.tl.pca(pca_obj, n_comps = 50, zero_center = True, svd_solver='arpack', random_state = seed, return_info = True, use_highly_variable = True)

adata.uns['pca'] = pca_obj.uns['pca']
adata.obsm['X_pca'] = pca_obj.obsm['X_pca']
adata.varm['PCs'] = pca_obj.varm['PCs']

# neighbourhood graph
sc.pp.neighbors(adata, n_pcs=50, random_state=seed)

# calculate umap
sc.tl.umap(adata, random_state=seed)

# Leiden clustering
sc.tl.leiden(adata, key_added='leiden', resolution=1)


# --- Save output

outs_dir = os.path.join(h5ad_dir, 'log_norm')

os.makedirs(outs_dir, exist_ok=True)

print(f"Saving .h5ad files to {outs_dir}")

ad.AnnData.write(adata, filename=os.path.join(outs_dir, f"{sample_id}_log_norm.h5ad"))

print("Script succesfully completed")
