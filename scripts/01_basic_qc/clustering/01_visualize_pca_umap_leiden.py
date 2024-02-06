#!/usr/bin/env python3

# --- Aim
# Visualize the results of PCA & UMAP


# --- Define environment
# env: pca-visium

# --- Load packages
import argparse
import scanpy as sc
import anndata as ad
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import date

import matplotlib.patches as mpatches

# Import custom functions
from sophiesfunctions.auto_crop import auto_crop_scanpy  # import my custom functions that are packaged locally as sophiesfunctions
from sophiesfunctions.save_multi_image import save_multi_image
from sophiesfunctions.color_functions import generate_color_variations
from sophiesfunctions.misc_functions import flatten
from sophiesfunctions.plotting_aids import adjust_label_positions

sc.set_figure_params(facecolor="white", figsize=(6, 6))
sc.settings.verbosity = 3


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & do basic QC",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--data_dir", help="Path to raw data")
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--unique_id", help="Unique sample ID")
parser.add_argument("--seed", help="Seed", type=int)

parser.add_argument("--row_max", help="Row max for cropping Visium data", type=int)
parser.add_argument("--col_max", help="Col max for cropping Visium data", type=int)

parser.add_argument("--min_count", help="Min count for QC", type=int)
parser.add_argument("--min_gene", help="Min gene for QC", type=int)
parser.add_argument("--min_spots", help="Min number of spots each gene is expressed in", type=int)
parser.add_argument("--use_data", help="Which data to use (SCT or log_norm)") # Not necessary for this script

parser.add_argument("--gene_file", help="File with gene names to plot") # Not necessary for this script

args = vars(parser.parse_args())


# --- Arguments
data_dir = args["data_dir"]
res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]

project_dir=args["project_dir"]

use_data=args["use_data"]


# --- Print which data using

print(f'Use data: {use_data}')

# Extract the section_id
sample_dir = os.path.join(project_dir, 'config')
sample_sheet_file = 'sample_sheet.csv'

sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))
sample_lut = sample_sheet.set_index('sample_id')

sample_name = sample_lut.loc[sample_id, 'sample_name']
sample_type = sample_lut.loc[sample_id, 'type']


# -- Load data

if use_data == "SCT":
    adata = ad.read(os.path.join(h5ad_dir, 'SCTransform', f'{sample_name}_SCT.h5ad'))
else:
    adata = ad.read(os.path.join(h5ad_dir, 'log_norm', f'{sample_name}_log_norm.h5ad'))


# --- Specify layer to use

if use_data == "SCT":
    use_layer = "SCT_data"
else:
    use_layer = "log_norm"


# --- Calculate crop coords for optimal presentation
coords = auto_crop_scanpy(adata, border_size=0.1)

# --- Make spot size a bit larger for plotting
spot_size = 1.5

# --- Check PCA results

sc.set_figure_params(scanpy=True, fontsize=15)

# --- Plot 1: Variance ratio
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
plt.suptitle(f'{sample_id} ({sample_type}) - PCA variance', y=1.1, fontweight='bold')

plt.savefig('test.pdf', bbox_inches='tight')

# Plot 2: Plot the genes driving each PC
# Note: this doesn't work with the SCT data because you can't save the PCs
if use_data == 'log_norm':
    dimensions_to_check=10
    plt.rcParams['figure.figsize']=(5, 5)
    sc.pl.pca_loadings(adata, components = range(1, dimensions_to_check+1, 1), n_points = 10)
    plt.suptitle(f'{sample_id} ({sample_type}) PCA loadings (top 10) - Use data: {use_data}', y=1, fontweight='bold')


# --- Check UMAP results

# --- Plot 3: Plot qc metrics in UMAP

plt.rcParams["figure.figsize"] = (5, 5)

# Calculate log(total_counts for plotting)
adata.obs['log_counts'] = np.log(adata.obs['total_counts'])

to_plot = ['total_counts', 'log_counts', 'n_genes_by_counts', 'pct_counts_mt']

sc.pl.umap(adata, wspace=0.3, color = to_plot, ncols=len(to_plot))
plt.suptitle(f'{sample_id} ({sample_type}) UMAP qc metrics - Use data: {use_data}', y=1.1, fontweight='bold')


# --- Plot 4: Leiden clusters on UMAP

sc.pl.umap(adata, color = 'leiden')
plt.suptitle(f'{sample_id} ({sample_type}) - Leiden clusters - Use data: {use_data}', y=1.1, fontweight='bold')


# --- Plot 5: Leiden clusters spatially

ncols = 3

fig, axs = plt.subplots(ncols=ncols, figsize=(5*ncols, 5), sharey=True)

n=0
sc.pl.spatial(adata, color=None, size=spot_size, bw=False, alpha_img=1, crop_coord=coords, title=f'{sample_id} - H&E', ax=axs[n])

n+=1
sc.pl.spatial(adata, color='leiden', bw=False, alpha_img=0.5, crop_coord=coords, title=f'{sample_id} - Leiden', ax=axs[n])

n+=1
sc.pl.spatial(adata, color='leiden', bw=False, size=spot_size, alpha_img=0.5, crop_coord=coords, title=f'{sample_id} - Leiden (big spots)', ax=axs[n])

plt.tight_layout()
plt.suptitle(f'{sample_id} ({sample_type}) - Leiden clusters spatially - Use data: {use_data}', y=1.1, fontweight='bold')

plt.savefig('test.pdf', bbox_inches='tight')


# --- Extract marker genes

rank_genes_key = 'rank_genes_leiden'
sc.tl.rank_genes_groups(adata, 'leiden', use_raw = False, layer=use_layer, method='wilcoxon', key_added=rank_genes_key)


# --- Plot 6: Heatmaps of the top 10 genes per cluster

sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, groupby='leiden', key=rank_genes_key, layer=use_layer)
plt.suptitle(f'{sample_id} ({sample_type}) - Top 10 genes per cluster - Use data: {use_data}', y=1.1, fontweight='bold')


# --- Plot 7: Dot plots for top 10 genes per cluster


sc.pl.rank_genes_groups_dotplot(adata, n_genes=10, groupby='leiden', key=rank_genes_key, layer=use_layer)
plt.suptitle(f'{sample_id} ({sample_type}) - Top 10 genes per cluster - Use data: {use_data}', y=1.1, fontweight='bold')


# --- Plot 8: 1 vs rest plots

sc.pl.rank_genes_groups(adata, n_genes=25, groupby='leiden', key=rank_genes_key, layer=use_layer)
plt.suptitle(f'{sample_id} ({sample_type}) - Leiden clusters 1 vs rest plots, 25 genes - Use data: {use_data}', y=1.1, fontweight='bold')


# --- Save output
# Put in a timestamped folder to avoid overwriting older plots
qc_dir = os.path.join(res_dir, 'clustering', 'initial_qc_clustering')
os.makedirs(qc_dir, exist_ok=True)

os.chdir(qc_dir)

filename = f"{sample_id}_PCA_UMAP_leiden.pdf"

today = date.today()
today = today.strftime("%Y%m%d")

out_dir = os.path.join(qc_dir, today)
os.makedirs(out_dir, exist_ok=True)

print(f'Saving plot as {filename} in {qc_dir}')
save_multi_image(os.path.join(out_dir, filename))

print("Script succesfully completed")