#!/usr/bin/env python3

# --- Aim
# Make plots for a list of gene for a single sample


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
import math

import re

from matplotlib.colors import LinearSegmentedColormap
from mycolorpy import colorlist as mcp


# Import custom functions
from sophiesfunctions.auto_crop import auto_crop_scanpy  # import my custom functions that are packaged locally as sophiesfunctions
from sophiesfunctions.save_multi_image import save_multi_image

sc.set_figure_params(facecolor="white", figsize=(8, 8))
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


# --- Get some info from the sample sheet

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


# --- Print which data using
print(f'Use data: {use_data}')

# --- Calculate crop coords for optimal presentation
coords = auto_crop_scanpy(adata, border_size=0.1)


# --- Make spot size a bit larger for plotting
spot_size = 1.5

# --- Read in gene list
gene_file = 'marker_genes.csv'
gene_dir = os.path.join(project_dir, 'resources', 'gene_lists')

gene_table = pd.read_csv(os.path.join(gene_dir, gene_file))

if 'gene' in gene_table.columns:
    gene_list = list(gene_table.gene)
else:
    gene_list = list(gene_table.iloc[:, 0]) # Assuming the genes are in the first column of this table

# Remove any potential nan values
gene_list = [x for x in gene_list if isinstance(x, str) or not np.isnan(x)]

# Select only the genes that are present in the anndata
to_plot = np.intersect1d(gene_list , adata.var.index).tolist()

missing = list(set(gene_list) - set(to_plot))

print(f'Genes not found: {missing}') 

# --- Divide genes into groups based on cell_type
gene_dict = {}
group_key = 'cell_type'

cell_types = gene_table[group_key].unique()

for cell_type in cell_types:
    gene_dict[cell_type] = list(gene_table[gene_table[group_key] == cell_type].iloc[:,0])

# Remove any nans
for cell_type in gene_dict:
    gene_dict[cell_type] = [x for x in gene_dict[cell_type] if isinstance(x, str) or not np.isnan(x)]

# Also remove any missing genes
gene_dict = {key: [value for value in values if value not in missing] for key, values in gene_dict.items()}


# --- Plot the genes

# Plotting
ncols = 4

# --- Plot 1: Overview spatial plot for all samples, plotted by cell type

# Define a colormap where all values below a cutoff (vmin = 0.1) are 50% transparent
vmin = 0.1

cmap = "YlOrRd"
color_list = mcp.gen_color(cmap=cmap, n=11)

color_list_extracted = color_list[1:]

n_bins = 100
cmap_name = 'my_list'
cmap_adj = LinearSegmentedColormap.from_list(cmap_name, color_list_extracted, N=n_bins)

cmap_adj.set_under(color=color_list[0], alpha=0.7)


# Make the font a bit bigger
sc.set_figure_params(scanpy=True, fontsize=18)

for cell_type in gene_dict:
    nrows = math.ceil((len(gene_dict[cell_type]) + 1)/ncols)
    if ncols < (len(gene_dict[cell_type])+1):
        plt.figure(figsize=(7*ncols, 6*nrows))
        plt.suptitle(f'Marker genes: {cell_type} - use_data = {use_data} - sample_type = {sample_type} - vmin={vmin}', y=1, fontweight='bold')
    else:
        plt.figure(figsize=(7*ncols, 5.5*nrows))
        plt.subplots_adjust(wspace=0.2)
        plt.suptitle(f'Marker genes: {cell_type} - use_data = {use_data} - sample_type = {sample_type} - vmin={vmin}', y=1.1, fontweight='bold')
    for n, gene in enumerate(gene_dict[cell_type]):
        if n == 0:
            ax1 = plt.subplot(nrows, ncols, n * 2 + 1)
            ax2 = plt.subplot(nrows, ncols, n * 2 + 2)
            # Plot
            sc.pl.spatial(adata, color=None, size=spot_size, bw=False, alpha_img=1, cmap=cmap_adj, layer=use_layer,
                crop_coord=coords, title=f'{sample_id} - H&E', ax=ax1, vmin=vmin)
            # add two new subplots iteratively
            ax = plt.subplot(nrows, ncols, n + 1)
            # Plot
            sc.pl.spatial(adata, color=gene, size=spot_size, bw=False, alpha_img=0.5, cmap=cmap_adj, layer=use_layer,
                crop_coord=coords, title=f'{sample_id} - {gene}', ax=ax2, vmin=vmin)
        else:
            n+=1
            ax = plt.subplot(nrows, ncols, n + 1)
            # Plot
            sc.pl.spatial(adata, color=gene, size=spot_size, bw=False, alpha_img=0.5, cmap=cmap_adj, layer=use_layer,
                crop_coord=coords, title=f'{sample_id} - {gene}', ax=ax, vmin=vmin)
    plt.subplots_adjust(hspace=0.3, wspace=0.1)



# --- Save output
# Put in a timestamped folder to avoid overwriting older plots
qc_dir = os.path.join(res_dir, 'marker_genes')
os.makedirs(qc_dir, exist_ok=True)

os.chdir(qc_dir)

filename = f"{sample_id}_marker_genes_spatial.pdf"

today = date.today()
today = today.strftime("%Y%m%d")

out_dir = os.path.join(qc_dir, 'per_sample', today)
os.makedirs(out_dir, exist_ok=True)

print(f'Saving plot as {filename} in {out_dir}')
save_multi_image(os.path.join(out_dir, filename))

print("Script succesfully completed")

