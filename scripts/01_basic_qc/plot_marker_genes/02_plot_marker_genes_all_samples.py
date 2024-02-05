#!/usr/bin/env python3

# --- Aim
# Plot marker genes for all samples as a heatmap

# --- Define environment
# env: pca-visium

# --- Load packages
import scanpy as sc
import anndata as ad
import pandas as pd

import numpy as np
import os

from datetime import date

import matplotlib.pyplot as plt

import argparse

# sc.logging.print_versions()
# sc.set_figure_params(facecolor="white", figsize=(8, 8))
# sc.settings.verbosity = 3

# Custom functions
from sophiesfunctions.save_multi_image import save_multi_image


# --- Arguments

# --- Argparse arguments
parser = argparse.ArgumentParser(description="Import Visium unfiltered data & plot",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--data_dir", help="Path to raw data")
parser.add_argument("-s", "--samples",
    nargs="*",  # 0 or more values expected => creates a list
    type=str,
    help="List of samples")
parser.add_argument("-r", "--res_dir",
                    help="Path to results directory (optional)")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files") # Not necessary for this script
parser.add_argument("--gene_file", help="Path to a csv file with genes, genes stored in first column") # Not necessary for this script
parser.add_argument("--use_data", help="Which data to use (SCT or log_norm)") # Not necessary for this script
parser.add_argument("--project_dir", help="Project directory (root)")

args = vars(parser.parse_args())

# Arguments
data_dir = args["data_dir"]
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]
sample_list = args["samples"]
project_dir = args["project_dir"]
use_data = args["use_data"]

print(f'Samples in python: {sample_list}')

# Make directories if they don't exist yet
if not os.path.exists(res_dir):
    os.makedirs(res_dir)

# --- Specify layer to use

if use_data == "SCT":
    use_layer = "SCT_data"
else:
    use_layer = "log_norm"

# --- Print which data using
print(f'Use data: {use_data}')


# --- Get some info from the sample sheet

# Extract the section_id
sample_dir = os.path.join(project_dir, 'config')
sample_sheet_file = 'sample_sheet.csv'

sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))
sample_lut = sample_sheet.set_index('sample_id')


# -- Load data

subtype_key = 'sample_type'
sample_id_key = 'sample_id'

# Multiple samples as a dictionary seems like an okay idea
adata_dict = {}

if use_data == 'SCT':
    for sample in sample_list:
        sample_name = sample_lut.loc[sample, 'sample_name']
        adata_dict[sample] = ad.read(os.path.join(h5ad_dir, 'SCTransform', f'{sample_name}_SCT.h5ad'))
        # Add type to obs for easier splitting
        adata_dict[sample].obs[subtype_key] = sample_lut.loc[sample, 'type']
        # Also add sample id for easier splitting
        adata_dict[sample].obs[sample_id_key] = sample
else:
    for sample in sample_list:
        sample_name = sample_lut.loc[sample, 'sample_name']
        adata_dict[sample] = ad.read(os.path.join(h5ad_dir, 'log_norm', f'{sample_name}_log_norm.h5ad'))
                # Add type to obs for easier splitting
        adata_dict[sample].obs[subtype_key] = sample_lut.loc[sample, 'type']
        # Also add sample id for easier splitting
        adata_dict[sample].obs[sample_id_key] = sample

# --- Order based on unique sample ID

ordered_keys = adata_dict.keys()

ordered_keys = sorted(ordered_keys)


# --- Convert to adata list & concatenate (easier for plotting)

adata_list = list(adata_dict.values())

all_adata = ad.concat(adata_list, merge='same')

# --- Prepare for plotting

# Split the adatas by subtype for plotting

all_subtypes = sample_sheet['type'].unique()
adata_per_subtype = {}

del adata_dict

for subtype in all_subtypes:
    adata_per_subtype[subtype] = all_adata[all_adata.obs[subtype_key] == subtype]

ncols = len(all_subtypes)

# Calculate the width of each subtype subplot based on how many samples per subtype

grid_spec = []

for subtype in all_subtypes:
    grid_spec.append(len(adata_per_subtype[subtype].obs[sample_id_key].unique()))


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

# Select only the genes that are present in the anndata objects
to_plot = np.intersect1d(gene_list , all_adata.var.index).tolist()

missing = list(set(gene_list) - set(to_plot))

print(f'Genes not found: {missing}') 

# Divide genes into groups based on cell_type

gene_dict = {}

group_key = 'cell_type'

cell_types = gene_table[group_key].unique()

for cell_type in cell_types:
    gene_dict[cell_type] = list(gene_table[gene_table[group_key] == cell_type].iloc[:,0])

# Remove any nans
for cell_type in gene_dict:
    gene_dict[cell_type] = [x for x in gene_dict[cell_type] if isinstance(x, str) or not np.isnan(x)]

# --- Remove genes that are missing in all samples

gene_dict = {key: [value for value in values if value not in missing] for key, values in gene_dict.items()}


# --- Plotting

# Plot 1: Heatmap
sc.pl.matrixplot(all_adata, var_names=gene_dict, groupby=[subtype_key, sample_id_key], layer=use_layer, vmax=5)
plt.suptitle(f'Marker genes - use_data = {use_data}, vmax = 5', y=1.15, fontweight='bold')

# Plot 2: Dot plot
sc.pl.dotplot(all_adata, var_names=gene_dict, groupby=[subtype_key, sample_id_key], layer=use_layer, vmax=5)
plt.suptitle(f'Marker genes - use_data = {use_data}, vmax = 5', y=1.15, fontweight='bold')


# --- Same plots above but now split per subtype:

ncols=len(all_subtypes)
nrows=1

# Adjust the gridspec ratios to make room for the legend on the last plot
grid_spec_adj = grid_spec.copy()
grid_spec_adj[-1] = grid_spec_adj[-1]*1.2
grid_spec_adj[-1] = grid_spec[-1]+5

# Calculate the height based on the number of genes to plot
height = 5/15* gene_table.shape[0] - len(missing)

# -- Plot 3: Heatmap split by subtype

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(3*ncols, height*nrows), gridspec_kw={'width_ratios': grid_spec_adj})

for n, subtype in enumerate(all_subtypes):
    plot = sc.pl.MatrixPlot(adata_per_subtype[subtype], var_names=gene_dict, groupby=[sample_id_key], layer=use_layer, ax=axs[n], vmin=0, vmax=3)
    plot.swap_axes()
    if n == (ncols-1):
        plot.legend(width=1.0)
    else:
        plot.legend(False)
    plot_axes = plot.get_axes()
    # move x-axis labels to the left
    xlabels = plot_axes['mainplot_ax'].get_xticklabels()
    plot_axes['mainplot_ax'].set_xticklabels(xlabels, rotation=45, ha='right', fontsize=9)
    if n != 0:
        plot_axes['mainplot_ax'].set_yticklabels([])
        plot_axes['mainplot_ax'].set_yticks([])
    plot_axes['mainplot_ax'].set_title(subtype, fontsize=9, fontweight='bold')

plt.subplots_adjust(wspace=0)
plt.suptitle(f'Marker genes - use_data = {use_data} - all spots', y=0.93, fontweight='bold')


# -- Plot 4: Dotplot

# Set some default parameters for all plots
sc.pl.DotPlot.DEFAULT_COLORMAP='Reds'
sc.pl.DotPlot.DEFAULT_DOT_MIN=0
sc.pl.DotPlot.DEFAULT_DOT_MAX=1

# Plotting

fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(3*ncols, height*nrows), gridspec_kw={'width_ratios': grid_spec_adj})

for n, subtype in enumerate(all_subtypes):
    plot = sc.pl.DotPlot(adata_per_subtype[subtype], var_names=gene_dict, groupby=[sample_id_key], layer=use_layer, ax=axs[n], vmin=0, vmax=3.5)
    # Genes on the y axis
    plot.swap_axes()
    # Remove legend or make it smaller (on the last plot)
    if n == (ncols-1):
        plot.legend(width=0.8, size_title='Fraction of spots\nin group (%)')
    else:
        plot.legend(False)
    # Extract the axes so we can manipulate these
    plot_axes = plot.get_axes()
    plot_axes['mainplot_ax'].margins(x=0)
    # Move x-axis labels to the left
    xlabels = plot_axes['mainplot_ax'].get_xticklabels()
    plot_axes['mainplot_ax'].set_xticklabels(xlabels, rotation=45, ha='right', fontsize=9)
    # Remove y-axis labels from all plots except the first one
    if n != 0:
        plot_axes['mainplot_ax'].set_yticklabels([])
        plot_axes['mainplot_ax'].set_yticks([])
    # Add titles, adjust the positioning a bit for Luminal B (HER2+) so it doesn't overlap
    plot_axes['mainplot_ax'].set_title(subtype, fontsize=9, fontweight='bold')

plt.subplots_adjust(wspace=0)
plt.suptitle(f'Marker genes - use_data = {use_data} - all spots', y=0.93, fontweight='bold')


# --- Save output
# Put in a timestamped folder to avoid overwriting older plots
qc_dir = os.path.join(res_dir, 'marker_genes')
os.makedirs(qc_dir, exist_ok=True)

os.chdir(qc_dir)

filename = "Marker_genes_all_samples_heatmap.pdf"

today = date.today()
today = today.strftime("%Y%m%d")

out_dir = os.path.join(qc_dir, 'all_samples', today)
os.makedirs(out_dir, exist_ok=True)

print(f'Saving plot as {filename} in {qc_dir}')
save_multi_image(os.path.join(out_dir, filename))

print("Script succesfully completed")