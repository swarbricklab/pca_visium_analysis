#!/usr/bin/env python3

# --- Aim
# Plot the raw feature matrix to see how much was captured outside of the tissue region
# For multiple samples at the same time (take adata dictionary object)
# Note: these plots are also very useful to determine which edges of the capture area are unreliable (have comaparatively high/low n/o genes all along the edge)
# --- Define environment
# env: pca-visium

# --- Load packages
import argparse
import scanpy as sc
import anndata as ad
from datetime import date
import numpy as np
import pandas as pd

import os

import matplotlib.pyplot as plt

from mycolorpy import colorlist as mcp
from matplotlib.colors import LinearSegmentedColormap


# Import custom functions
from sophiesfunctions.auto_crop import auto_crop_scanpy  # import my custom functions that are packaged locally as sophiesfunctions
from sophiesfunctions.save_multi_image import save_multi_image

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

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
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]
sample_list = args["samples"]
project_dir = args['project_dir']

print(f'Samples in python: {sample_list}')

# --- Get some info from the sample sheet

# Extract the section_id
sample_dir = os.path.join(project_dir, 'config')
sample_sheet_file = '20240716_sample_sheet.csv'

sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))
sample_lut = sample_sheet.set_index('section_name')


# -- Load data
# Multiple samples as a dictionary seems like an okay idea
adata_dict = {}

for sample in sample_list:
    sample_name = sample_lut.loc[sample, 'sample_name']
    adata_dict[sample] = ad.read(os.path.join(h5ad_dir, 'raw_feature_bc_matrix', f'{sample_name}_raw_feature_bc_matrix.h5ad'))


# --- Order based on unique sample ID

ordered_keys = adata_dict.keys()
ordered_keys = sorted(ordered_keys)

# --- Calculate crop coords for optimal presentation
coords = {}
for sample in adata_dict.keys():
    coords[sample] = auto_crop_scanpy(adata_dict[sample], border_size=0.1)

# Define a custom color map
color_list = mcp.gen_color(cmap="OrRd", n=10)
color_list_added = ['lightgray'] + color_list + ['black']

n_bins = 100
cmap_name = 'my_list'

cmap_test = LinearSegmentedColormap.from_list(cmap_name, color_list_added, N=n_bins)


# --- Calculate log total counts for each sample
for sample in adata_dict.keys():
    total_counts_copy = adata_dict[sample].obs['total_counts'] 
    total_counts_copy[total_counts_copy == 0] = 1e-15 # Deal with the 0s
    adata_dict[sample].obs['log_counts'] = np.log(total_counts_copy)

# --- Plot 1: for each sample, plot the total_counts per spot, log_counts, genes per spot, H&E and cytassist image

for sample in ordered_keys:
    # Check if there is a cytassist image
    cytassist = False
    sample_name = sample_lut.loc[sample, 'sample_name']
    if 'cytassist' in adata_dict[sample].uns['spatial'][sample_name]['images']:
        cytassist = True
        ncol = 6
        width_ratios = [1,1,1,1,1,0.9]
    else:
        ncol = 5
        width_ratios = [1,1,1,1,1]
    fig, axs = plt.subplots(1, ncol, figsize=(7*ncol, 5), gridspec_kw={'wspace': 0, 'width_ratios': width_ratios})
    sc.pl.spatial(adata_dict[sample], color=['in_tissue'], size=1.5, color_map=cmap_test, wspace=0, crop_coord=coords[sample], ax=axs[0])
    sc.pl.spatial(adata_dict[sample], color=['total_counts'], size=1.5, color_map=cmap_test, wspace=0, crop_coord=coords[sample], ax=axs[1])
    sc.pl.spatial(adata_dict[sample], color=['log_counts'], size=1.5, color_map=cmap_test, wspace=0, crop_coord=coords[sample], vmin=0, ax=axs[2])
    sc.pl.spatial(adata_dict[sample], color=['n_genes_by_counts'], size=1.5, color_map=cmap_test, wspace=0, crop_coord=coords[sample], ax=axs[3])
    sc.pl.spatial(adata_dict[sample], size=1.5, wspace=0, crop_coord=coords[sample], ax=axs[4], title=f"{sample} - H&E")
    if cytassist:
        axs[5].imshow(adata_dict[sample].uns['spatial'][sample_name]['images']['cytassist'])
        plt.grid(False)
        plt.axis('off')
        plt.title('cytassist')
    sample_type = sample_lut.loc[sample, 'type']
    plt.suptitle(f'{sample} ({sample_type}) unfiltered plots', fontsize=18, fontweight='bold', y=1.05, x=0.5, ha='center')


# --- Save output
# Put in a timestamped folder to avoid overwriting older plots
qc_dir = os.path.join(res_dir, 'basic_qc')
os.makedirs(qc_dir, exist_ok=True)

filename = f"All_samples_unfiltered_plots.pdf"

today = date.today()
today = today.strftime("%Y%m%d")

out_dir = os.path.join(qc_dir, 'all_samples', today)
os.makedirs(out_dir, exist_ok=True)

print(f'Saving plot as {filename} in {qc_dir}')
save_multi_image(os.path.join(out_dir, filename))

print("Script succesfully completed")
