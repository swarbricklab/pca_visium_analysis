#!/usr/bin/env python3

# --- Aim
# Use cell2location for deconvolution
# Plot cell type mapping results for quick check (does it look as expected)


# --- Define environment
# env: cell2location

# --- Load packages
import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl

from datetime import datetime

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

import pandas as pd
import os

# --- Argparse arguments                
parser = argparse.ArgumentParser(description="Import variables & train Cell2Location model",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--projectName", help="Main project name")
parser.add_argument("--repo", help="Repo name")
parser.add_argument("--exp", help="Exp name")
parser.add_argument("--analysis", help="Analysis step name")
parser.add_argument("--section_name", help="Section name")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files - Visium raw data")
parser.add_argument("--model_dir", help="Path to directory for output files from model training")

args = vars(parser.parse_args())

# Create directory structure using argparse arguments
projectName = args["projectName"]
repo = args["repo"]
exp = args["exp"]
analysis = args["analysis"]
section_name = args["section_name"]
h5ad_dir = args["h5ad_dir"]   # /share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/data/anndata_objects/
model_dir = args["model_dir"]   # /share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/data/anndata_objects/

# Define directory paths
repoDir = os.path.join('/share/ScratchGeneral/evaapo/projects/', projectName, repo)
resultDir = os.path.join(repoDir, "results", exp, analysis, section_name)
objectDir = os.path.join(resultDir, "objects")
figDir = os.path.join(resultDir, "figures")
tabDir = os.path.join(resultDir, "tables")
                        
# Ensure output directories exist
os.makedirs(repoDir, exist_ok=True)
os.makedirs(resultDir, exist_ok=True)
os.makedirs(objectDir, exist_ok=True)
os.makedirs(figDir, exist_ok=True)
os.makedirs(tabDir, exist_ok=True)


# inDir for results from previous step
inDir = os.path.join("/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/02_cell2location/06a_cell_type_mapping_major", section_name, "objects")

# --- Extract some info from the sample sheet

# Extract the sample_id
sample_dir = os.path.join(repoDir, 'config')
sample_sheet_file = 'sample_sheet.csv'

sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))

# Filter the DataFrame based on the condition sample_id == sample_id
filtered_df = sample_sheet[sample_sheet['section_name'] == section_name]

# Extract the sample_name from the first row (assuming there's only one match)
file_id = filtered_df['sample_name'].iloc[0]


# --- Load the cell2location output

# Load results
adata_file = f"{file_id}_cell2location.h5ad"
adata_vis = sc.read_h5ad(os.path.join(inDir, adata_file))


# --- Visualizing cell abundance in spatial coordinates

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

# Extract cell type names
cell_types = adata_vis.uns['mod']['factor_names']

with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=cell_types,
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )

plt.savefig(os.path.join(figDir, 'cell_type_abundance_spatial.pdf'))
plt.close('all')

# Now we use cell2location plotter that allows showing multiple cell types in one panel
from cell2location.plt import plot_spatial

# select up to 6 clusters
clust_labels = ['CAFs', 'SMCs', 'Epithelial', 'PNS_Glial']
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=adata_vis,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
        show_img=True,
        # 'fast' (white background) or 'dark_background'
        style='fast',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=6,
        colorbar_position='right'
    )

plt.savefig(os.path.join(figDir, 'selected_cell_types.pdf'))
plt.close('all')


print('Script succesfully completed')