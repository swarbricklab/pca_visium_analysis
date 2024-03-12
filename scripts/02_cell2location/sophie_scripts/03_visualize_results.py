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

parser = argparse.ArgumentParser(description="Import Visium data & save as .h5ad file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--raw_data_dir", help="Path to raw data")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--min_count", help="Minimal number of counts per spot", type=int)
parser.add_argument("--min_gene", help="Minimal number of genes per spot", type=int)

args = vars(parser.parse_args())


# --- Arguments
sample_id=args["sample_id"]
project_dir=args["project_dir"]


# --- Extract some info from the sample sheet

# Extract the section_id
sample_dir = os.path.join(project_dir, 'config')
sample_sheet_file = 'sample_sheet.csv'

sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))
external_id = sample_sheet[sample_sheet['section_id'] == sample_id]['external_id'].values[0]

joined_id = f'{external_id}-{sample_id}'


# --- Load the cell2location output

res_dir = os.path.join(project_dir, 'results')

qc_dir = os.path.join(res_dir, 'cell2location', sample_id)
os.makedirs(qc_dir, exist_ok=True)

# Find most recent results
subdirs = [d for d in os.listdir(qc_dir) if os.path.isdir(os.path.join(qc_dir, d))]

# Find most recent results
# Filter subdirectories to keep only those with date format YYYYmmdd
date_format = '%Y%m%d'
valid_subdirs = []

for subdir in subdirs:
    try:
        datetime.strptime(subdir, date_format)
        valid_subdirs.append(subdir)
    except ValueError:
        pass

# Sort the valid subdirectories by date and get the most recent one
if valid_subdirs:
    most_recent_dir = max(valid_subdirs, key=lambda x: datetime.strptime(x, date_format))
    most_recent_dir_path = os.path.join(qc_dir, most_recent_dir)
    print("Most recent results:", most_recent_dir_path)
else:
    print("No valid directories found.")


out_dir = most_recent_dir_path
os.makedirs(out_dir, exist_ok=True)

# Also make an output directory for the plots
out_dir_plots = os.path.join(out_dir, 'plots')
os.makedirs(out_dir_plots, exist_ok=True)

# Load results
adata_file = f"{joined_id}_cell2location.h5ad"
adata_vis = sc.read_h5ad(os.path.join(out_dir, 'model_output', adata_file))


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

plt.savefig(os.path.join(out_dir_plots, 'cell_type_abundance_spatial.pdf'))
plt.close('all')

# Now we use cell2location plotter that allows showing multiple cell types in one panel
from cell2location.plt import plot_spatial

# select up to 6 clusters
clust_labels = ['CAFs myCAF-like', 'B cells Naive', 'T cells CD4+']
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

plt.savefig(os.path.join(out_dir_plots, 'selected_cell_types.pdf'))
plt.close('all')


print('Script succesfully completed')