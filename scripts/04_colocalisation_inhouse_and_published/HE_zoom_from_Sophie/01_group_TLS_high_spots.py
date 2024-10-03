#!/usr/bin/env python3

# --- Aim
# Zoom in on regions that have a high TLS score
# Plot high-res H&E to check that immune aggregates are visible


# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import scanpy as sc
import anndata as ad
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
from datetime import datetime
import squidpy as sq

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
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--project_dir", help="Project directory (root)")
parser.add_argument("--sample_id", help="Sample ID")

parser.add_argument("--min_count", help="Min count for QC", type=int)
parser.add_argument("--min_gene", help="Min gene for QC", type=int)
parser.add_argument("--min_spots", help="Minimal spots per gene", type=int)

args = vars(parser.parse_args())


# --- Arguments

res_dir = args["res_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]


# --- Get some info from the sample sheet

# Extract the external_id
sample_dir = os.path.join(project_dir, 'config')
sample_sheet_file = 'sample_sheet.csv'

sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))
external_id = sample_sheet[sample_sheet['section_id'] == sample_id]['external_id'].values[0]
tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]
slide_size = sample_sheet[sample_sheet['section_id'] == sample_id]['slide_size'].values[0]

joined_id = f'{external_id}-{sample_id}'

# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Create a sample - subtype dictionary
subtype = tissue_metadata[tissue_metadata['tissue_id'] == tissue_id]['clinical_subtype'].values[0]

# --- Load the adata object

use_data = 'log_norm'

adata = ad.read_h5ad(os.path.join(h5ad_dir, 'raw', f'{joined_id}_raw.h5ad'))


# --- Calculate crop coords for optimal presentation
coords = auto_crop_scanpy(adata, border_size=0.1)


# --- Set spot size for large & small capture area
# Plotting the spots larger to make it easier to see

if slide_size == '11 mm':
    spot_size = 1.8
else:
    spot_size = 1.3


# --- Load TLS results

TLS_res_dir = os.path.join(res_dir, 'TLS_model')

# Find most recent results
subdirs = [d for d in os.listdir(TLS_res_dir) if os.path.isdir(os.path.join(TLS_res_dir, d))]

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
    most_recent_dir_path = os.path.join(TLS_res_dir, most_recent_dir)
    print("Most recent results:", most_recent_dir_path)
else:
    print("No valid directories found.")

# Read in the results & add to anndata

TLS_file_name = f'{sample_id}_TLS_model.csv'
TLS_results = pd.read_csv(os.path.join(most_recent_dir_path, TLS_file_name), index_col=0)
adata.obs['TLS_model'] = TLS_results


# --- Prepare for saving the plots

# Put in a timestamped folder to avoid overwriting older plots
qc_dir = os.path.join(res_dir, 'TLS_model')
os.makedirs(qc_dir, exist_ok=True)

os.chdir(qc_dir)

today = date.today()
today = today.strftime("%Y%m%d")


# --- Plot 1: Plot the TLS model score + spots above a certain threshold

# Define custom colormap for plotting TLS model results
# Set everything under certain threshold as the lowest color of a color map

cmap = "YlOrRd"
color_list = mcp.gen_color(cmap=cmap, n=11)
color_list_extracted = ['lightgray'] + color_list[2:]

n_bins = 100
cmap_name = 'my_list'
cmap_adj = LinearSegmentedColormap.from_list(cmap_name, color_list_extracted, N=n_bins)

cmap_adj.set_under(color='lightgray', alpha=0.4)

# Also define a threshold for selecting TLS-high spots
threshold = 0.01
# Subset anndata to TLS score above threshold
adata_sub = adata[adata.obs['TLS_model'] > threshold]

# Plotting
sc.set_figure_params(fontsize=16)

nrows = 1
ncols = 2

fig, axs = plt.subplots(nrows, ncols, figsize=(8*ncols, 6.5*nrows))

# TLS model results
sc.pl.spatial(adata, color='TLS_model', size=spot_size, cmap=cmap_adj, bw=False, alpha_img=1, crop_coord=coords, title=f'TLS score', vmin=0, ax=axs[0])

# Spots above threshold
adata_sub.obs['in_tissue'] = adata_sub.obs['in_tissue'].astype('category')
sc.pl.spatial(adata_sub, color='in_tissue', size=spot_size, palette=['yellow'], bw=False, alpha_img=1, crop_coord=coords, title=f'Spots >{threshold}', vmin=0, 
            ax=axs[1])

# Remove the legend
axs[1].get_legend().set_visible(False)

plt.suptitle(f'{sample_id} ({subtype}) - TLS score + TLS-high spots', y=1.05, fontweight='bold', x=0.5)


# --- If there are no spots over the threshold: need to stop here

threshold_adj = str(threshold).split('.')[-1]

if len(adata_sub) == 0 or len(adata_sub) == 1:
    out_dir = os.path.join(qc_dir, 'TLS_plots_per_sample', f'threshold_{threshold_adj}')
    os.makedirs(out_dir, exist_ok=True)
    # Save plots
    out_dir_plots = os.path.join(out_dir, 'plots', 'groups', today)
    os.makedirs(out_dir_plots, exist_ok=True)
    filename = f'{sample_id}_TLS_groups.pdf'
    save_multi_image(os.path.join(out_dir_plots, filename))
    # Also save the groupings as a .csv
    out_dir_csv = os.path.join(out_dir, 'csv_files', today)
    os.makedirs(out_dir_csv, exist_ok=True)
    to_save_df = adata_sub.obs[['TLS_model']]
    to_save_df['group'] = ''
    to_save_df['group_singleton'] = ''
    filename = f'{sample_id}_TLS_groups.csv'
    to_save_df.to_csv(os.path.join(out_dir_csv, filename))
    print('Script succesfully completed')


# --- Define groups of neighboring spots

# Get a scalefactor for converting pixel coords to um
spot_diameter_pixel = adata.uns['spatial'][sample_id]['scalefactors']['spot_diameter_fullres']
spot_diameter_um = 55 # Regular Visium spots: 55 µm 
scalefactor_um = spot_diameter_pixel/spot_diameter_um

# Regular Visium spots are 100 um apart
spot_dist_um = 100
spot_dist_pixel = spot_dist_um * scalefactor_um

# Decide which radius I want to use for grouping spots
num_spots = 1 # How many spots away is part of same group
radius_max = spot_dist_pixel*num_spots

# Use squidpy to define which TLS-high spots are neighbors
sq.gr.spatial_neighbors(adata_sub, coord_type='generic', radius=(0, radius_max))

adata_sub.obs['in_tissue'] = adata_sub.obs['in_tissue'].astype('category')

# --- Plot 2: Plot scanpy's neighborhood graph result

if slide_size == '11 mm':
    squidpy_size = 4
    squidpy_edge = 0.5
else:
    squidpy_size = 30
    squidpy_edge = 1

sc.set_figure_params(fontsize=15, figsize=(6,6))

sq.pl.spatial_scatter(
    adata_sub,
    color="in_tissue",
    connectivity_key="spatial_connectivities",
    edges_color="black",
    shape=None,
    edges_width=squidpy_edge,
    size=squidpy_size,
    palette='Dark2'
)

plt.suptitle(f'{sample_id} ({subtype}) - neighborhood graph of spots with TLS score >{threshold})', y=1.1, fontweight='bold', x=0.5)


# --- Now annotate the spots as groups

# Extract connected spots as groups
connected_groups = []
for i, neighbors in enumerate(adata_sub.obsp['spatial_connectivities'].toarray()):
    group = [j for j, is_neighbor in enumerate(neighbors) if is_neighbor]
    group.append(i)  # Include the spot itself in the group
    connected_groups.append(group)

# squidpy's spatial_connectivities just saves which spots are neighbors
# so there will be overlap between the sublists (if a spot is a neighbor of multiple other spots, it appears multiple times)
# therefore merge the sublists until all values are unique

pooled = [set(subList) for subList in connected_groups]
merging = True
while merging:
    merging=False
    for i,group in enumerate(pooled):
        merged = next((g for g in pooled[i+1:] if g.intersection(group)),None)
        if not merged: continue
        group.update(merged)
        pooled.remove(merged)
        merging = True

# Convert the list of sets to a list of lists
list_of_lists = [list(s) for s in pooled]

# Generate group names
groups = [f'group{i+1}' for i in range(len(list_of_lists))]

# Flatten the list of lists and pair each item with its corresponding group
data_list = []
for group, sublist in zip(groups, list_of_lists):
    for item in sublist:
        data_list.append({'index': item, 'group': group})

# Create a DataFrame
df = pd.DataFrame(data_list)

# Sort DataFrame based on the 'index' column
df_sorted = df.sort_values(by='index')

# Transfer the barcode id to the index from the anndata object
df_sorted.index = adata_sub.obs_names

# Tack the group annotation onto the anndata object
adata_sub.obs['group'] = df_sorted['group']


# --- Plot 3: Plot the grouping of high TLS spots

sc.set_figure_params(fontsize=12, figsize=(6,6))

sc.pl.spatial(adata_sub, color='group', size=spot_size, bw=False, alpha_img=1, crop_coord=coords, title=f'{sample_id} - groups of spots >{threshold}', )

plt.suptitle(f'{sample_id} ({subtype}) - grouping of spots with TLS score >{threshold})', y=1.05, fontweight='bold', x=0.5)


# --- Remove groups consisting of single spots

# Find groups that appear only once
singleton_groups = adata_sub.obs['group'].value_counts()[adata_sub.obs['group'].value_counts() == 1].index

# Make a temporary df
temp_df = adata_sub.obs[['group']]

# Set singleton groups to 'singleton'
temp_df['group'] = np.where(temp_df['group'].isin(singleton_groups), 'singleton', temp_df['group'])

# Save singletons as separat df
singletons = temp_df[temp_df['group'] == 'singleton']

# Remove the singletons for now
temp_df = temp_df[temp_df['group'] != 'singleton']

# Update group names based on size
group_counts = temp_df['group'].value_counts()
sorted_groups = group_counts.sort_values(ascending=False).index
group_map = {group: f"group{i+1}" for i, group in enumerate(sorted_groups)}

temp_df['group'] = temp_df['group'].map(group_map)

# Add results to a new column
adata_sub.obs['TLS_group'] = temp_df['group']

# Add groups + singletons to another column
group_singleton = pd.concat([temp_df, singletons])

adata_sub.obs['group_singleton'] = group_singleton


# --- Plot 4: Groups after removing singletons

sc.pl.spatial(adata_sub, color='group_singleton', size=spot_size, bw=False, alpha_img=1, crop_coord=coords, title=f'{sample_id} - TLS-high spots', )

plt.suptitle(f'{sample_id} ({subtype}) - grouping of spots with TLS score >{threshold} - singletons removed', y=1.05, fontweight='bold', x=0.5)


# --- Save everything so far

threshold_adj = str(threshold).split('.')[-1]

out_dir = os.path.join(qc_dir, 'TLS_plots_per_sample', f'threshold_{threshold_adj}')
os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'plots', 'groups', today)
os.makedirs(out_dir_plots, exist_ok=True)

filename = f'{sample_id}_TLS_groups.pdf'

save_multi_image(os.path.join(out_dir_plots, filename))


# Also save the groupings as a .csv
out_dir_csv = os.path.join(out_dir, 'csv_files', today)
os.makedirs(out_dir_csv, exist_ok=True)

to_save_df = adata_sub.obs[['TLS_model', 'group', 'group_singleton', 'TLS_group']]

# Rename columns
to_save_df.columns = ['TLS_model', 'all_groups', 'group_singleton', 'group_singleton_removed']

filename = f'{sample_id}_TLS_groups.csv'
to_save_df.to_csv(os.path.join(out_dir_csv, filename))


# --- Finish

print("Script successfully completed")
