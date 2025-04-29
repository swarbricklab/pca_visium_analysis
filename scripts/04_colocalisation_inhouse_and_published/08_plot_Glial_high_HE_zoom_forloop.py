#!/usr/bin/env python3

# --- Aim
# Zoom in on regions that have a high signature score
# Plot high-res H&E to check that nerves are visible

# --- Define environment
# env: spatial-pilot

# --- Load packages
import argparse
import scanpy as sc
from matplotlib_scalebar.scalebar import ScaleBar
import anndata as ad
import os
import pandas as pd
import matplotlib.pyplot as plt
from datetime import date
from datetime import datetime
import PIL
import glob
import numpy as np

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as patches
from mycolorpy import colorlist as mcp


# Import custom functions
from sophiesfunctions.auto_crop import auto_crop_scanpy  # import my custom functions that are packaged locally as sophiesfunctions

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3


# --- Directories

project_dir = "/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/"
h5ad_dir = "/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/data/anndata_objects/log_norm/"
res_dir = "/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/04_colocalisation_inhouse_and_published/"
score_dir = "/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/04_colocalisation_inhouse_and_published/05_zoom_in_HandE/csv_files/PNS_Glial_score_per_sample/20240730/"

# --- Get some info from the sample sheet

# Extract the external_id
sample_dir = os.path.join(project_dir, 'config')
sample_sheet_file = 'sample_sheet.csv'
sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))

# Get unique section names for batches 2 and 3
section_names = sample_sheet[sample_sheet['batch'].isin([2, 3])]['section_name'].unique()

section_name = "A1_P11_C1"

slide_size = sample_sheet[sample_sheet['section_name'] == section_name]['slide_size'].values[0]

# --- Load the adata object

adata = ad.read_h5ad(os.path.join(h5ad_dir, f'{section_name}_log_norm_meta.h5ad'))


# --- Calculate crop coords for optimal presentation
coords = auto_crop_scanpy(adata, border_size=0.1)


# --- Set spot size for large & small capture area
# Plotting the spots larger to make it easier to see

if slide_size == '11 mm':
    spot_size = 1.8
else:
    spot_size = 1.3


# --- Also load the group results

threshold = 0.5

# threshold_adj = str(threshold).split('.')[-1]

TLS_group_dir = os.path.join(score_dir)

# # Find most recent results
# subdirs = [d for d in os.listdir(TLS_group_dir) if os.path.isdir(os.path.join(TLS_group_dir, d))]
# 
# valid_subdirs = []
# for subdir in subdirs:
#     try:
#         datetime.strptime(subdir, date_format)
#         valid_subdirs.append(subdir)
#     except ValueError:
#         pass
# 
# # Sort the valid subdirectories by date and get the most recent one
# if valid_subdirs:
#     most_recent_dir = max(valid_subdirs, key=lambda x: datetime.strptime(x, date_format))
#     most_recent_dir_path = os.path.join(TLS_group_dir, most_recent_dir)
#     print("Most recent results:", most_recent_dir_path)
# else:
#     print("No valid directories found.")

# Read in the results & add to anndata
TLS_file_name = f'{section_name}_PNS_glial_groups.csv'

group_key = 'group_singleton_removed'

TLS_results = pd.read_csv(os.path.join(TLS_group_dir, TLS_file_name), index_col=0)

adata.obs[group_key] = TLS_results[group_key]

# --- Define custom colormap for plotting TLS model results
# Set everything under certain threshold as the lowest color of a color map

cmap = "YlOrRd"
color_list = mcp.gen_color(cmap=cmap, n=11)
color_list_extracted = ['lightgray'] + color_list[2:]

n_bins = 100
cmap_name = 'my_list'
cmap_adj = LinearSegmentedColormap.from_list(cmap_name, color_list_extracted, N=n_bins)

cmap_adj.set_under(color='lightgray', alpha=0.4)


# Subset anndata to PNS_glial score above threshold
adata_sub = adata[adata.obs['PNS_glial_score'] > threshold]


# --- Zoom into the groups: import the original high-resolution scan

# Because the resolution of 'hires' isn't great, everything becomes blurry once you zoom in
# Therefore use the original image for plotting the zoomed in plots

# Find the original high-res image

# new samples:
# /share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_eva/data/processed/spaceranger_outs/PCaP10_C1_P3_C1/outs/spatial
#  tissue_hires_image.png

# old samples:
# /share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_eva/data/processed/spaceranger_outs/A1_P11_C1/outs/spatial
# tissue_hires_image.png


image_dir = os.path.join(project_dir, 'data', 'images')
image_file_path = glob.glob(os.path.join(image_dir, f'{section_name}*.tif'))[0]

# Need to set this or it won't let me import my image
PIL.Image.MAX_IMAGE_PIXELS = 24078974976

# Load the original high-res scan
img = plt.imread(image_file_path)

# Define key for adding to anndata
original_hires_key = 'original_hires'

# -- Calculate the scale factor
# Divide the height of the original hires image by the height of the downsampled hires image
size_fact = img.shape[0]/adata_sub.uns['spatial'][section_name]['images']['hires'].shape[0]

# Extract the downsampled hires scale factor
hires_fact = adata_sub.uns['spatial'][section_name]['scalefactors']['tissue_hires_scalef']

# Calculate the scale factor for the original hires image and add to the anndata object
adata_sub.uns['spatial'][section_name]['scalefactors'][f'tissue_{original_hires_key}_scalef'] = size_fact * hires_fact

# -- Add scale factor to original anndata as well
adata.uns['spatial'][section_name]['scalefactors'][f'tissue_{original_hires_key}_scalef'] = size_fact * hires_fact


# --- Plot 1: Zoom into groups and plot on top of original high-res H&E

sc.set_figure_params(facecolor="white", figsize=(8, 8))

adata.obs['group_singleton_removed'] = adata_sub.obs['group_singleton_removed']

all_groups = adata_sub.obs['group_singleton_removed'].astype('category').unique()
all_groups = all_groups.categories

# Make sure they're in the correct order
# Sort the groups by number
all_groups = sorted(all_groups, key=lambda x: int(x[5:]))  # Sort by the number part of the group name
print(all_groups)

# Drop 'singleton' from the list
# all_groups = [x for x in all_groups if x != 'singleton']

ncols = 6
nrows = 1

sc.set_figure_params(fontsize=18)

# Get the scalefactor for the hires image
# To be able to plot a rectangle marking the zoom area
scalef = adata.uns['spatial'][section_name]['scalefactors']['tissue_hires_scalef']

# Specify vmax for TLS model
vmax = 0.5

# --- Prepare for saving the plots

# Put in a timestamped folder to avoid overwriting older plots
qc_dir = os.path.join(res_dir, 'PNS_glial_score')
os.makedirs(qc_dir, exist_ok=True)

os.chdir(qc_dir)

today = date.today()
today = today.strftime("%Y%m%d")

out_dir = os.path.join(qc_dir, 'pnCAF_plots_per_sample')
os.makedirs(out_dir, exist_ok=True)

# Save plots
out_dir_plots = os.path.join(out_dir, 'zoom', today)
os.makedirs(out_dir_plots, exist_ok=True)

# Specify a filename for saving the plot
filename = f'{section_name}_PNS_glial_groups_zoom.pdf'

# --- Plotting

adata_sub.obs['in_tissue'] = adata_sub.obs['in_tissue'].astype('category')

# Calculate microns per pixel for the scalebar
# Get the diameter of the spots in the fullres image
spot_diameter_fulres = adata.uns['spatial'][section_name]['scalefactors']['spot_diameter_fullres']

# Spots are 55 µm
spot_diameter = 55

microns_per_pixel = spot_diameter/spot_diameter_fulres

# Function to crop the image with safe coordinates
# Bc. sometimes the calculated crop coordinates might be outside of the image area (especially if the original image has been cropped)
def crop_image(img, coords):
    # Clip the coordinates to ensure they are within the image boundaries
    x_min = max(coords[0], 0)
    x_max = min(coords[1], img.shape[1])
    y_min = max(coords[2], 0)
    y_max = min(coords[3], img.shape[0])
    # Crop the image using the clipped coordinates
    crop_img = img[y_min:y_max, x_min:x_max]
    return crop_img, (x_min, x_max, y_min, y_max)

# Plot the subplots for each group, saving and closing each plot as you go (hopefully easier on the memory)
with PdfPages(os.path.join(out_dir_plots, filename)) as pdf:
    for i, group in enumerate(all_groups):
        fig, axs = plt.subplots(nrows, ncols, figsize=(8*ncols, 6*nrows))
        # Subset the anndata
        adata_group_sub = adata_sub[adata_sub.obs['group_singleton_removed'] == group]
        # Make sure border_size is equal despite differences in area
        desired_border = spot_diameter_fulres * 3
        coords_group = auto_crop_scanpy(adata_group_sub, border_size=0)
        width = coords_group[1] - coords_group[0]
        border_size=desired_border/width
        # Calculate new crop coord based on the group
        coords_group = auto_crop_scanpy(adata_group_sub, border_size=border_size)
        # -- Plot the full tissue with a square marking the zoom area (TLS model results)
        n=0
        ax=axs[n]
        # Plot the group on H&E hires image
        sc.pl.spatial(adata, color='PNS_glial_score', groups=[group], bw=False, alpha_img=1, crop_coord=coords, title=f'{section_name} - {group}',
                    na_in_legend=False, ax=ax, cmap=cmap_adj, vmin=0, vmax=vmax)
        # Plot box marking the zoom area
        # Adjust coordinates for the hires image
        coords_group_adj = tuple(int(x * scalef) for x in coords_group)
        x1, x2, y1, y2 = coords_group_adj
        width = x2 - x1
        height = y2 - y1
        rect = patches.Rectangle((x1, y1), width, height, linewidth=2, edgecolor='black', facecolor='none', linestyle=(0, (0.5,0.5)))
        # Add the rectangle to the plot
        ax.add_patch(rect)
        # -- Plot the full tissue with a square marking the zoom area
        n+=1
        ax=axs[n]
        # Plot the group on H&E hires image
        sc.pl.spatial(adata_group_sub, color='group_singleton_removed', groups=[group], bw=False, alpha_img=1, crop_coord=coords, title=f'{section_name} - {group}',
                    na_in_legend=False, ax=ax,)
        # Remove legend
        ax.get_legend().remove()
        # Plot box marking the zoom area
        rect = patches.Rectangle((x1, y1), width, height, linewidth=2, edgecolor='black', facecolor='none', linestyle=(0, (0.5,0.5)))
        # Add the rectangle to the plot
        ax.add_patch(rect)
        # -- Prepare an anndata that contains only the cropped area we want to plot
        # This will drastically speed things up
        # Crop the H&E image
        # Check if any of the coords are under 0, need to set to 0 if they are
        crop_img, adjusted_coords = crop_image(img, coords_group)
        # Calculate the height and width of the cropped image
        height_crop = adjusted_coords[3] - adjusted_coords[2]
        width_crop = adjusted_coords[1] - adjusted_coords[0]
        # Also check if either x0 or y0 are below 0
        xmin = 0
        if adjusted_coords[0] < xmin:
            xmin = adjusted_coords[0]
        ymin = 0
        if adjusted_coords[2] < ymin:
            ymin = adjusted_coords[2]
        # Calculate the width & height of the crop
        height_crop = coords_group[3] - coords_group[2]
        width_crop = coords_group[1] - coords_group[0]
        # Make a copy of anndata & add the crop in the hires slot
        adata_crop = adata.copy()
        adata_crop.uns['spatial'][section_name]['images'][original_hires_key] = crop_img
        # Shift the spot coordinates
        spot_coords = adata_crop.obsm['spatial']
        shifted_coords = spot_coords - np.array([adjusted_coords[0], adjusted_coords[2]])
        adata_crop.obsm['spatial'] = shifted_coords
        # -- Plot the spots belonging to the group, zoomed in
        n+=1
        ax=axs[n]
        sc.pl.spatial(adata_crop, color='group_singleton_removed', groups=[group], bw=False, alpha_img=1, crop_coord=(xmin, width_crop, ymin, height_crop), title=f'{section_name} - {group}',
                    na_in_legend=False, ax=ax, img_key=original_hires_key)
        # Remove legend
        ax.get_legend().remove()
        # Add the scalebar
        scalebar = ScaleBar(microns_per_pixel, units='µm', location='lower right', box_alpha=0.5)
        ax.add_artist(scalebar)
        # -- Plot the PNS_glial score within the zoom area
        n+=1
        ax=axs[n]
        sc.pl.spatial(adata_crop, color='PNS_glial_score', bw=False, alpha_img=1, crop_coord=(xmin, width_crop, ymin, height_crop), ax=ax, cmap=cmap_adj, vmin=0, vmax=vmax, title=f'{group} - PNS_glial score',
                    img_key=original_hires_key)
        # Add the scalebar
        scalebar = ScaleBar(microns_per_pixel, units='µm', location='lower right', box_alpha=0.5)
        ax.add_artist(scalebar)
        # -- Plot the H&E with nothing over it in the zoom area
        n+=1
        ax=axs[n]
        sc.pl.spatial(adata_crop, color=None, bw=False, alpha_img=1, crop_coord=(xmin, width_crop, ymin, height_crop), ax=ax, cmap=cmap_adj, vmin=0, title=f'{group} - H&E',
                    img_key=original_hires_key)
        # Add the scalebar
        scalebar = ScaleBar(microns_per_pixel, units='µm', location='lower right', box_alpha=0.5)
        ax.add_artist(scalebar)
        # -- Final one: open circle spots
        # Subset
        adata_crop_sub = adata_crop[adata_crop.obs['group_singleton_removed'] == group]
        adata_crop_sub.obs['in_tissue'] = adata_crop_sub.obs['in_tissue'].astype('category')
        n+=1
        ax=axs[n]
        sc.pl.spatial(adata_crop_sub, color='in_tissue', bw=False, alpha_img=1, crop_coord=(xmin, width_crop, ymin, height_crop), ax=ax, title=f'{group} - spots',
                    img_key=original_hires_key)
        # Remove legend
        ax.get_legend().remove()
        # Plot spots as open circles
        for collection in ax.collections:
            collection.set_facecolor((0, 0, 0, 0))
            collection.set_edgecolor((0, 0, 0, 0.7))
            collection.set_linewidth(2)
        # Add the scalebar
        scalebar = ScaleBar(microns_per_pixel, units='µm', location='lower right', box_alpha=0.5)
        ax.add_artist(scalebar)
        plt.suptitle(f'{section_name} - grouping of spots with PNS_glial score >{threshold}) - {group}', y=1.1, fontweight='bold', x=0.5)
        # Save the figure to the PdfPages object
        pdf.savefig(bbox_inches='tight')
        # Close the plot
        plt.close()


# --- Finish

print("Script successfully completed")
