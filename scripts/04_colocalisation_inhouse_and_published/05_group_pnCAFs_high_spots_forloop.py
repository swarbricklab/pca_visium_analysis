#!/usr/bin/env python3

# --- Aim
# Zoom in on regions that have a high signature score
# Plot high-res H&E to check that nerves are visible

# need to run this in 'scripts', conda env fails in root directory ('pca_visium_analysis')

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

from matplotlib.backends.backend_pdf import PdfPages  # Import PdfPages

# Import custom functions
from sophiesfunctions.auto_crop import auto_crop_scanpy  # import my custom functions that are packaged locally as sophiesfunctions
from sophiesfunctions.save_multi_image import save_multi_image

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

# --- Directories

project_dir = "/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/"
h5ad_dir = "/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/data/anndata_objects/log_norm/"
res_dir = "/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/04_colocalisation_inhouse_and_published/"

# --- Get some info from the sample sheet

# Extract the external_id
sample_dir = os.path.join(project_dir, 'config')
sample_sheet_file = 'sample_sheet.csv'
sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))

# Get unique section names for batches 2 and 3
section_names = sample_sheet[sample_sheet['batch'].isin([2, 3])]['section_name'].unique()

# List files in h5ad_dir
files_in_h5ad_dir = os.listdir(h5ad_dir)

# Filter section names to only keep those present in files with 'log_norm_meta.h5ad' suffix
filtered_section_names = [name for name in section_names if f"{name}_log_norm_meta.h5ad" in files_in_h5ad_dir]
# Exclude 'A1_P25' and filter section names based on available files
# filtering out these section names manually - the for loop is failing, perhaps bc they dont have enough spots above threshold??
filtered_section_names = [name for name in section_names if name != 'A1_P25' and name != 'B1_P15_C2' and name != 'D1_P4_C2' and f"{name}_log_norm_meta.h5ad" in files_in_h5ad_dir]

# Print the filtered section names
print(filtered_section_names)

# separate pdf for each section_name
# spaces removed bc I couldn't run the script interactively otherwise
# output directory is a bit messy when scores above threshold are not detected

for section_name in filtered_section_names:
    # Extract external_id, tissue_id, and slide_size
    slide_size = sample_sheet[sample_sheet['section_name'] == section_name]['slide_size'].values[0]
    # Load the adata object
    adata = ad.read_h5ad(os.path.join(h5ad_dir, f'{section_name}_log_norm_meta.h5ad'))
    # Calculate crop coords for optimal presentation
    coords = auto_crop_scanpy(adata, border_size=0.1)
    # Set spot size based on slide size
    spot_size = 1.8 if slide_size == '11 mm' else 1.3
    # Prepare for saving the plots
    today = date.today().strftime("%Y%m%d") 
    qc_dir = os.path.join(res_dir, '05_zoom_in_HandE', 'figures', 'pnCAFs_score_per_sample', today)
    os.makedirs(qc_dir, exist_ok=True)
    # Define PDF file path
    pdf_path = os.path.join(qc_dir, f'{section_name}_combined_plots.pdf')
    with PdfPages(pdf_path) as pdf:
        # Plot 1: Plot the pnCAFs_score model score + spots above a certain threshold
        cmap = "YlOrRd"
        color_list = mcp.gen_color(cmap=cmap, n=11)
        color_list_extracted = ['lightgray'] + color_list[2:]
        n_bins = 100
        cmap_name = 'my_list'
        cmap_adj = LinearSegmentedColormap.from_list(cmap_name, color_list_extracted, N=n_bins)
        cmap_adj.set_under(color='lightgray', alpha=0.4)
        threshold = 0.5
        adata_sub = adata[adata.obs['pnCAFs_score'] > threshold]
        sc.set_figure_params(fontsize=16)
        nrows = 1
        ncols = 2
        fig, axs = plt.subplots(nrows, ncols, figsize=(8*ncols, 6.5*nrows))
        # TLS model results
        sc.pl.spatial(adata, color='pnCAFs_score', size=spot_size, cmap=cmap_adj, bw=False, alpha_img=1, crop_coord=coords, title=f'pnCAFs_score', vmin=0, ax=axs[0])
        # Spots above threshold
        adata_sub.obs['in_tissue'] = adata_sub.obs['in_tissue'].astype('category')
        sc.pl.spatial(adata_sub, color='in_tissue', size=spot_size, palette=['yellow'], bw=False, alpha_img=1, crop_coord=coords, title=f'Spots >{threshold}', vmin=0, ax=axs[1])
        # Remove the legend
        axs[1].get_legend().set_visible(False)
        plt.suptitle(f'{section_name} - pnCAFs_score + pnCAFs-high spots', y=1.05, fontweight='bold', x=0.5)
        # Save the first plot to the PDF
        pdf.savefig(fig)
        plt.close(fig)
        # Save plots if there are no spots over the threshold
        threshold_adj = str(threshold).split('.')[-1]
        if len(adata_sub) == 0 or len(adata_sub) == 1:
            out_dir = os.path.join(qc_dir, f'threshold_{threshold_adj}')
            os.makedirs(out_dir, exist_ok=True)
            # Also save the groupings as a .csv
            out_dir_csv = os.path.join(res_dir, '05_zoom_in_HandE', 'csv_files', 'pnCAFs_score_per_sample', today)
            os.makedirs(out_dir_csv, exist_ok=True)
            to_save_df = adata_sub.obs[['pnCAFs_score']]
            to_save_df['group'] = ''
            to_save_df['group_singleton'] = ''
            filename = f'{section_name}_pnCAFs_groups.csv'
            to_save_df.to_csv(os.path.join(out_dir_csv, filename))
            print(f'Successfully completed for {section_name}')  # Updated print statement
            continue  # Skip to the next section_name
        # Define groups of neighboring spots
        spot_diameter_pixel = adata.uns['spatial'][section_name]['scalefactors']['spot_diameter_fullres']
        spot_diameter_um = 55  # Regular Visium spots: 55 µm
        scalefactor_um = spot_diameter_pixel / spot_diameter_um
        spot_dist_um = 100
        spot_dist_pixel = spot_dist_um * scalefactor_um
        num_spots = 1  # How many spots away is part of the same group
        radius_max = spot_dist_pixel * num_spots
        sq.gr.spatial_neighbors(adata_sub, coord_type='generic', radius=(0, radius_max))
        adata_sub.obs['in_tissue'] = adata_sub.obs['in_tissue'].astype('category')
        # Plot 2: Plot scanpy's neighborhood graph result
        squidpy_size = 4 if slide_size == '11 mm' else 30
        squidpy_edge = 0.5 if slide_size == '11 mm' else 1
        sc.set_figure_params(fontsize=15, figsize=(6,6))
        fig, ax = plt.subplots(figsize=(6,6))
        sq.pl.spatial_scatter(
            adata_sub,
            color="in_tissue",
            connectivity_key="spatial_connectivities",
            edges_color="black",
            shape=None,
            edges_width=squidpy_edge,
            size=squidpy_size,
            palette='Dark2',
            ax=ax
        )
        plt.suptitle(f'{section_name} - neighborhood graph of spots with pnCAFs score >{threshold})', y=1.1, fontweight='bold', x=0.5)
        # Save the second plot to the PDF
        pdf.savefig(fig)
        plt.close(fig)
        # Now annotate the spots as groups
        connected_groups = []
        for i, neighbors in enumerate(adata_sub.obsp['spatial_connectivities'].toarray()):
            group = [j for j, is_neighbor in enumerate(neighbors) if is_neighbor]
            group.append(i)  # Include the spot itself in the group
            connected_groups.append(group)
        pooled = [set(subList) for subList in connected_groups]
        merging = True
        while merging:
            merging = False
            for i, group in enumerate(pooled):
                merged = next((g for g in pooled[i+1:] if g.intersection(group)), None)
                if not merged:
                    continue
                group.update(merged)
                pooled.remove(merged)
                merging = True
        list_of_lists = [list(s) for s in pooled]
        groups = [f'group{i+1}' for i in range(len(list_of_lists))]
        data_list = []
        for group, sublist in zip(groups, list_of_lists):
            for item in sublist:
                data_list.append({'index': item, 'group': group})
        df = pd.DataFrame(data_list)
        df_sorted = df.sort_values(by='index')
        df_sorted.index = adata_sub.obs_names
        adata_sub.obs['group'] = df_sorted['group']
        # Plot 3: Plot the grouping of high TLS spots
        sc.set_figure_params(fontsize=12, figsize=(6,6))
        fig, ax = plt.subplots(figsize=(6,6))
        sc.pl.spatial(adata_sub, color='group', size=spot_size, bw=False, alpha_img=1, crop_coord=coords, title=f'{section_name} - groups of spots >{threshold}', ax=ax)
        plt.suptitle(f'{section_name} - grouping of spots with pnCAFs_score >{threshold})', y=1.05, fontweight='bold', x=0.5)
        # Save the third plot to the PDF
        pdf.savefig(fig)
        plt.close(fig)
 		# Remove groups consisting of single spots
        singleton_groups = adata_sub.obs['group'].value_counts()[adata_sub.obs['group'].value_counts() == 1].index
        temp_df = adata_sub.obs[['group']]
        temp_df['group'] = np.where(temp_df['group'].isin(singleton_groups), 'singleton', temp_df['group'])
        singletons = temp_df[temp_df['group'] == 'singleton']
        temp_df = temp_df[temp_df['group'] != 'singleton']
        group_counts = temp_df['group'].value_counts()
        sorted_groups = group_counts.sort_values(ascending=False).index
        group_map = {group: f"group{i+1}" for i, group in enumerate(sorted_groups)}
        temp_df['group'] = temp_df['group'].map(group_map)
        adata_sub.obs['pnCAFs_group'] = temp_df['group']
        group_singleton = pd.concat([temp_df, singletons])
        adata_sub.obs['group_singleton'] = group_singleton
        # Plot 4: Groups after removing singletons
        fig, ax = plt.subplots(figsize=(6,6))
        sc.pl.spatial(adata_sub, color='group_singleton', size=spot_size, bw=False, alpha_img=1, crop_coord=coords, title=f'{section_name} - Groups (with singletons)', ax=ax)
        plt.suptitle(f'{section_name} - Groups and singletons after removal', y=1.05, fontweight='bold', x=0.5)
        # Save the fourth plot to the PDF
        pdf.savefig(fig)
        plt.close(fig)
    out_dir_csv = os.path.join(res_dir, '05_zoom_in_HandE', 'csv_files', 'pnCAFs_score_per_sample', today)
    os.makedirs(out_dir_csv, exist_ok=True)
    to_save_df = adata_sub.obs[['pnCAFs_score', 'group', 'group_singleton', 'pnCAFs_group']]
    to_save_df.columns = ['pnCAFs_score', 'all_groups', 'group_singleton', 'group_singleton_removed']
    filename = f'{section_name}_pnCAFs_groups.csv'
    to_save_df.to_csv(os.path.join(out_dir_csv, filename))
    print(f'Successfully completed for {section_name}')