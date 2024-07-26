#!/usr/bin/env python3

# --- Aim
# Plot basic qc metrics (calculated using scripts in the original dataset repo)
# Generates 1 pdf with all plots as separate pages


# --- Define environment
# env: pca-visium

# --- Load packages
import argparse
import scanpy as sc
import anndata as ad
from datetime import date
import numpy as np
import pandas as pd

import matplotlib.patheffects as pe

import os

import matplotlib.pyplot as plt
import seaborn as sns

# Import custom functions
from sophiesfunctions.auto_crop import auto_crop_scanpy  # import my custom functions that are packaged locally as sophiesfunctions
from sophiesfunctions.save_multi_image import save_multi_image

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & do basic QC",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
# parser.add_argument("--raw_data_dir", help="Path to raw data")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad input/output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--min_count", help="Minimal number of counts per spot", type=int)
parser.add_argument("--min_gene", help="Minimal number of genes per spot", type=int)
parser.add_argument("--min_spots", help="Minimal spots per gene", type=int)

# parser.add_argument("--outs_path", help="Path to spaceranger 'outs'")
# parser.add_argument("--use_data", help="Which data to use (log_norm or SCT)")

args = vars(parser.parse_args())

# --- Arguments
# raw_data_dir = args["raw_data_dir"]
project_dir=args["project_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]

min_count=args["min_count"]
min_gene=args["min_gene"]
min_spots = args["min_spots"]


# Print the arguments
print(f'.h5ad dir: {h5ad_dir}')
print(f'Sample id: {sample_id}')

print(f'Min count: {min_count}')
print(f'Min gene: {min_gene}')


# --- Get some info from the sample sheet

# Extract the section_id
sample_dir = os.path.join(project_dir, 'config')
sample_sheet_file = 'sample_sheet_mengxiao.csv'

sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))
sample_lut = sample_sheet.set_index('sample_id')

sample_name = sample_lut.loc[sample_id, 'sample_name']
sample_type = sample_lut.loc[sample_id, 'type']


# --- Load the adata object

adata = ad.read(os.path.join(h5ad_dir, 'filtered', f'{sample_name}_filtered.h5ad'))


# --- Plot QC metrics

# Enable gridlines behind the plot for all plots
plt.rcParams['axes.axisbelow'] = True

# Calculate cropping coords for optimal presentation
coords = auto_crop_scanpy(adata, border_size=0.1)

# --- Plot 1: Plot UMIs / genes per spot violin
median_counts = adata.obs['total_counts'].median().astype(int)
median_gene = adata.obs['n_genes_by_counts'].median().astype(int)

nrows = 1
ncols = 3

fig, axs = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows))

# -- Counts per spot
n = 0
ax = axs[n]
sc.pl.violin(adata, ['total_counts'], jitter=0.2, show=False, size=1.5, height=3, ax=ax)
ax.set_title('Counts per spot', fontweight='bold')

# Add median line
ax.axhline(median_counts, color="lime")
x_coord = ax.get_xlim()[1] - 0.01 * (ax.get_xlim()[1] - ax.get_xlim()[0])
y_buffer = 0.01 * (ax.get_ylim()[1] - ax.get_ylim()[0])
ax.text(x_coord, median_counts + y_buffer, median_counts, fontsize=10, ha='right', va='bottom', path_effects=[pe.withStroke(linewidth=3, foreground='white')])

# Add proposed filter threshold
ax.axhline(min_count, color="red")
ax.text(x_coord, min_count + y_buffer, min_count, fontsize=10, ha='right', va='bottom', path_effects=[pe.withStroke(linewidth=3, foreground='white')])

# -- Genes per spot
n+=1
ax = axs[n]
sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.2, show=False, size=1.5, ax=ax)
ax.axhline(min_gene, color="red")
ax.axhline(median_gene, color="lime")
ax.set_title('Genes per spot', fontweight='bold')

# Add median line
x_coord = ax.get_xlim()[1] - 0.01 * (ax.get_xlim()[1] - ax.get_xlim()[0])
y_buffer = 0.01 * (ax.get_ylim()[1] - ax.get_ylim()[0])
ax.text(x_coord, median_gene + y_buffer, median_gene, fontsize=10, ha='right', va='bottom', path_effects=[pe.withStroke(linewidth=3, foreground='white')])

# Add proposed filter threshold
ax.axhline(min_gene, color="red")
ax.text(x_coord, min_gene + y_buffer, min_gene, fontsize=10, ha='right', va='bottom', path_effects=[pe.withStroke(linewidth=3, foreground='white')])


# -- Percentage mitochondrial counts
n+=1
ax = axs[n]
sc.pl.violin(adata, ['pct_counts_mt'], jitter=0.2, show=False, size=1.5, ax=ax)
ax.set_title('% mitochondrial counts', fontweight='bold')

# -- Beautify the plots
for ax in axs:
    for collection in ax.collections: # Adjust the colors of the violin
            # Get the RGBA values of the fill color
            fill_color = collection.get_facecolor()
            darkened_color = fill_color.copy()
            darkened_color[:,:-1] *= 0.8   # multiply RGB values by 0.8 to darken
            collection.set_edgecolor(darkened_color) 
            fill_color[0][-1] = 0.8 # make the fill color a bit transparent
            collection.set_facecolor(fill_color)

plt.suptitle(f'{sample_id} ({sample_type}) - QC metrics (green line = median, red line = potential filtering threshold)', y=1, fontweight='bold', x=0.5)
plt.tight_layout()


# --- Plot 2: Plot all of the above in space as well

sc.set_figure_params(scanpy=True, fontsize=20)

to_plot = [None, 'total_counts', 'n_genes_by_counts', 'pct_counts_mt']
spot_size=1.5

ncols = len(to_plot)

nrows = 1

plt.figure(figsize=(8*ncols, 6*nrows))

for n, plot in enumerate(to_plot):
    # add a new subplot iteratively
    ax = plt.subplot(nrows, ncols, n + 1)
    if n==0:
        sc.pl.spatial(adata, color=None, crop_coord=coords, size=1.2, alpha_img=1, title=f'{sample_id}', ax=ax)
    else:
        sc.pl.spatial(adata, color=[plot], crop_coord=coords, size=spot_size, alpha_img=0.5, ax=ax, cmap='magma_r')

plt.suptitle(f'{sample_id} ({sample_type}) - QC metrics spatial', y=1.1, fontweight='bold', x=0.5)
plt.subplots_adjust(wspace=0)


# --- Plot 3: Plot UMIs / genes per spot histogram

sc.set_figure_params(scanpy=True, fontsize=15)

fig, axs = plt.subplots(2, 4, figsize=(27, 10), gridspec_kw={'wspace': 0.3, 'hspace': 0.3})
axs = axs.ravel()

# -- Counts histograms

n=0
ax = axs[n]
sns.histplot(adata.obs['total_counts'], kde=False, ax=ax)
ax.set_title('Total counts', fontweight='bold')
ax.axvline(median_counts, color="lime")
ax.text(x=median_counts+0.02*ax.get_xlim()[1], y=0.75*ax.get_ylim()[1], s=f"Median: {median_counts}",
        path_effects=[pe.withStroke(linewidth=3, foreground='white')])

n+=1
ax = axs[n]
sns.histplot(adata.obs['total_counts'][adata.obs['total_counts'] < 3000],
             kde=False, bins=60, ax=ax)
ax.set_title('Total counts < 3000', fontweight='bold')
ax.axvline(min_count, color="red")
if median_counts < ax.get_xlim()[1]:
    ax.axvline(median_counts, color="lime")

n+=1
ax = axs[n]
sns.histplot(adata.obs['total_counts'][adata.obs['total_counts'] < 20000],
             kde=False, bins=60, ax=ax)
ax.set_title('Total counts < 20,000', fontweight='bold')
ax.axvline(min_count, color="red")
if median_counts < ax.get_xlim()[1]:
    ax.axvline(median_counts, color="lime")

n+=1
ax = axs[n]
sns.histplot(adata.obs['total_counts'][adata.obs['total_counts'] > 25000],
             kde=False, bins=60, ax=ax)
ax.set_title('Total counts > 25,000', fontweight='bold')

# -- Genes histograms

n+=1
ax = axs[n]
sns.histplot(adata.obs['n_genes_by_counts'], kde=False, bins=60,
             ax=ax).set(title='Total genes')
ax.set_title('Total genes', fontweight='bold')
ax.axvline(median_gene, color="lime")
ax.text(x=median_gene+0.02*ax.get_xlim()[1], y=0.75*ax.get_ylim()[1], s=f"Median: {median_gene}",
        path_effects=[pe.withStroke(linewidth=3, foreground='white')])

n+=1
ax = axs[n]
sns.histplot(adata.obs['n_genes_by_counts'][adata.obs['n_genes_by_counts'] < 2000],
             kde=False, bins=60, ax=ax)
ax.set_title('Total genes < 2000', fontweight='bold')
ax.axvline(min_gene, color="red")
if median_gene < ax.get_xlim()[1]:
    ax.axvline(median_gene, color="lime")

n+=1
ax = axs[n]
sns.histplot(adata.obs['n_genes_by_counts'][adata.obs['n_genes_by_counts'] < 4000],
             kde=False, bins=60, ax=ax)
ax.set_title('Total genes < 4000', fontweight='bold')
ax.axvline(min_gene, color="red")
if median_gene < ax.get_xlim()[1]:
    ax.axvline(median_gene, color="lime")

n+=1
ax = axs[n]
sns.histplot(adata.obs['n_genes_by_counts'][adata.obs['n_genes_by_counts'] > 6000],
             kde=False, bins=60, ax=ax)
ax.set_title('Total genes > 6000', fontweight='bold')

plt.suptitle(f'{sample_id} ({sample_type}) - QC metrics histograms: total counts & total genes per spot', y=1, fontweight='bold', x=0.5)


# --- Plot 4: Number of genes & UMI correlations

nrows = 1
ncols = 2

# Mark spots below the filter thresholds
adata.obs['low_qc'] = np.where((adata.obs['total_counts'] <= min_count) | (adata.obs['n_genes_by_counts'] <= min_gene), 'low counts/genes', 'pass')

# Need to do this to make sure the pass isn't colored red if there is none below pass
adata.obs['low_qc'] = pd.Categorical(adata.obs['low_qc'], categories=['pass', 'low counts/genes'], ordered=True)


## Check n_gene, mt and UMI correlation
fig, axs = plt.subplots(nrows, ncols, figsize=(7*ncols, 5*nrows))

n=0
ax = axs[n]
# n_gene
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='low_qc', palette=['lightgray', 'red'], size=20, ax=ax)
ax.set_title('Genes vs counts', fontweight='bold')
# Rotate xlabels
xlabels = ax.get_xticklabels()
ax.set_xticklabels(xlabels, rotation=45, ha='right')
ax.axvline(x=min_count, linestyle='--', linewidth=1, alpha=0.5, color='black')
ax.axhline(y=min_gene, linestyle='--', linewidth=1, alpha=0.5, color='black')

# mt
n+=1
ax = axs[n]
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', color='low_qc', palette=['lightgray', 'red'], size=20, ax=axs[1])
ax.set_title('pct_mt vs counts', fontweight='bold')
# Rotate xlabels
xlabels = ax.get_xticklabels()
ax.set_xticklabels(xlabels, rotation=45, ha='right')
plt.suptitle(f'{sample_id} ({sample_type}) - counts correlations', y=1, fontweight='bold', x=0.5)

plt.tight_layout()


# --- Plot 5: Histogram of the number of spots each gene is expressed in

sc.set_figure_params(scanpy=True, fontsize=12)

ncols = 3

fig, axs = plt.subplots(nrows, ncols, figsize=(5*ncols, 3*nrows), gridspec_kw={'wspace': 0.3})

median_spots = adata.var['n_cells_by_counts'].median().astype(int)

n=0
ax = axs[n]
sns.histplot(adata.var['n_cells_by_counts'], kde=False, bins=60, ax=ax)
ax.set_title('Number of spots per gene', fontweight='bold')
ax.axvline(median_spots, color="lime")
ax.text(x=median_spots+0.02*ax.get_xlim()[1], y=0.75*ax.get_ylim()[1], s=f"Median: {median_spots}")

n+=1
ax = axs[n]
sns.histplot(adata.var['n_cells_by_counts'][adata.var['n_cells_by_counts'] < 200], kde=False, bins=60, ax=ax)
ax.set_title('<200 spots', fontweight='bold')
if median_spots < ax.get_xlim()[1]:
    ax.axvline(median_spots, color="lime")

n+=1
ax = axs[n]
sns.histplot(adata.var['n_cells_by_counts'][adata.var['n_cells_by_counts'] < 25], kde=False, bins=25, ax=ax)
ax.set_title('<25 spots', fontweight='bold')
ax.axvline(min_spots, color="red")
ax.text(x=min_spots+0.02*ax.get_xlim()[1], y=0.75*ax.get_ylim()[1], s=f"Min spots: {min_spots}")


plt.suptitle(f'{sample_id} ({sample_type}) - Number of spots each gene is detected in', y=1.15, fontweight='bold', x=0.5)


# --- Plot 6: Spatial plot of low quality spots

sc.set_figure_params(scanpy=True, fontsize=20)

nrows = 1
ncols = 3

# Flag spots with low counts
adata.obs['low_counts'] = np.where(adata.obs['total_counts'] <= min_count, 'low_counts', 'pass')
adata.obs['low_counts'] = pd.Categorical(adata.obs['low_counts'], categories=['pass', 'low_counts'], ordered=True)

# Flag spots with low genes
adata.obs['low_gene'] = np.where(adata.obs['n_genes_by_counts'] <= min_gene, 'low_gene', 'pass')
adata.obs['low_gene'] = pd.Categorical(adata.obs['low_gene'], categories=['pass', 'low_gene'], ordered=True)

# Check if there's a cytassist image
sample_name = list(adata.uns['spatial'].keys())[0]
if 'cytassist' in adata.uns['spatial'][sample_name]['images'].keys():
    ncols = 2
    nrows = 2

# Plotting
fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8*ncols, 8*nrows))
axs = axs.ravel()

# Plot low counts
n=0
ax = axs[n]

sc.pl.spatial(adata, color=["low_counts"], alpha_img = 0.7, palette = ['lightgrey', 'red'], crop_coord=coords, ax=ax)
ax.set_title(ax.title.get_text(), fontweight='bold')

# Plot low genes
n+=1
ax = axs[n]

sc.pl.spatial(adata, color=["low_gene"], alpha_img = 0.7, palette = ['lightgrey', 'red'], crop_coord=coords, ax=ax)
ax.set_title(ax.title.get_text(), fontweight='bold')

# Also plot the H&E vs cytassist image
n+=1
ax = axs[n]

sc.pl.spatial(adata, alpha_img = 1, crop_coord=coords, ax=ax)
ax.set_title(sample_id, fontweight='bold')

# Check if there is a cytassist image

if 'cytassist' in adata.uns['spatial'][sample_name]['images'].keys():
    n+=1
    ax = axs[n]
    ax.imshow(adata.uns['spatial'][sample_name]['images']['cytassist'])
    plt.grid(False)
    plt.axis('off')
    ax.set_title('Cytassist image', fontweight='bold')

plt.suptitle(f'{sample_id} ({sample_type}) - Flag low QC spots spatially', y=1, fontweight='bold', x=0.5)
plt.tight_layout()


# --- Plot 7: Box plot with top genes in the sample

sc.set_figure_params(scanpy=True, fontsize=12)

sc.pl.highest_expr_genes(adata, n_top=20)

plt.suptitle(f'{sample_id} ({sample_type}) - Top genes', y=1, fontweight='bold', x=0.5)


# --- Save output
# Put in a timestamped folder to avoid overwriting older plots
qc_dir = os.path.join('/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/', 'basic_qc')
os.makedirs(qc_dir, exist_ok=True)

filename = f"{sample_id}_qc_plots.pdf"

today = date.today()
today = today.strftime("%Y%m%d")

out_dir = os.path.join(qc_dir, 'per_sample', today)
os.makedirs(out_dir, exist_ok=True)

print(f'Saving plot as {filename} in {qc_dir}')
save_multi_image(os.path.join(out_dir, filename))

print("Script succesfully completed")