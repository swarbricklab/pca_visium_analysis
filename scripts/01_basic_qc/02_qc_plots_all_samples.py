#!/usr/bin/env python3

# --- Aim
# Plot basic qc metrics (calculated using scripts in the original dataset repo)
# Generates 1 pdf with all plots as separate pages


# --- Define environment
# env: spatial-brca-preview

# --- Load packages
import argparse
import scanpy as sc
import anndata as ad
from datetime import date
import numpy as np
import pandas as pd

import os

import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib import ticker as mticker

import scanpy.plotting.palettes as sc_palettes
import matplotlib.patches as mpatches

# Import custom functions
from sophiesfunctions.auto_crop import auto_crop_scanpy  # import my custom functions that are packaged locally as sophiesfunctions
from sophiesfunctions.save_multi_image import save_multi_image

sc.logging.print_header()
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
project_dir = args["project_dir"]
data_dir = args["data_dir"]
h5ad_dir = args["h5ad_dir"]
res_dir = args["res_dir"]
sample_list = args["samples"]
use_data = args['use_data']

print(f'Samples in python: {sample_list}')

# -- Load data
# Multiple samples as a dictionary seems like an okay idea
adata_dict = {}

for sample in sample_list:
    adata_dict[sample] = ad.read(os.path.join(h5ad_dir, 'raw', f'{sample}_raw.h5ad'))

# Add some info for plotting
# And make a sample patient LUT
sample_donor_lut = {}

for sample in sample_list:
    adata_dict[sample].obs['assay_id'] = adata_dict[sample].uns['metadata']['sample_id']
    adata_dict[sample].obs['donor_id'] = adata_dict[sample].uns['metadata']['patient_id']
    sample_donor_lut[sample] = adata_dict[sample].uns['metadata']['patient_id']
    adata_dict[sample].obs['sample_type'] = adata_dict[sample].uns['metadata']['sample_type']
    adata_dict[sample].obs['subtype'] = 'PCa'

# --- Get some info from the sample sheet

# Extract the section_id
sample_dir = os.path.join(project_dir, 'config')
sample_sheet_file = 'sample_sheet.csv'

sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))
sample_lut = sample_sheet.set_index('sample_id')


# --- Order based on unique sample ID

ordered_keys = adata_dict.keys()
ordered_keys = sorted(ordered_keys)


# --- Extract adjacent benign samples - want to highlight these in the plots

highlight_samples = []

for sample in ordered_keys:
    if sample == '20333':
        sample = '20033'
    if sample_lut.loc[sample, 'type'] == 'adj_benign':
        highlight_samples.append(sample)


# --- Calculate log total counts for each sample (for plotting later)
for sample in adata_dict.keys():
    total_counts_copy = adata_dict[sample].obs['total_counts'] 
    total_counts_copy[total_counts_copy == 0] = 1e-15 # Deal with the 0s
    adata_dict[sample].obs['log10_counts'] = np.log10(total_counts_copy)

# --- Convert to adata list & concatenate (easier for plotting)

adata_list = list(adata_dict.values())

all_adata = ad.concat(adata_list, merge='same')

# --- Plotting

# Set seaborn style to have a grid
sns.set_style("ticks", {'axes.grid' : True})

# Split the adatas by subtype for plotting
all_subtypes = all_adata.obs.subtype.unique()
adata_per_subtype = {}

for subtype in all_subtypes:
    adata_per_subtype[subtype] = all_adata[all_adata.obs.subtype == subtype]

ncols = len(all_subtypes)

# --- Plot 1: violin plot with counts (split by subtype)

sc.set_figure_params(scanpy=True, fontsize=12)


# Calculate the median
median_counts = {}

for sample in adata_dict.keys():
    median_counts[sample] = adata_dict[sample].obs['total_counts'].median()


median_counts_subtype = {}
for subtype in all_subtypes:
    medians = []
    samples = adata_per_subtype[subtype].obs['assay_id'].astype('category')
    for sample in samples.cat.categories:
        medians.append(median_counts[sample])
    median_counts_subtype[subtype] = medians

# Calculate the width of each subtype subplot based on how many samples per subtype

grid_spec = []

for subtype in all_subtypes:
    grid_spec.append(len(adata_per_subtype[subtype].obs['assay_id'].unique()))

# Get some colors for plotting - color by patient id

all_patients = []

for subtype in all_subtypes:
    tmp = pd.Categorical(adata_per_subtype[subtype].obs['donor_id'])
    all_patients += tmp.categories.tolist()

no_patients = len(all_patients)

colors = sc_palettes.default_20[:no_patients] # TODO: fix this so it works for more than 20 colors

patient_colors = {}

for n, patient in enumerate(all_patients):
    patient_colors[patient] = colors[n]

# Match the samples to their color

plotting_colors = {}

for n, subtype in enumerate(all_subtypes):
    samples = adata_per_subtype[subtype].obs['assay_id'].astype('category')
    tmp_colors = []
    for sample in samples.cat.categories:
        key = sample_donor_lut[sample]
        if key in patient_colors:
            tmp_colors += [patient_colors[key]]
    plotting_colors[subtype] = tmp_colors


# Make the actual plot

fig, axs = plt.subplots(ncols=ncols, figsize=(5*ncols, 5), sharey=True, gridspec_kw={'width_ratios': grid_spec})

if ncols == 1:
    axs = np.reshape(axs, 1)

for n, subtype in enumerate(all_subtypes):
    sc.pl.violin(adata_per_subtype[subtype], 
                keys = ['total_counts'],
                groupby = 'assay_id', stripplot=False, ax = axs[n], title = subtype, linewidth=1, palette = plotting_colors[subtype])
    for collection in axs[n].collections: # Adjust the colors of the violin
        # Get the RGBA values of the fill color
        fill_color = collection.get_facecolor()
        darkened_color = fill_color.copy()
        darkened_color[:,:-1] *= 0.8   # multiply RGB values by 0.8 to darken
        collection.set_edgecolor(darkened_color) 
        fill_color[0][-1] = 0.8 # make the fill color a bit transparent
        collection.set_facecolor(fill_color)
    # make the subtype the title of each subplot
    axs[n].set_title('Counts per spot', fontweight='bold')
    # remove the x-axis label
    axs[n].set_xlabel('')
    # move x-axis labels to the left
    xlabels = axs[n].get_xticklabels()
    axs[n].set_xticklabels(xlabels, rotation=45, ha='right')
    # plot median dots
    medians = median_counts_subtype[subtype]
    pos = range(len(medians))
    for tick, median in zip(pos, medians):
        axs[n].scatter(tick, median, marker='D', color='black', s=10, zorder=3)
    if n > 0:
    #     # remove y-axis
    #     axs[n].set_yticks([])
    #     axs[n].set_yticklabels([])
        axs[n].set_ylabel('')
    #     # axs[n].spines['left'].set_visible(False)

xmin = axs[n].get_xlim()[0]
xmax = axs[n].get_xlim()[1]
for sample in highlight_samples:
    index = ordered_keys.index(sample)
    # Find the corresponding x-tick position
    # Calculate relative positions based on the x-tick position
    xmin_relative = (index - 0.5) / (len(ordered_keys) - 1)
    xmax_relative = (index + 0.5) / (len(ordered_keys) - 1)
    axs[n].axhline(y=0, xmin=xmin_relative, xmax=xmax_relative, color='black', linewidth=2)

# Get all x-tick positions, including minor ticks
xtick_positions = axs[n].get_xaxis().get_ticklocs()

for sample in highlight_samples:
    index = ordered_keys.index(sample)
    x_values = [index-0.5, index+0.5]
    y_values = [axs[n].get_ylim()[0], axs[n].get_ylim()[0]]
    axs[0].axline(x_values, y_values, linewidth=5, color='black')






    # Find the corresponding x-tick position
    xtick_position = xtick_positions[index]
    # Calculate relative positions based on the x-tick position
    xmin_relative = (xtick_position - 0.5) / (len(xtick_positions) - 1)
    xmax_relative = (xtick_position + 0.5) / (len(xtick_positions) - 1)
    axs[n].axhline(y=0, xmin=xmin_relative, xmax=xmax_relative, color='gray', linewidth=2)


xtick_positions = axs[0].get_xticks()
for sample in highlight_samples:
    xmin = axs[n].get_xlim()[0]
    xmax = axs[n].get_xlim()[1]
    index = ordered_keys.index(sample)
    # Find the corresponding x-tick position
    xtick_position = xtick_positions[index]
    index = ordered_keys.index(sample)
    axs[n].axhline(y=axs[n].get_ylim()[0], xmin=index - 0.5, xmax=index + 0.5, color='red', linewidth=5)


axs[n].axhline(axs[n].get_ylim()[0], xmin=0, xmax=(index + 0.5)/10, color='blue', linewidth=5)

axs[n].axhline(axs[n].get_ylim()[0], xmin=0, xmax=0.1, color='blue', linewidth=5)


axs[n].get_ylim()[0]

axs[n].get_xlim()


plt.suptitle(f'Number of counts per spot', fontweight='bold', y=1.02, x=0.5, ha='center')

plt.savefig('test.pdf', bbox_inches='tight')

# --- Plot 2: violin plot log(total counts), split by subtype

# Calculate the median

median_log_counts = median_counts_subtype.copy()

for subtype in all_subtypes:
    for i in range(0, len(median_log_counts[subtype])):
        median_log_counts[subtype][i] = np.log10(median_log_counts[subtype][i])

fig, axs = plt.subplots(ncols=ncols, figsize=(5*ncols, 5), sharey=True, gridspec_kw={'width_ratios': grid_spec})

if ncols == 1:
    axs = np.reshape(axs, 1)

for n, subtype in enumerate(all_subtypes):
    sc.pl.violin(adata_per_subtype[subtype], 
                keys = ['log10_counts'],
                groupby = 'assay_id', stripplot=False, ax = axs[n], title = subtype, linewidth=1, palette = plotting_colors[subtype])
    for collection in axs[n].collections: # Adjust the colors of the violin
        # Get the RGBA values of the fill color
        fill_color = collection.get_facecolor()
        darkened_color = fill_color.copy()
        darkened_color[:,:-1] *= 0.8   # multiply RGB values by 0.8 to darken
        collection.set_edgecolor(darkened_color) 
        fill_color[0][-1] = 0.8 # make the fill color a bit transparent
        collection.set_facecolor(fill_color)
    # make the subtype the title of each subplot
    axs[n].set_title('Counts per spot', fontweight='bold')
    # remove the x-axis label
    axs[n].set_xlabel('')
    # move x-axis labels to the left
    xlabels = axs[n].get_xticklabels()
    axs[n].set_xticklabels(xlabels, rotation=45, ha='right')
    # set up y-axis log scale layout
    axs[n].set_ylim(bottom=1)
    ymin, ymax = axs[n].get_ylim()
    ymin = 0
    tick_range = np.arange(np.floor(ymin), ymax)
    axs[n].yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)], minor=True)
    axs[n].yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
    # axs[n].set_yticks(y_tick_labels, y_ticks)
    # plot median dots
    medians = median_log_counts[subtype]
    pos = range(len(medians))
    for tick, median in zip(pos, medians):
        axs[n].scatter(tick, median, marker='D', color='black', s=10, zorder=3)
    if n > 0:
    #     # remove y-axis
    #     axs[n].set_yticks([])
    #     axs[n].set_yticklabels([])
        axs[n].set_ylabel('')
    #     # axs[n].spines['left'].set_visible(False)

plt.suptitle(f'Log(number of counts) per spot', fontweight='bold', y=1.02, x=0.5, ha='center')

plt.savefig('test.pdf', bbox_inches='tight')



# --- Plot 2: violin plot with genes per spot, split by subtype

# Calculate the median

median_genes = {}

for sample in adata_dict.keys():
    median_genes[sample] = adata_dict[sample].obs['n_genes_by_counts'].median()


median_genes_subtype = {}
for subtype in all_subtypes:
    medians = []
    samples = adata_per_subtype[subtype].obs['assay_id'].astype('category')
    for sample in samples.cat.categories:
        medians.append(median_genes[sample])
    median_genes_subtype[subtype] = medians

# Make the actual plot

fig, axs = plt.subplots(ncols=ncols, figsize=(5*ncols, 5), sharey=True, gridspec_kw={'width_ratios': grid_spec})

if ncols == 1:
    axs = np.reshape(axs, 1)

for n, subtype in enumerate(all_subtypes):
    sc.pl.violin(adata_per_subtype[subtype], 
                keys = ['n_genes_by_counts'],
                groupby = 'assay_id', stripplot=False, ax = axs[n], title = subtype, linewidth=1, palette=plotting_colors[subtype])
    for collection in axs[n].collections: # Adjust the colors of the violin
        # Get the RGBA values of the fill color
        fill_color = collection.get_facecolor()
        darkened_color = fill_color.copy()
        darkened_color[:,:-1] *= 0.8   # multiply RGB values by 0.8 to darken
        collection.set_edgecolor(darkened_color) 
        fill_color[0][-1] = 0.8 # make the fill color a bit transparent
        collection.set_facecolor(fill_color)
    axs[n].set_title('Genes per spot', fontweight='bold')
    # move x-axis tick labels to the left
    xlabels = axs[n].get_xticklabels()
    axs[n].set_xticklabels(xlabels, rotation=45, ha='right')
    # remove the x-axis label
    axs[n].set_xlabel('')
    # plot median dots
    medians = median_genes_subtype[subtype]
    pos = range(len(medians))
    for tick, median in zip(pos, medians):
        axs[n].scatter(tick, median, marker='D', color='black', s=10, zorder=3)
    if n > 0:
    #     # remove y-axis
    #     axs[n].set_yticks([])
    #     axs[n].set_yticklabels([])
        axs[n].set_ylabel('')
    #     # axs[n].spines['left'].set_visible(False)

plt.suptitle(f'Number of genes per spot', fontweight='bold', y=1.02, x=0.5, ha='center')

plt.savefig('test.pdf', bbox_inches='tight')


# --- Plot 3: violin plot with pct_mt per spot, split by subtype

# Calculate the median

median_mt = {}

for sample in adata_dict.keys():
    median_mt[sample] = adata_dict[sample].obs['pct_counts_mt'].median()

median_mt_subtype = {}
for subtype in all_subtypes:
    medians = []
    samples = adata_per_subtype[subtype].obs['assay_id'].astype('category')
    for sample in samples.cat.categories:
        medians.append(median_mt[sample])
    median_mt_subtype[subtype] = medians

# Make the actual plot

fig, axs = plt.subplots(ncols=ncols, figsize=(5*ncols, 5), sharey=True, gridspec_kw={'width_ratios': grid_spec})

if ncols == 1:
    axs = np.reshape(axs, 1)

for n, subtype in enumerate(all_subtypes):
    sc.pl.violin(adata_per_subtype[subtype], 
                keys = ['pct_counts_mt'],
                groupby = 'assay_id', stripplot=False, ax = axs[n], title = subtype, linewidth=1, palette=plotting_colors[subtype])
    for collection in axs[n].collections: # Adjust the colors of the violin
        # Get the RGBA values of the fill color
        fill_color = collection.get_facecolor()
        darkened_color = fill_color.copy()
        darkened_color[:,:-1] *= 0.8   # multiply RGB values by 0.8 to darken
        collection.set_edgecolor(darkened_color) 
        fill_color[0][-1] = 0.8 # make the fill color a bit transparent
        collection.set_facecolor(fill_color)
    # make the subtype the title of each subplot
    axs[n].set_title('% MT counts', fontweight='bold')
    # move x-axis tick labels to the left
    xlabels = axs[n].get_xticklabels()
    axs[n].set_xticklabels(xlabels, rotation=45, ha='right')
    # remove the x-axis label
    axs[n].set_xlabel('')
    # plot median dots
    medians = median_mt_subtype[subtype]
    pos = range(len(medians))
    for tick, median in zip(pos, medians):
        axs[n].scatter(tick, median, marker='D', color='black', s=10, zorder=3)
    if n > 0:
    #     # remove y-axis
    #     axs[n].set_yticks([])
    #     axs[n].set_yticklabels([])
          axs[n].set_ylabel('')
    #     # axs[n].spines['left'].set_visible(False)

plt.suptitle(f'Percentage of mitochondrial counts per spot', fontweight='bold', y=1.02, x=0.5, ha='center')

plt.savefig('test.pdf', bbox_inches='tight')


# --- Plot 4: Plot the number of spots per gene, split by subtype

# Calculate the median

median_spots_genes = {}

for sample in adata_dict.keys():
    median_spots_genes[sample] = adata_dict[sample].var['n_cells_by_counts'].median()


median_spots_genes_subtype = {}
for subtype in all_subtypes:
    medians = []
    samples = adata_per_subtype[subtype].obs['assay_id'].astype('category')
    for sample in samples.cat.categories:
        medians.append(median_spots_genes[sample])
    median_spots_genes_subtype[subtype] = medians

# Extract a dataframe with n_cells_by_counts for every sample, concat per subtype

n_cells_by_counts_dict = {}

for sample in adata_dict:
    # subtype = adata_dict[sample].uns['metadata']['subtype']
    subtype = 'PCa'
    tmp = adata_dict[sample].var['n_cells_by_counts']
    add_sample_col = pd.Series(adata_dict[sample].uns['metadata']['sample_id'], name = 'sample').repeat(len(tmp))
    tmp_df = pd.concat([tmp.reset_index(drop=True), add_sample_col.reset_index(drop=True)], axis=1)
    tmp_df['sample'] = tmp_df['sample'].astype('category')
    if not subtype in n_cells_by_counts_dict.keys():
        n_cells_by_counts_dict[subtype] = tmp_df
    else:
        n_cells_by_counts_dict[subtype] = pd.concat([n_cells_by_counts_dict[subtype].reset_index(drop=True), tmp_df.reset_index(drop=True)], axis=0)


# Make the actual plot

fig, axs = plt.subplots(ncols=ncols, figsize=(5*ncols, 5), sharey=True, gridspec_kw={'width_ratios': grid_spec})

if ncols == 1:
    axs = np.reshape(axs, 1)

for n, subtype in enumerate(all_subtypes):
    sns.violinplot(data = n_cells_by_counts_dict[subtype], 
                x = 'sample', y = 'n_cells_by_counts', stripplot=False, ax = axs[n], title = subtype, linewidth=1, palette=plotting_colors[subtype], 
                inner = None, order=sorted(n_cells_by_counts_dict[subtype]['sample'].unique()), scale='width')
    for collection in axs[n].collections: # Adjust the colors of the violin
        # Get the RGBA values of the fill color
        fill_color = collection.get_facecolor()
        darkened_color = fill_color.copy()
        darkened_color[:,:-1] *= 0.8   # multiply RGB values by 0.8 to darken
        collection.set_edgecolor(darkened_color) 
        fill_color[0][-1] = 0.8 # make the fill color a bit transparent
        collection.set_facecolor(fill_color)
    axs[n].set_title('Spots per gene', fontweight='bold')
    # move x-axis tick labels to the left
    xlabels = axs[n].get_xticklabels()
    axs[n].set_xticklabels(xlabels, rotation=45, ha='right')
    # remove the x-axis label
    axs[n].set_xlabel('')
    # plot median dots
    medians = median_spots_genes_subtype[subtype]
    pos = range(len(medians))
    for tick, median in zip(pos, medians):
        axs[n].scatter(tick, median, marker='D', color='black', s=10, zorder=3)
    if n > 0:
    #     # remove y-axis
    #     axs[n].set_yticks([])
    #     axs[n].set_yticklabels([])
          axs[n].set_ylabel('')
    #     # axs[n].spines['left'].set_visible(False)

plt.suptitle(f'Number of spots per gene', fontweight='bold', y=1.02, x=0.5, ha='center')
plt.savefig('test.pdf', bbox_inches='tight')


# --- Plot the total spots on tissue, total counts, total genes per sample

barplot_dict = {}

for sample in adata_dict:
    subtype = 'PCa'
    # subtype = adata_dict[sample].uns['metadata']['subtype']
    in_tissue = [len(adata_dict[sample].obs['in_tissue'])]
    n_counts = [sum(adata_dict[sample].obs['total_counts'])]
    median_counts = [adata_dict[sample].obs['total_counts'].median()]
    n_genes = [len(adata_dict[sample].var_names.unique())]
    median_genes = [adata_dict[sample].obs['n_genes_by_counts'].median()]
    median_pct_counts_mt = [adata_dict[sample].obs['pct_counts_mt'].median()]
    add_sample_col = [adata_dict[sample].uns['metadata']['sample_id']]
    pct_low_q = [sum(adata_dict[sample].obs['low_qc'] != 'pass')/adata_dict[sample].obs.shape[0]*100]
    tmp_df = pd.DataFrame({'sample': add_sample_col, 'in_tissue': in_tissue, 'n_counts': n_counts, 'median_counts': median_counts,
                            'n_genes': n_genes, 'median_genes': median_genes, 'median_pct_counts_mt': median_pct_counts_mt, 'pct_low_q': pct_low_q})
    tmp_df['sample'] = tmp_df['sample'].astype('category')
    if not subtype in barplot_dict.keys():
        barplot_dict[subtype] = tmp_df
    else:
        barplot_dict[subtype] = pd.concat([barplot_dict[subtype].reset_index(drop=True), tmp_df.reset_index(drop=True)], axis=0)

# --- Plot 5: Plot the total number of counts (barplot)

fig, axs = plt.subplots(ncols=ncols, figsize=(5*ncols, 5), sharey=True, gridspec_kw={'width_ratios': grid_spec})

if ncols == 1:
    axs = np.reshape(axs, 1)

for n, subtype in enumerate(all_subtypes):
    sns.barplot(data = barplot_dict[subtype], 
                x = 'sample', y = 'n_counts', ax = axs[n], linewidth=1, palette=plotting_colors[subtype], 
                order=sorted(barplot_dict[subtype]['sample'].unique()))
    for bar in axs[n].patches: # Adjust the colors of the barplot
        # Get the RGBA values of the fill color
        fill_color = bar.get_facecolor()
        fill_color = list(fill_color)
        darkened_color = fill_color.copy()
        for i in range(len(darkened_color)-1):
            darkened_color[i] *= 0.8
        fill_color[-1] = 0.8 # make the fill color a bit transparent
        bar.set_color(tuple(fill_color))
        bar.set_edgecolor(tuple(darkened_color))
    # make the subtype the title of each subplot
    axs[n].set_title('Total counts', fontweight='bold')
    # move x-axis tick labels to the left
    xlabels = axs[n].get_xticklabels()
    axs[n].set_xticklabels(xlabels, rotation=45, ha='right')
    # remove the x-axis label
    axs[n].set_xlabel('')
    if n > 0:
    #     # remove y-axis
    #     axs[n].set_yticks([])
    #     axs[n].set_yticklabels([])
          axs[n].set_ylabel('')
    #     # axs[n].spines['left'].set_visible(False)

plt.suptitle(f'Total number of counts', fontweight='bold', y=1.02, x=0.5, ha='center')
plt.savefig('test.pdf', bbox_inches='tight')


# --- Plot 6: Plot the number of spots 'in tissue' (barplot)

fig, axs = plt.subplots(ncols=ncols, figsize=(5*ncols, 5), sharey=True, gridspec_kw={'width_ratios': grid_spec})

if ncols == 1:
    axs = np.reshape(axs, 1)

for n, subtype in enumerate(all_subtypes):
    sns.barplot(data = barplot_dict[subtype], 
                x = 'sample', y = 'in_tissue', ax = axs[n], linewidth=1, palette=plotting_colors[subtype], 
                order=sorted(barplot_dict[subtype]['sample'].unique()))
    for bar in axs[n].patches: # Adjust the colors of the barplot
        # Get the RGBA values of the fill color
        fill_color = bar.get_facecolor()
        fill_color = list(fill_color)
        darkened_color = fill_color.copy()
        for i in range(len(darkened_color)-1):
            darkened_color[i] *= 0.8
        fill_color[-1] = 0.8 # make the fill color a bit transparent
        bar.set_color(tuple(fill_color))
        bar.set_edgecolor(tuple(darkened_color))
    axs[n].set_title('Spots in tissue', fontweight='bold', loc='left') # Move title so it doesn't overlap
    # move x-axis tick labels to the left
    xlabels = axs[n].get_xticklabels()
    axs[n].set_xticklabels(xlabels, rotation=45, ha='right')
    # remove the x-axis label
    axs[n].set_xlabel('')
    if n > 0:
    #     # remove y-axis
    #     axs[n].set_yticks([])
    #     axs[n].set_yticklabels([])
          axs[n].set_ylabel('')
    #     # axs[n].spines['left'].set_visible(False)

plt.suptitle(f"Number of spots 'in tissue'", fontweight='bold', y=1.02, x=0.5, ha='center')


# --- Plot 7: Plot the percentage of low quality spots (barplot)

fig, axs = plt.subplots(ncols=ncols, figsize=(5*ncols, 5), sharey=True, gridspec_kw={'width_ratios': grid_spec})

if ncols == 1:
    axs = np.reshape(axs, 1)

for n, subtype in enumerate(all_subtypes):
    sns.barplot(data = barplot_dict[subtype], 
                x = 'sample', y = 'pct_low_q', ax = axs[n], linewidth=1, palette=plotting_colors[subtype], 
                order=sorted(barplot_dict[subtype]['sample'].unique()))
    for bar in axs[n].patches: # Adjust the colors of the barplot
        # Get the RGBA values of the fill color
        fill_color = bar.get_facecolor()
        fill_color = list(fill_color)
        darkened_color = fill_color.copy()
        for i in range(len(darkened_color)-1):
            darkened_color[i] *= 0.8
        fill_color[-1] = 0.8 # make the fill color a bit transparent
        bar.set_color(tuple(fill_color))
        bar.set_edgecolor(tuple(darkened_color))
    # make the subtype the title of each subplot
    axs[n].set_title('Pct low quality spots', fontweight='bold')
    # move x-axis tick labels to the left
    xlabels = axs[n].get_xticklabels()
    axs[n].set_xticklabels(xlabels, rotation=45, ha='right')
    # remove the x-axis label
    axs[n].set_xlabel('')
    if n > 0:
    #     # remove y-axis
    #     axs[n].set_yticks([])
    #     axs[n].set_yticklabels([])
          axs[n].set_ylabel('')
    #     # axs[n].spines['left'].set_visible(False)

plt.suptitle(f"Percentage of low quality spots", fontweight='bold', y=1.02, x=0.5, ha='center')


# --- Plot 8: Dotplots with the median number of genes vs number of counts, also pct mt vs counts

barplot_df = pd.concat(barplot_dict.values(), ignore_index=True)

# Get plotting colors in the right order
plotting_colors_all = []
for sample in list(barplot_df['sample']):
    patient_id = sample_donor_lut[sample]
    color = patient_colors[patient_id]
    plotting_colors_all.append(color)

# Create a list of handles and labels for the legend
handles = []
labels = []
for cat, color in patient_colors.items():
    handles.append(mpatches.Patch(color=color))
    labels.append(cat)

# Plotting
fig, axs = plt.subplots(ncols=2, figsize=(4*2, 4))

sns.scatterplot(barplot_df, x='median_counts', y='median_genes', hue='sample', palette=plotting_colors_all, ax=axs[0], legend=False, linewidth=0, alpha=0.8)
axs[0].set_title('Genes vs counts per spot', fontweight='bold')

sns.scatterplot(barplot_df, x='median_counts', y='median_pct_counts_mt', hue='sample', palette=plotting_colors_all, ax=axs[1], legend=False, linewidth=0, alpha=0.8)
axs[1].set_title('Percentage MT vs counts per spot', fontweight='bold')

plt.suptitle(f"Median number of genes / percentage mitochondrial counts per spot vs counts", fontweight='bold', y=1.02, x=0.5, ha='center')

plt.tight_layout()

# set the marker for each handle to 'o' (circle)
plt.legend(handles=handles, labels=labels, bbox_to_anchor=(1.05, 1), loc='upper left')

plt.savefig('test.pdf', bbox_inches='tight')

# --- Save output
# Put in a timestamped folder to avoid overwriting older plots
qc_dir = os.path.join(res_dir, 'basic_qc')
os.makedirs(qc_dir, exist_ok=True)

filename = f"All_samples_qc_plots.pdf"

today = date.today()
today = today.strftime("%Y%m%d")

out_dir = os.path.join(qc_dir, 'all_samples', today)
os.makedirs(out_dir, exist_ok=True)

print(f'Saving plot as {filename} in {qc_dir}')
save_multi_image(os.path.join(out_dir, filename))

print("Script succesfully completed")
