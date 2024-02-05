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

# --- Get some info from the sample sheet

# Extract the section_id
sample_dir = os.path.join(project_dir, 'config')
sample_sheet_file = 'sample_sheet.csv'

sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))
sample_lut = sample_sheet.set_index('sample_id')


# -- Load data
# Multiple samples as a dictionary seems like an okay idea
adata_dict = {}

for sample in sample_list:
    sample_name = sample_lut.loc[sample, 'sample_name']
    adata_dict[sample] = ad.read(os.path.join(h5ad_dir, 'raw', f'{sample_name}_raw.h5ad'))

# Add some info for plotting
# And make a sample patient LUT
sample_donor_lut = {}

for sample in sample_list:
    adata_dict[sample].obs['assay_id'] = sample
    adata_dict[sample].obs['donor_id'] = sample_lut.loc[sample, 'patient_id']
    sample_donor_lut[sample] = sample_lut.loc[sample, 'patient_id']
    adata_dict[sample].obs['sample_type'] = sample_lut.loc[sample, 'type']
    adata_dict[sample].obs['subtype'] = 'PCa'


# --- Order based on unique sample ID

ordered_keys = adata_dict.keys()
ordered_keys = sorted(ordered_keys)


# --- Extract adjacent benign samples - want to highlight these in the plots

highlight_samples = []

for sample in ordered_keys:
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

# Get the positions of samples that should be highlighted
highlight_samples_indices = []
for sample in highlight_samples:
            index = ordered_keys.index(sample)
            highlight_samples_indices.append(index)


# Define a function for plotting the custom violin
# --- Step 1: prepare for plotting the violin: calculate the medians, gridspec and prepare the figure
def prepare_custom_violin(adata_per_subtype, plot_key, ncols=ncols, groupby='assay_id', all_subtypes=all_subtypes, obs=True):
    # Calculate the median (to be able to plot the median on the violin)
    median_counts = {}
    # Deal with variables stored in .obs vs .var
    if obs:
        for sample in adata_dict.keys():
            median_counts[sample] = adata_dict[sample].obs[plot_key].median()
    else:
        for sample in adata_dict.keys():
            median_counts[sample] = adata_dict[sample].var[plot_key].median()
    median_counts_subtype = {}
    for subtype in all_subtypes:
        medians = []
        samples = adata_per_subtype[subtype].obs[groupby].astype('category')
        for sample in samples.cat.categories:
            medians.append(median_counts[sample])
        median_counts_subtype[subtype] = medians
    # Calculate the width of each subtype subplot based on how many samples per subtype
    grid_spec = []
    for subtype in all_subtypes:
        grid_spec.append(len(adata_per_subtype[subtype].obs[groupby].unique()))
    fig, axs = plt.subplots(ncols=ncols, figsize=(5*ncols, 5), sharey=True, gridspec_kw={'width_ratios': grid_spec})
    # If ncols = 1: need to reshape to make axs iterable
    if ncols == 1:
        axs = np.reshape(axs, 1)
    return median_counts_subtype, axs

# --- Step 2: Plot the actual violin
def plot_custom_violin(adata_per_subtype, plot_key, median_counts_subtype, axs, groupby='assay_id', title=None, all_subtypes=all_subtypes):
    # Do the actual plotting
    for n, (subtype, ax) in enumerate(zip(all_subtypes, axs)):
        sc.pl.violin(adata_per_subtype[subtype], keys=[plot_key], groupby=groupby, stripplot=False, ax=ax, title=subtype, linewidth=1, palette=plotting_colors[subtype])
        for i, collection in enumerate(ax.collections):  # Adjust the colors of the violin
            fill_color = collection.get_facecolor()
            darkened_color = fill_color.copy()
            darkened_color[:, :-1] *= 0.8  # multiply RGB values by 0.8 to darken
            collection.set_edgecolor(darkened_color)
            fill_color[0][-1] = 0.8  # make the fill color a bit transparent
            collection.set_facecolor(fill_color)
            if i in highlight_samples_indices:  # If it's an adj_benign sample: mark this sample with slightly lighter color & dashed outline
                collection.set_linestyle('dashed')
                fill_color[0][-1] = 0.6  # make the fill color a bit transparent
                collection.set_facecolor(fill_color)
        if title is None:
            ax.set_title(subtype, fontweight='bold')
        else:
            ax.set_title(title, fontweight='bold')
        ax.set_xlabel('')  # remove the x-axis label
        xlabels = ax.get_xticklabels()  # move x-axis labels to the left
        ax.set_xticklabels(xlabels, rotation=45, ha='right')
        medians = median_counts_subtype[subtype]  # plot median dots
        pos = range(len(medians))
        for tick, median in zip(pos, medians):
            ax.scatter(tick, median, marker='D', color='black', s=10, zorder=3)
        if n > 0:
            ax.set_ylabel('')
    plt.tight_layout()

plot_key = 'total_counts'

# Prepare for plotting
medians, axs = prepare_custom_violin(adata_per_subtype, plot_key=plot_key)

# Plot the Violin
plot_custom_violin(adata_per_subtype, plot_key=plot_key, median_counts_subtype=medians, axs=axs, title='Counts per spot')


# Prepare for adding a legend (cancer vs adj_benign)
labels = list(sample_lut['type'].unique())
colors = ['black'] * len(labels)
handles = []

for n, color in enumerate(colors):
    handles.append(mpatches.Patch(color=color))

# Make function to add the legend with the correct formatting
def add_custom_legend():
    legend = plt.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left', handlelength=1.2, ncol=1)
    # Make the legend boxes look the same as the barplot
    for n, legpatch in enumerate(legend.get_patches()):
            fill_color = legpatch.get_facecolor()
            fill_color = list(fill_color)
            fill_color[-1] = 0.6 # make the fill color a bit transparent
            legpatch.set_color(tuple(fill_color))
            legpatch.set_edgecolor('black')
            if n > 0:
                legpatch.set_linestyle('dashed')
                fill_color[-1] = 0.4 # make the fill color a bit transparent
                legpatch.set_color(tuple(fill_color))
                legpatch.set_edgecolor((0.2,0.2,0.2))

# Add the legend
add_custom_legend()

plt.suptitle(f'Number of counts per spot', fontweight='bold', y=1.02, x=0.5, ha='center')


# --- Plot 2: violin plot log(total counts), split by subtype

plot_key = 'log10_counts'

# Prepare for plotting
medians, axs = prepare_custom_violin(adata_per_subtype, plot_key=plot_key)

# Plot the Violin
plot_custom_violin(adata_per_subtype, plot_key=plot_key, median_counts_subtype=medians, axs=axs, title='Counts per spot')

# Add custom legend
add_custom_legend()

plt.suptitle(f'Log(number of counts) per spot', fontweight='bold', y=1.02, x=0.5, ha='center')


# --- Plot 3: violin plot with genes per spot, split by subtype

plot_key = 'n_genes_by_counts'

# Prepare for plotting
medians, axs = prepare_custom_violin(adata_per_subtype, plot_key=plot_key)

# Plot the Violin
plot_custom_violin(adata_per_subtype, plot_key=plot_key, median_counts_subtype=medians, axs=axs, title='Genes per spot')

# Add custom legend
add_custom_legend()

plt.suptitle(f'Number of genes per spot', fontweight='bold', y=1.02, x=0.5, ha='center')


# --- Plot 4: violin plot with pct_mt per spot, split by subtype

plot_key = 'pct_counts_mt'

# Prepare for plotting
medians, axs = prepare_custom_violin(adata_per_subtype, plot_key=plot_key)

# Plot the Violin
plot_custom_violin(adata_per_subtype, plot_key=plot_key, median_counts_subtype=medians, axs=axs, title='% MT counts')

# Add custom legend
add_custom_legend()

plt.suptitle(f'Percentage of mitochondrial counts per spot', fontweight='bold', y=1.02, x=0.5, ha='center')


# --- Plot 5: Plot the number of spots per gene, split by subtype

plot_key = 'n_cells_by_counts'

# Prepare for plotting
# Set obs to false because this info is stored in .var
medians, axs = prepare_custom_violin(adata_per_subtype, plot_key=plot_key, obs=False)

# Define a function for extracting a df from something in var 
def extract_column_var(plot_key, subtype_in_uns=True, sample_id_key='sample_id'):
    # Extract a dataframe with n_cells_by_counts for every sample, concat per subtype
    n_cells_by_counts_dict = {}
    for sample in adata_dict:
        if subtype_in_uns:
            subtype = adata_dict[sample].uns['metadata']['subtype']
        else:
            subtype = 'PCa'
        tmp = adata_dict[sample].var[plot_key]
        add_sample_col = pd.Series(adata_dict[sample].uns['metadata'][sample_id_key], name = 'sample').repeat(len(tmp))
        tmp_df = pd.concat([tmp.reset_index(drop=True), add_sample_col.reset_index(drop=True)], axis=1)
        tmp_df['sample'] = tmp_df['sample'].astype('category')
        if not subtype in n_cells_by_counts_dict.keys():
            n_cells_by_counts_dict[subtype] = tmp_df
        else:
            n_cells_by_counts_dict[subtype] = pd.concat([n_cells_by_counts_dict[subtype].reset_index(drop=True), tmp_df.reset_index(drop=True)], axis=0)
    return n_cells_by_counts_dict

extracted_df_dict = extract_column_var(plot_key, subtype_in_uns=False)

# Define a function for plotting a violin plot with seaborn (i.e. starting from a df rather than anndata object)

def plot_violin_seaborn(extracted_df_dict, plot_key, medians_per_subtype, axs, title=None, all_subtypes=all_subtypes):
    # Make the actual plot
    for n, (subtype, ax) in enumerate(zip(all_subtypes, axs)):
        sns.violinplot(data = extracted_df_dict[subtype], 
                    x = 'sample', y = plot_key, stripplot=False, ax = ax, title = subtype, linewidth=1, palette=plotting_colors[subtype], 
                    inner = None, order=sorted(extracted_df_dict[subtype]['sample'].unique()), scale='width')
        for i, collection in enumerate(ax.collections):  # Adjust the colors of the violin
            fill_color = collection.get_facecolor()
            darkened_color = fill_color.copy()
            darkened_color[:, :-1] *= 0.8  # multiply RGB values by 0.8 to darken
            collection.set_edgecolor(darkened_color)
            fill_color[0][-1] = 0.8  # make the fill color a bit transparent
            collection.set_facecolor(fill_color)
            if i in highlight_samples_indices:  # If it's an adj_benign sample: mark this sample with slightly lighter color & dashed outline
                collection.set_linestyle('dashed')
                fill_color[0][-1] = 0.6  # make the fill color a bit transparent
                collection.set_facecolor(fill_color)
        if title is None:
            ax.set_title(subtype, fontweight='bold')
        else:
            ax.set_title(title, fontweight='bold')
        # move x-axis tick labels to the left
        xlabels = ax.get_xticklabels()
        ax.set_xticklabels(xlabels, rotation=45, ha='right')
        # remove the x-axis label
        ax.set_xlabel('')
        # plot median dots
        medians = medians_per_subtype[subtype]
        pos = range(len(medians))
        for tick, median in zip(pos, medians):
            ax.scatter(tick, median, marker='D', color='black', s=10, zorder=3)
        if n > 0:
            axs[n].set_ylabel('')

# Plot the violin
plot_violin_seaborn(extracted_df_dict, plot_key=plot_key, medians_per_subtype=medians, axs=axs, title='Spots per gene')

# Add custom legend
add_custom_legend()

plt.suptitle(f'Number of spots per gene', fontweight='bold', y=1.02, x=0.5, ha='center')


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

medians, axs = prepare_custom_violin(adata_per_subtype, plot_key='n_cells_by_counts', obs=False)

# Define a function for creating barplots
def plot_seaborn_barplot(barplot_dict, plot_key, axs, title=None, all_subtypes=all_subtypes):
    for n, (subtype, ax) in enumerate(zip(all_subtypes, axs)):
        sns.barplot(data = barplot_dict[subtype], 
                x='sample', y=plot_key, ax = ax, linewidth=1, palette=plotting_colors[subtype], 
                order=sorted(barplot_dict[subtype]['sample'].unique()))
        for i, bar in enumerate(ax.patches):  # Adjust the colors of the bar
            # Get the RGBA values of the fill color
            fill_color = bar.get_facecolor()
            fill_color = list(fill_color)
            darkened_color = fill_color.copy()
            for j in range(len(darkened_color)-1):
                darkened_color[j] *= 0.8
            fill_color[-1] = 0.8 # make the fill color a bit transparent
            bar.set_color(tuple(fill_color))
            bar.set_edgecolor(tuple(darkened_color))
            sample_index = sorted(barplot_dict[subtype]['sample'].unique())[i]
            if sample_index in highlight_samples:  # If it's an adj_benign sample: mark this sample with slightly lighter color & dashed outline
                bar.set_linestyle('dashed')
                fill_color[-1] = 0.6  # make the fill color a bit transparent
                bar.set_color(fill_color)
        if title is None:
            ax.set_title(subtype, fontweight='bold')
        else:
            ax.set_title(title, fontweight='bold')
        # move x-axis tick labels to the left
        xlabels = ax.get_xticklabels()
        ax.set_xticklabels(xlabels, rotation=45, ha='right')
        # remove the x-axis label
        ax.set_xlabel('')
        if n > 0:
            ax.set_ylabel('')


# Plot the barplot
plot_key = 'n_counts'
plot_seaborn_barplot(barplot_dict, plot_key, axs, title='Total counts')

# Add custom legend
add_custom_legend()

plt.suptitle(f'Total number of counts', fontweight='bold', y=1.02, x=0.5, ha='center')


# --- Plot 6: Plot the number of spots 'in tissue' (barplot)

# Prepare  for plotting
medians, axs = prepare_custom_violin(adata_per_subtype, plot_key='n_cells_by_counts', obs=False)

# Plot the barplot
plot_key = 'in_tissue'
plot_seaborn_barplot(barplot_dict, plot_key, axs, title='Spots in tissue')

# Add custom legend
add_custom_legend()

plt.suptitle(f"Number of spots 'in tissue'", fontweight='bold', y=1.02, x=0.5, ha='center')


# --- Plot 7: Plot the percentage of low quality spots (barplot)

# Prepare  for plotting
medians, axs = prepare_custom_violin(adata_per_subtype, plot_key='n_cells_by_counts', obs=False)

# Plot the barplot
plot_key = 'pct_low_q'
plot_seaborn_barplot(barplot_dict, plot_key, axs, title='Spots in tissue')

# Add custom legend
add_custom_legend()

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
