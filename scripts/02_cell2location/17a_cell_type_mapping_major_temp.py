#!/usr/bin/env python3

# --- Aim
# Use cell2location for deconvolution
# Recycle a pre-built model and do the cell type mapping only

# On a single sample

# this script is exactly the same as scripts 03a and 06a
# copying it to maintain script-results directory structure


# --- Define environment
# env: cell2location

# --- Load packages
import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

from datetime import date
from datetime import datetime

import cell2location

# Set seed for reproducibility
import scvi
scvi.settings.seed = 2023

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

import pandas as pd
import os

import anndata as ad

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

# --- Extract some info from the sample sheet

# Extract the sample_id
sample_dir = os.path.join(repoDir, 'config')
sample_sheet_file = 'sample_sheet.csv'

sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))

# Filter the DataFrame based on the condition sample_id == sample_id
filtered_df = sample_sheet[sample_sheet['section_name'] == section_name]

# Extract the sample_name from the first row (assuming there's only one match)
file_id = filtered_df['sample_name'].iloc[0]

# --- Load the visium data

input_adata_dir = os.path.join(h5ad_dir, 'raw')
adata_vis = ad.read(os.path.join(input_adata_dir, f'{file_id}_raw.h5ad'))

# --- Visium preprocessing

# # Switch out gene_ids for ensembl ids
adata_vis.var['SYMBOL'] = adata_vis.var_names
# adata_vis.var.set_index('gene_ids', drop=True, inplace=True) # comment out this line, Sophie initially used ensembl ids (from BCa mini atlas cellxgene obj)

# find mitochondria-encoded (MT) genes
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]


# --- Load the model

adata_file = "sc.h5ad"
adata_ref = sc.read_h5ad(os.path.join(model_dir, adata_file))
mod = cell2location.models.RegressionModel.load(model_dir, adata_ref)

# --- Record all parameters

# Record parameters
N_cells_per_location=10
detection_alpha=20
batch_size=None
max_epochs=30000
train_size=1
use_gpu=True

parameters_log = {'N_cells_per_location': N_cells_per_location, 'detection_alpha': detection_alpha,'path_to_model': model_dir, 'max_epochs': max_epochs, 
                  'batch_size': str(batch_size), 'train_size':train_size, 'use_gpu': use_gpu}

# Convert dictionary to DataFrame

params = pd.DataFrame(list(parameters_log.items()), columns=['parameter', 'value'])

params.to_csv(os.path.join(tabDir, 'parameters.csv'), index=False)

# --- Export cell type signatures

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()

inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]


# --- Spatial mapping

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key=None)

# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=N_cells_per_location,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=detection_alpha
)

mod.view_anndata_setup()


# --- Train cell2location

mod.train(max_epochs=max_epochs,
          # train using full data (batch_size=None)
          batch_size=batch_size,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=train_size,
          use_gpu=use_gpu,
         )

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training']);

plt.savefig(os.path.join(figDir, 'training_elbo_loss.pdf'))
plt.close()

# --- Export estimated cell abundance

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Save model
mod.save(os.path.join(objectDir, 'model_output'), overwrite=True)

# Save anndata object with results
adata_file = f"{section_name}_cell2location.h5ad"
adata_vis.write(os.path.join(objectDir, adata_file))

# Can load again like this
# adata_vis = sc.read_h5ad(os.path.join(out_dir, 'model_output', adata_file))
# mod = cell2location.models.Cell2location.load(os.path.join(out_dir, 'model_output'), adata_vis)

# Also save the cell type mapping as a csv file for ease of access
# 5% quantile of the posterior distribution, representing the value of cell abundance that the model has high confidence in (aka ‘at least this amount is present’)
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
abundance_df = adata_vis.obs[adata_vis.uns['mod']['factor_names']]
abundance_df.to_csv(os.path.join(objectDir, 'q05_cell_abundance_w_sf.csv'))

# --- Assess mapping quality

mod.plot_QC()
plt.savefig(os.path.join(figDir, 'mapping_quality.pdf'))
plt.close()


# --- Finish

print("Script succesfully completed")