#!/usr/bin/env python3

# --- Aim
# Train cell2location model to be used on many samples

# --- Define environment
# env: cell2location

# mapping can be performed with scripts 03a... or 06a... (they are the same script, so I won't duplicate it again)

# --- Load packages
import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import os

from datetime import date

# this line should go before importing cell2location
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'
import cell2location

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

# Some imports for qc plots
import numpy as np
import matplotlib
from scvi import REGISTRY_KEYS
from scipy.sparse import issparse


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import variables & train Cell2Location model",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--projectName", help="Main project name")
parser.add_argument("--repo", help="Repo name")
parser.add_argument("--exp", help="Exp name")
parser.add_argument("--analysis", help="Analysis step name")
parser.add_argument("--sample_id", help="Sample ID")

args = vars(parser.parse_args())

# Create directory structure using argparse arguments
projectName = args["projectName"]
repo = args["repo"]
exp = args["exp"]
analysis = args["analysis"]
sample_id = args["sample_id"]

# Define directory paths
repoDir = os.path.join('/share/ScratchGeneral/evaapo/projects/', projectName, repo)
resultDir = os.path.join(repoDir, "results", exp, analysis)
rObjectDir = os.path.join(resultDir, "rObjects")
figDir = os.path.join(resultDir, "figures")
tabDir = os.path.join(resultDir, "tables")
qcDir = os.path.join(resultDir, 'model')


# Ensure output directories exist
os.makedirs(repoDir, exist_ok=True)
os.makedirs(resultDir, exist_ok=True)
os.makedirs(rObjectDir, exist_ok=True)
os.makedirs(figDir, exist_ok=True)
os.makedirs(tabDir, exist_ok=True)
os.makedirs(qcDir, exist_ok=True)


# today = date.today()
# today = today.strftime("%Y%m%d")
# out_dir = os.path.join(qc_dir, today)
# os.makedirs(out_dir, exist_ok=True)

# --- Load the single-cell ref

# Use the miniatlas downloaded from cellxgene
pca_atlas_dir = os.path.join(repoDir, 'results', '02_cell2location', '01_h5ad_conversion', 'rObjects')
pca_atlas_file = 'PCa_merged_filtered.h5ad'

adata_ref = sc.read(os.path.join(pca_atlas_dir, pca_atlas_file))

# Extract the raw data (cell2location wants raw data)
adata_ref = adata_ref.raw.to_adata()

# --- Filter genes

from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_ref = adata_ref[:, selected].copy()

plt.savefig(os.path.join(figDir, 'genes_filtered.pdf'))

plt.close()

# --- Prepare for regression model

# Arguments
labels_key = 'celltype_minor_mal'
donor_id = 'sample_id'
batch_key = 'orig.ident' # this is the capture reaction (PCa1:PCa16)

# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key=batch_key,
                        # cell type, covariate used for constructing signatures
                        labels_key=labels_key,
                        # multiplicative technical effects (donor effect)
                        categorical_covariate_keys=[donor_id]
                       )


# create the regression model
# Note: ref needs to be raw (i.e. untransformed, unnormalized) for this to work
# If it's not raw this step will throw an error
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

# view anndata_setup as a sanity check
mod.view_anndata_setup()


# --- Training the model
mod.train(max_epochs=250, use_gpu=True)

# Determine if model needs more training
mod.plot_history(20)

plt.savefig(os.path.join(figDir, 'model_training.pdf'))
plt.close()

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# Save model
mod.save(os.path.join(figDir, f"reference_signatures"), overwrite=True)

# Save anndata object with results
adata_file = f"reference_signatures/sc.h5ad"
adata_ref.write(os.path.join(figDir, adata_file))


# --- Examine QC plots

# Note: this plot doesn't seem to be plotted correctly, but not sure if I can do anything about it
# Should be two separate plots but they're plotted on top of each other

# mod.plot_QC()
# plt.savefig(os.path.join(figDir, 'model_qc.pdf'))
# plt.close()

# Actually, take the source code to stop them being plotted on top of each other

# Some default params
summary_name = 'means'
use_n_obs = 1000
scale_average_detection = True

# Plot 1: Reconstruction accuracy
# Should be roughly diagonal

if getattr(mod, "samples", False) is False:
    raise RuntimeError("self.samples is missing, please run self.export_posterior() first")

if use_n_obs is not None:
    ind_x = np.random.choice(
        mod.adata_manager.adata.n_obs, np.min((use_n_obs, mod.adata.n_obs)), replace=False
    )
else:
    ind_x = None

mod.expected_nb_param = mod.module.model.compute_expected(
    mod.samples[f"post_sample_{summary_name}"], mod.adata_manager, ind_x=ind_x
)
x_data = mod.adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)[ind_x, :]

if issparse(x_data):
    x_data = np.asarray(x_data.toarray())

mod.plot_posterior_mu_vs_data(mod.expected_nb_param["mu"], x_data)

plt.savefig(os.path.join(qcDir, 'model_qc-1.pdf'))

plt.close()


# Plot 2: Estimated reference expression signatures comp. to average expression in cluster
# Should look like a perfect diagonal

inf_aver = mod.samples[f"post_sample_{summary_name}"]["per_cluster_mu_fg"].T
if scale_average_detection and ("detection_y_c" in list(mod.samples[f"post_sample_{summary_name}"].keys())):
    inf_aver = inf_aver * mod.samples[f"post_sample_{summary_name}"]["detection_y_c"].mean()

aver = mod._compute_cluster_averages(key=REGISTRY_KEYS.LABELS_KEY)
aver = aver[mod.factor_names_]

plt.hist2d(
    np.log10(aver.values.flatten() + 1),
    np.log10(inf_aver.flatten() + 1),
    bins=50,
    norm=matplotlib.colors.LogNorm(),
)
plt.xlabel("Mean expression for every gene in every cluster")
plt.ylabel("Estimated expression for every gene in every cluster")
plt.show()

plt.savefig(os.path.join(qcDir, 'model_qc-2.pdf'))

plt.close()

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

# Save this as a table
inf_aver.to_csv(os.path.join(tabDir, 'reference_cell_type_signatures.csv'))

# --- Finish

print("Script succesfully completed")