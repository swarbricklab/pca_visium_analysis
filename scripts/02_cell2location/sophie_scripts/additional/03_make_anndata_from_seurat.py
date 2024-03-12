#!/usr/bin/env python3

# --- Aim
# Make anndata objects from Marcos' Seurat objects
# So we're working on the same data
# Based on the extracted .csv files

# Marcos workflow:
#   NormalizeData() |>
#   FindVariableFeatures() |> 
#   ScaleData() |> 
#   RunPCA()


# --- Define environment
# env: preprocessing

# --- Load packages
import argparse
import scanpy as sc
import pandas as pd
import anndata as ad

import numpy as np
import os

sc.settings.verbosity = 3

# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & save as .h5ad file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--raw_data_dir", help="Path to raw data")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--min_count", help="Minimal number of counts per spot", type=int)
parser.add_argument("--min_gene", help="Minimal number of genes per spot", type=int)

args = vars(parser.parse_args())


# --- Arguments
raw_data_dir = args["raw_data_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]
project_dir=args["project_dir"]

min_count=args["min_count"]
min_gene=args["min_gene"]

# --- Extract some info from the sample sheet

# Extract the section_id
sample_dir = os.path.join(project_dir, 'config')
sample_sheet_file = 'sample_sheet.csv'

sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))
external_id = sample_sheet[sample_sheet['section_id'] == sample_id]['external_id'].values[0]

# Also extract the alternative external id (external id that isn't the p number)
external_id_alt = sample_sheet[sample_sheet['section_id'] == sample_id]['external_id_alt'].values[0]

joined_id = f'{external_id}-{sample_id}'


# --- Make an anndata object from the raw data

# Create dictionary to store metadata
metadata_dict = {}

metadata_dict['section_id'] = sample_id
metadata_dict['external_id'] = external_id

# Patient id
metadata_dict['donor_id'] = sample_sheet[sample_sheet['section_id'] == sample_id]['donor_id'].values[0].astype(str)

# Visium experiment id
metadata_dict['vis_exp_id'] = sample_sheet[sample_sheet['section_id'] == sample_id]['vis_exp_id'].values[0]


# --- Specify the path to spaceranger outs

joined_id = f'{external_id}-{sample_id}'
outs_path = os.path.join(raw_data_dir, joined_id, 'outs')

# -- Load data

v2_spatial_file = outs_path + '/spatial/tissue_positions.csv'
v1_spatial_file = outs_path + '/spatial/tissue_positions_list.csv'

# if V1 spatial position file does not exist, create a new one
if os.path.exists(v1_spatial_file) == False:
    v2_spatial_df = pd.read_csv(v2_spatial_file)
    v2_spatial_df.to_csv(path_or_buf = v1_spatial_file, header=None,index=False)

adata = sc.read_visium(
    path=outs_path,
    source_image_path=os.path.join(outs_path, 'spatial'),
    library_id=sample_id
)

adata.var_names_make_unique()

# Add the metadata to the anndata object for convenience
adata.uns['metadata'] = metadata_dict


# --- Load Marcos data

path_to_objects = os.path.join(project_dir, 'data', 'processed', 'seurat_objects', 'post_qc', external_id_alt)

counts_file = f'{external_id_alt}_matrix_counts.csv'
norm_file = f'{external_id_alt}_matrix_data.csv'
scale_file = f'{external_id_alt}_matrix_scale_data.csv'


# --- Update adata.X

counts = pd.read_csv(os.path.join(path_to_objects, counts_file), index_col=0)

# Define a function for prepping extracted matrix for anndata

def prep_matrix_for_adata(df, adata):
    # Transpose
    df = df.T
    # Adjust barcode notation (convert . to -)
    df = df.rename(lambda x: x.replace('.', '-'))
    # Extract lists of genes
    keep_genes = list(df.columns)
    # Deal with genes where .1 needs to become -1 for them to match the anndata object
    genes_not_in_adata = [item for item in keep_genes if item not in list(adata.var_names)]
    # Replace .1 with -1 for the genes that don't match
    new_columns = [col.replace('.1', '-1') if col in genes_not_in_adata else col for col in df.columns]
    df.columns = new_columns
    # Return the adjusted matrix
    return df

counts = prep_matrix_for_adata(counts, adata)

# Extract lists of genes & spots to keep
keep_barcodes = list(counts.index)
keep_genes = list(counts.columns)

# Subset anndata
adata = adata[keep_barcodes, keep_genes]

# Add counts layer to adata.X
adata.X = counts


# --- Add normalized values

norm = pd.read_csv(os.path.join(path_to_objects, norm_file), index_col=0)

# Prepare norm for adding to adata
norm = prep_matrix_for_adata(norm, adata)

adata.layers['log_norm'] = norm

# --- Add scaled values

scaled = pd.read_csv(os.path.join(path_to_objects, scale_file), index_col=0)

# Prepare scaled for adding to adata
scaled = prep_matrix_for_adata(scaled, adata)

# Put in obsm instead of in layers because this is the only highly variable genes (2000 genes)
adata.obsm['highly_var_scaled'] = scaled


# --- While we're at it, also record which genes are highly variable
# Obtained using FindVariableFeatures()

variable_features = list(scaled.columns)

adata.var['highly_variable'] = np.where(adata.var.reset_index()['index'].isin(variable_features), True, False)


# --- Save output

outs_dir = os.path.join(h5ad_dir, 'post_qc_marcos')

os.makedirs(outs_dir, exist_ok=True)

print(f"Saving .h5ad files to {outs_dir}")

ad.AnnData.write(adata, filename=os.path.join(outs_dir, f"{joined_id}_post_qc.h5ad"))

print("Script succesfully completed")
