#!/usr/bin/env python3

# --- Aim
# Import visium data into anndata object & and do some basic qc
# Save results as .h5ad file


# --- Define environment
# env: preprocessing

# --- Load packages
import argparse
import scanpy as sc
import pandas as pd
import anndata as ad

import numpy as np
import os

import matplotlib.pyplot as plt


sc.settings.verbosity = 3

# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & save as .h5ad file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
# parser.add_argument("--raw_data_dir", help="Path to raw data")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad input/output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--min_count", help="Minimal number of counts per spot", type=int)
parser.add_argument("--min_gene", help="Minimal number of genes per spot", type=int)
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

# Print the arguments
print(f'.h5ad dir: {h5ad_dir}')
print(f'Sample id: {sample_id}')

print(f'Min count: {min_count}')
print(f'Min gene: {min_gene}')

# Make directories if they don't exist yet
os.makedirs(h5ad_dir, exist_ok=True)


# --- Start
# read in the existing raw.h5ad object (needs to be renamed - this is the filtered unprocessed counts matrix)

in_file = os.path.join(h5ad_dir, 'raw', sample_id + "_raw.h5ad")
adata = ad.read(in_file)

# --- Extract some info from the sample sheet

# Extract the section_id
sample_dir = os.path.join(project_dir, 'config')
sample_sheet_file = 'sample_sheet_mengxiao.csv'

sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))
section_id = sample_sheet[sample_sheet['sample_name'] == sample_id]['sample_id'].values[0]

# Create dictionary to store metadata
metadata_dict = {}
metadata_dict['sample_id'] = section_id
metadata_dict['sample_name'] = sample_id

# Patient id
metadata_dict['patient_id'] = sample_sheet[sample_sheet['sample_name'] == sample_id]['patient_id'].values[0]  #.astype(str)

# Sample type
metadata_dict['sample_type'] = sample_sheet[sample_sheet['sample_name'] == sample_id]['type'].values[0]

# Add the metadata to the anndata object for convenience
adata.uns['metadata'] = metadata_dict

# --- Initial QC

# Annotate the group of mitochondrial genes as 'mt'
adata.var['mt'] = adata.var_names.str.startswith(('MT-'))
print(f"Sample {sample_id}: Mitochondrial genes:", adata.var['mt'].sum())

# Annotate hemoglobin genes as hb
adata.var['hb'] = adata.var_names.str.contains(("^HB[AB]"))
print(f"Sample {sample_id}: Hemoglobin genes:", adata.var['hb'].sum())

# There are no ribosomal genes in the probe set

sc.pp.calculate_qc_metrics(
    adata, qc_vars=['mt','hb'], 
    percent_top=None, log1p=False, inplace=True)


# --- Flag spots with low QC values (but don't filter)

adata.obs['low_qc'] = np.where((adata.obs['total_counts'] <= min_count) | (adata.obs['n_genes_by_counts'] <= min_gene), 'low counts/genes', 'pass')

# Need to do this to make sure the pass isn't colored red if there is none below pass
adata.obs['low_qc'] = pd.Categorical(adata.obs['low_qc'], categories=['pass', 'low counts/genes'], ordered=True)

# Flag spots with low counts
adata.obs['low_counts'] = np.where(adata.obs['total_counts'] <= min_count, 'low_counts', 'pass')
adata.obs['low_counts'] = pd.Categorical(adata.obs['low_counts'], categories=['pass', 'low_counts'], ordered=True)

# Flag spots with low genes
adata.obs['low_gene'] = np.where(adata.obs['n_genes_by_counts'] <= min_gene, 'low_gene', 'pass')
adata.obs['low_gene'] = pd.Categorical(adata.obs['low_gene'], categories=['pass', 'low_gene'], ordered=True)

# Flag spots with 0 counts
adata.obs['no_counts'] = np.where(adata.obs.total_counts == 0.0, 'no counts', 'pass')
adata.obs['no_counts'] = pd.Categorical(adata.obs['total_counts'], categories=['pass', 'no counts'], ordered=True)


# --- Save output

outs_dir = os.path.join(h5ad_dir, 'filtered')

os.makedirs(outs_dir, exist_ok=True)

print(f"Saving .h5ad files to {outs_dir}")

ad.AnnData.write(adata, filename=os.path.join(outs_dir, f"{sample_id}_filtered.h5ad"))

print("Script succesfully completed")