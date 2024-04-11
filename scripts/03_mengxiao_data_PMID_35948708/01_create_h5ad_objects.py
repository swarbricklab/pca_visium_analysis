# date: 2024-04-09
# conda env: cell2location

# read in Mengxiao's data (filtered counts matrices), add images, and save as h5ad objects

import os
import pandas as pd
import json
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

# directory structure
projectName="PCa_Visium"
repo="pca_visium_analysis"
exp="03_mengxiao_data_PMID_35948708"
analysis="01_create_h5ad_objects"

# # --- Argparse arguments                
# parser = argparse.ArgumentParser(description="import variables",
#                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# 
# parser.add_argument("--sample_id", help="Sample ID")
# sample_id = args["sample_id"]

# Define directory paths
repoDir = os.path.join('/share/ScratchGeneral/evaapo/projects/', projectName, repo)
resultDir = os.path.join(repoDir, "results", exp, analysis)
figDir = os.path.join(resultDir, "figures")
h5ad_dir = os.path.join(repoDir, 'data', 'anndata_objects') # save objects here

# Ensure output directories exist
os.makedirs(resultDir, exist_ok=True)
os.makedirs(figDir, exist_ok=True)
os.makedirs(h5ad_dir, exist_ok=True)

# Define the directory path for the current sample
# sample_ids = ["H1_2", "H1_4", "H1_5", "H2_1", "H2_2", "H2_5", "V1_2"]

sample_id = "H1_2"
inDir = os.path.join('/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/resources/published_data/PMID_35948708/Count_matrices/Patient_1/Visium_with_annotation/', sample_id)

# Read the data
adata_mengxiao = sc.read_visium(path=inDir, load_images=False)

adata_mengxiao.var_names_make_unique()

# Read tissue positions
tissue_positions = pd.read_csv(os.path.join(inDir, 'tissue_positions_list.csv'), header=None, index_col=0)
sorted_names = adata_mengxiao.obs_names
tissue_positions = tissue_positions.reindex(sorted_names)

# Create coordinates array
coordinates = tissue_positions[[4, 5]].to_numpy()
adata_mengxiao.obsm['spatial'] = coordinates

# Load the image
img = plt.imread(os.path.join(inDir, sample_id + '_tissue_hires_image.png'))

# Create a figure and plot the image
fig, ax = plt.subplots(nrows=1, ncols=1)
plt.imshow(img)
plt.savefig(os.path.join(figDir, sample_id + '_img.pdf'))

# Add image to the anndata object
library_id = list(adata_mengxiao.uns['spatial'].keys())[0]
image_key = 'hires'
adata_mengxiao.uns['spatial'][library_id]['images'] = {image_key: img}

# Read scale factors
json_file = 'scalefactors_json.json'
with open(os.path.join(inDir, json_file), 'r') as file:
    data = json.load(file)

print(data) 

adata_mengxiao.uns['spatial'][library_id]['scalefactors'] = data

# Plot spatial
sc.pl.spatial(adata_mengxiao, color='PCA3')
plt.savefig(os.path.join(figDir, sample_id + '_scale_factors_test.pdf'), bbox_inches='tight')

# Transform the image
img = np.fliplr(img)
img = np.rot90(img)
adata_mengxiao.uns['spatial'][library_id]['images'][image_key] = img

# Plot spatial with transformed image
sc.pl.spatial(adata_mengxiao, color='PCA3', vmax='p95')  # scale.factor = 1, with original image
plt.savefig(os.path.join(figDir, sample_id + '_final_position.pdf'), bbox_inches='tight')

# add histopath annotations
histo_annotations = pd.read_csv(os.path.join(inDir, sample_id + "_Final_Consensus_Annotations.csv"), header=None, index_col=0)
adata_mengxiao.obs_names
adata_mengxiao.obs['Histology'] = histo_annotations

# Plot with histopath annotations
sc.pl.spatial(adata_mengxiao, color='Histology', vmax='p95')  # scale.factor = 1, with original image
plt.savefig(os.path.join(figDir, sample_id + '_histology.pdf'), bbox_inches='tight')

adata_file = sample_id + "_raw.h5ad"
adata_mengxiao.write(os.path.join(h5ad_dir, 'raw', adata_file))