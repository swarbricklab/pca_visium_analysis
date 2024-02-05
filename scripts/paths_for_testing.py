# This scripts contains paths for quick interactive testing of scripts on specific samples

import os

# --- Arguments
project_dir = "/share/ScratchGeneral/sopvan/projects/pca_visium_analysis"

data_dir = os.path.join(project_dir, 'data')
res_dir = os.path.join(project_dir, "results")
h5ad_dir = os.path.join(project_dir, "data/anndata_objects")
sample_id = "19617-2"

sample_list = ['20111-1', '20384-2',
'20216-1',
'19617-2',
'20153-1',
'20153-2',
'20033',
'20111-2',
'20130-1',
'20130-2']

seed=42

min_count=1000
min_gene=100

min_spots=3

use_data='log_norm'             # which data normalisation to use for clustering (SCT or log_norm)

today = '20230309'