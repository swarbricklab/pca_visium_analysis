# This scripts contains paths for quick interactive testing of scripts on specific samples

import os

# --- Arguments
project_dir = "/share/ScratchGeneral/sopvan/datasets/pca_visium_eva"

data_dir = os.path.join(project_dir, 'data')
res_dir = os.path.join(project_dir, "results")
h5ad_dir = os.path.join(project_dir, "data/processed/anndata_objects")
sample_id = "H1_5"
raw_data_dir=os.path.join(project_dir, 'data', 'processed', 'spaceranger_outs')



min_count=1000
min_gene=100

min_spots=3

use_data='log_norm'             # which data normalisation to use for clustering (SCT or log_norm)

today = '20240415'


# outs_path should contain tissue_positions_list.csv
# output obj is saved to: outs_dir = os.path.join(h5ad_dir, 'raw')





# Required
parser.add_argument("--raw_data_dir", help="Path to raw data")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--min_count", help="Minimal number of counts per spot", type=int)
parser.add_argument("--min_gene", help="Minimal number of genes per spot", type=int)
parser.add_argument("--outs_path", help="Path to spaceranger 'outs'")
parser.add_argument("--use_data", help="Which data to use (log_norm or SCT)")


args = vars(parser.parse_args())


# ---------------------------------------------------------
seed=42

# --- Arguments
# raw_data_dir = args["raw_data_dir"]
# raw_data_dir = "/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/data/anndata_objects/raw/" # - DO NOT HAVE RAW DATA
# in_dir - read in existing objects
# h5ad_dir=args["h5ad_dir"]
raw_data_dir = h5ad_dir
# res_dir = os.path.join(project_dir, "results")
res_dir = ""
# sample_id=args["sample_id"]
# project_dir=args["project_dir"]
# min_count=args["min_count"]
min_count=1000
# min_gene=args["min_gene"]

# ---
min_spots=3
use_data='log_norm'             # which data normalisation to use for clustering (SCT or log_norm)
today = '20240415'

project_dir = "/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/"
in_dir = "/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/data/anndata_objects/raw"
h5ad_dir = "/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/data/anndata_objects/"

sample_id = "H1_5"
min_count=1000
min_gene=100
