# Parameters for submission scripts (submit_universal.sh & submit_dict.sh)

# --- Directories

cwd=$(pwd)

# Get the project directory (root)
# Use dirname to get one directory up (assuming the scripts folder is one down from the root)
project_dir="/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/"
h5ad_dir="$project_dir/data/anndata_objects" # Path to the h5ad files (anndata objects)

## res_dir="$project_dir/results"

log_dir="$project_dir/logs"
mkdir -p $log_dir

# Other arguments (for qc)
min_count=1000
min_gene=200
min_spots=3

# Environment to use
environment="pca-visium"

# Seed
seed=42

# Which data to use - different normalization methods: SCTransform or scanpy's log norm
use_data="log_norm" # which data to use (either 'SCT' or 'log_norm')

# Visium sample sheet: contains the names of the visium samples
config_dir="$project_dir/config"
visium_sample_sheet="$config_dir/sample_sheet_mengxiao.csv"