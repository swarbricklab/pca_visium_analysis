# Parameters for submission scripts (submit_universal.sh & submit_dict.sh)

# --- Directories

# 01 PROJECT
projectName="PCa_Visium"
# 02 REPO NAME
repo="pca_visium_analysis"
# 03 EXP CODE
exp="02_cell2location"
# 04 ANALYSIS
analysis="04_visualize_results" # change name of analysis
# 05 SAMPLE ID
# sample_id="PCa_atlas"

# 08 ENVIRONMENT
environment="cell2location"

# Paths to downstream scripts
h5ad_dir="/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/data/anndata_objects/"
model_dir="/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/02_cell2location/02_train_model/figures/reference_signatures/"

# Environment to use
environment="cell2location"