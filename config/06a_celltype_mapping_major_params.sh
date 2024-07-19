# Parameters for submission scripts (submit_universal.sh & submit_dict.sh)

# --- Directories

# 01 PROJECT
projectName="PCa_Visium"
# 02 REPO NAME
repo="pca_visium_analysis"
# 03 EXP CODE
#exp="02_cell2location"
exp="02_cell2location"
# 04 ANALYSIS
analysis="06a_cell_type_mapping_major"
# 05 SAMPLE ID
# sample_id="PCa_atlas"

# 08 ENVIRONMENT
environment="cell2location"

# Paths to downstream scripts
h5ad_dir="/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/data/anndata_objects/"
model_dir="/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/${exp}/05_train_model_major/figures/reference_signatures/"

# Visium sample sheet: contains the names of the visium samples
# config_dir="/$project_dir/config"
config_dir="/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/config"
visium_sample_sheet="$config_dir/sample_sheet.csv"