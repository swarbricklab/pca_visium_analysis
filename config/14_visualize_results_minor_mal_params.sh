# Parameters for submission scripts (submit_universal.sh & submit_dict.sh)

# --- Directories

# 01 PROJECT
projectName="PCa_Visium"
# 02 REPO NAME
repo="pca_visium_analysis"
# 03 EXP CODE
# exp="02_cell2location"
exp="03_mengxiao_data_PMID_35948708"
# 04 ANALYSIS
# analysis="14_visualize_results_minor_mal" # change name of analysis
analysis="10_visualize_results_minor_mal" # change name of analysis
# 05 SAMPLE ID
# sample_id="PCa_atlas"

# 08 ENVIRONMENT
environment="cell2location"

# Paths to downstream scripts
h5ad_dir="/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/data/anndata_objects/"
# model_dir="/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/${exp}/12_train_model_minor_mal/figures/reference_signatures/"
model_dir="/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/${exp}/08_train_model_minor_mal/figures/reference_signatures/"