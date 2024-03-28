# Parameters for submission scripts (submit_universal.sh & submit_dict.sh)

# --- Directories

# 01 PROJECT
projectName="PCa_Visium"
# 02 REPO NAME
repo="pca_visium_analysis"
# 03 EXP CODE
exp="02_cell2location"
# 04 ANALYSIS
analysis="12_train_model_minor_mal"
# 05 SAMPLE ID
sample_id="PCa_atlas"
#  PATH TO SINGLE CELL REFERENCE

## --- Other Args (for QC) 
## 06 MIN COUNT
#min_count=1000
## 07 MIN GENE
#min_gene=200
## min_spots=3

# 08 ENVIRONMENT
environment="cell2location"

# Paths to downstream scripts
h5ad_dir="/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/data/anndata_objects/"
model_dir="/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/02_cell2location/12_train_model_minor_mal/figures/reference_signatures/"

# 09 DATA TO USE (Which data to use - different normalization methods: SCTransform or scanpy's log norm)
# use_data="log_norm"

# make a log directory
logDir="/share/ScratchGeneral/evaapo/projects/${projectName}/${repo}/results/${exp}/${analysis}/logs/"
mkdir -p $logDir

# seed
seed=42

### Visium sample sheet: contains the names of the visium samples
##config_dir="$project_dir/config"
##visium_sample_sheet="$config_dir/sample_sheet.csv"