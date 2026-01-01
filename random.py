# --- Arguments
project_dir = '/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis'
sample_id = 'P25'
use_data = 'log_norm'

# data_dir = args["data_dir"]
res_dir = os.path.join(project_dir, 'results')
h5ad_dir = os.path.join(project_dir, 'data', 'anndata_objects')
log_dir = os.path.join(project_dir, 'logs')
config_dir = os.path.join(project_dir, 'config')
visium_sample_sheet = os.path.join(config_dir, 'sample_sheet.csv')


# Other arguments (for qc)
min_count=1000
min_gene=200
min_spots=3