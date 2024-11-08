# 02_qc_plots_all_samples.py run manually for Batch 3 PCa Visium due to issues with calling scanpy in conda env
# 03_plot_all_unfiltered.py too

project_dir = "/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/"
h5ad_dir = "/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/data/anndata_objects/"
res_dir = "/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/"
sample_list = ['PCaP4_C1_P13_C2', 'PCaP10_C1_P3_C1']
# sample_list = ['P4-1', 'P13-2', 'P10-1', 'P3-1']
use_data = "log_norm"