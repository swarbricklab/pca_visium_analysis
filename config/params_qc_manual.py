# 02_qc_plots_all_samples.py run manually for Batch 3 PCa Visium due to issues with calling scanpy in conda env
# 03_plot_all_unfiltered.py too

project_dir = "/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/"
h5ad_dir = "/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/data/anndata_objects/"
res_dir = "/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/"
sample_list = ['PCa20130_C1_20272_C2', 'PCa20153_C1_20128_C1']
# sample_list = ['20130-1', '20272-2', '20153-1', '20128-1']
use_data = "log_norm"