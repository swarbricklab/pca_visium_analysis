# pca_visium_analysis  
Repo for analysis of the PCa Visium dataset  

Data analysed:  
- In-house Visium samples  
  - Batch 1: 2 samples - Visium FFPE V1  
  - Batch 2: 8 samples - Visium FFPE V2  
- Published data from Joakim's lab (PMID 35948798)    


Brief info on directories below:  

**config/**  
- Contains parameter files for each script & other ad hoc files  

**resources/**  
- Contains published PCa spatial data from Joakim's Nature paper (PMID 35948798)  

**results/**  
- Script outputs are saved to corresponding results dir  


**scripts/**  

01_preprocessing  
- Initial pre-processing (without filtering spots), QC & visualisation
  
02_cell2location
- 01_h5ad_conversion.R: Convert Seurat object to .h5ad  
- 02_train_model.py: Train cell2location model on PCa scRNA-seq reference (cell type minor)  
- 03a_cell_type_mapping.py: Map cell types (cell type minor)  
- 04_visualize_results.py: Visualize results (cell type minor)  
- 05_train_model_major.py: Train cell2location model on PCa scRNA-seq reference (cell type major)  
- 06a_cell_type_mapping_major.py: Map cell types (cell type major)  
- 07_visualize_results_major.py: Visualize results (cell type major)
- 08_visualise_proportions_minor.R
- 09_visualise_proportions_major.R
- 10_colocalisation_major.R: summarise cell2location results, & calculate spot-level cell-cell correlations (grouped by sample id & histological feature; cell type major)  
- 11_colocalisation_minor.R: summarise cell2location results, & calculate spot-level cell-cell correlations (grouped by sample id & histological feature; cell type minor)  
