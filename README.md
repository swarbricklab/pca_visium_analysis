TO DO:  

- [ ] Further deidentify Cansto IDs - need to be removed from scripts and .dvc files!
- [ ] Create a new repo after updates have been made & remove scripts/data not relevant to paper (i.e. analysis of data from PMID 35948798)  
- [ ] Add CELLxGENE links (when available)  
- [ ] Add EGA accession number (when available)
- [ ] Reference scRNA-seq GitHub repo  


# pca_visium_analysis  
This repository contains code for Visium data processing, downstream analysis, and figure generation for the study titled "Profiling of epithelial functional states and fibroblast phenotypes in hormone therapy-naive localized prostate cancer." 

The Visium data was processed using the scanpy package (v1.9.1) and spatial deconvolution was performed using Cell2Location (v0.1.3).  

## Data analysed:  
- In-house Visium samples  
  - Batch 1: 2 samples - Visium FFPE V1 (2 x 6.5x6.5mm reactions)
  - Batch 2: 8 samples - Visium FFPE V1 (8 x 6.5x6.5mm reactions)
  - Batch 3: 4 samples - Visium FFPE V2 (2 x 11x11mm reactions)
- Published data from Joakim's lab (PMID 35948798)

## Brief info on directories below:  

**config/**  
- Contains parameter files for each script, sample sheet & other ad hoc files

**data/**  
- Contains images, anndata objects (filtered, log_norm, raw, raw_beature_bc_matrix), per-spot histopathology annotation (only for samples in the paper)

**environment/**  
- Contains EDF and ESF files for pca_visium  

**resources/**  
- Contains gene lists (e.g., Cepo derived signatures for pnCAFs and Glial cells), and published data (i.e. PMID 35948798)  

**results/**  
- Script outputs are saved to corresponding results dir
  - Includes results related to:
    - Basic processing and QC (found in basic_qc/, clustering/, and marker_genes/)   
    - Cell2Location analysis (found in 02_cell2location/) - the contents of the subdirectories should be self-explanatory.
      - **NOTA BENE** - most analyses related to the PCa manuscript can be found here!
    - Cell2Location colocalisation analyses (inhouse and published PMID 35948798 - not relevant to PCa manuscript), H&E plots, and pnCAF and Glial scanpy scores (found in 04_colocalisation_inhouse_and_published)  

**scripts/**  

01_basic_qc/  
- Initial pre-processing (without filtering spots), QC & visualisation
  
02_cell2location/
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
- 12_train_model_minor_mal.py: Train cell2location model on PCa scRNA-seq reference (cell type minor mal - i.e. a combination of cell type minor annotation and malignant spectrum annotation for Epithelial cells)  
- 13a_cell_type_mapping_minor_mal.py: Map cell types (cell type minor mal)  
- 14_visualize_results_minor_mal.py: Visualize results (cell type minor mal)   
- 15_colocalisation_minor_mal.R: summarise cell2location results, & calculate spot-level cell-cell correlations (grouped by sample id & histological feature; cell type minor mal)   
- 16_train_model_major_temp.py: Train cell2location model on PCa scRNA-seq reference (cell type major temp - i.e. use minor annotation for CAFs and major annotation for all other cell types)  
- 17a_cell_type_mapping_major_temp.py:  Map cell types (cell type major temp)   
- 18_visualize_results_major_temp.py: Visualize results (cell type major temp)   
- 19_colocalisation_major_temp.R: summarise cell2location results, & calculate spot-level cell-cell correlations (grouped by sample id & histological feature; cell type major temp)  

