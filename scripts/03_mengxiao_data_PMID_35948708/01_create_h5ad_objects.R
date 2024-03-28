# SCRIPT INFORMATION -------------------------------------------------------

# author: Eva Apostolov
# date: 2024-03-27
# last edit: 

# script name: 01_create_h5ad_objects.R
# conda env: convert_data

# NOTES ---------------------------------------------------------------------

# I have downloaded Mengxiao's PCa Visium data (PMID 35948708) from Mendeley Data
# I saved this manually to the HPC
# I need to reconstruct objects (easier done in Seurat) and then convert them to .h5ad objects to run cell2location on

# SETUP ----------------------------------------------------------------------

library(Seurat)
library(SeuratDisk)
library(tibble)
library(dplyr)

# ARGS ------
temp_args <- commandArgs(trailingOnly = TRUE)
# 01 PROJECT
projectName <- temp_args[1]
# 02 REPO
repo <- temp_args[2]
# 03 EXP CODE
exp <- temp_args[3]
# 04 ANALYSIS
analysis <- temp_args[4]
# 05 SAMPLE ID
sample_id <- temp_args[5]

# set up directory structure
projectDir=paste0("/share/ScratchGeneral/evaapo/projects/",projectName,"/")
resultDir=paste0(projectDir, "/", repo, "/results/", exp, "/", analysis, "/")
rObjectDir=paste0(resultDir, "rObjects/")
figDir=paste0(resultDir, "figures/")
# manually created
annDataDir=paste0(projectDir, "pca_visium_analysis/data/PMID_35948708/filtered_feature_bc_matrix/")

# ensure output directories exist
system(paste0("mkdir -p ",resultDir))
system(paste0("mkdir -p ",rObjectDir))
system(paste0("mkdir -p ",figDir))
system(paste0("mkdir -p ",annDataDir))

# # one time thing:
# # restructure directories, i.e. copy all images to counts_inDir so they are found in corresponding sample folders
# sample_ids <- c("H1_2", "H1_4", "H1_5", "H2_1", "H2_2", "H2_5", "V1_2")
# 
# for(sample_id in sample_ids){
#   print(sample_id)
# 
#   counts_inDir = paste0("/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/resources/published_data/PMID_35948708/Count_matrices/Patient_1/Visium_with_annotation/", sample_id, "/")
# 
#   # Assuming images_inDir is already defined, if not, define it accordingly.
#   images <- list.files(images_inDir, pattern=".png", full.names=TRUE)
# 
#   image_file <- grep(sample_id, images, value=TRUE)
#   
#   # Copying the image to counts_inDir
#   file.copy(image_file, counts_inDir)
# }


# paths to objects - from patient 1 only for now
counts_inDir = paste0("/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/resources/published_data/PMID_35948708/Count_matrices/Patient_1/Visium_with_annotation/", sample_id, "/")

# sample ids (i.e. section_ids)
# sample_ids <- c("H1_1", "H1_2", "H1_4", "H1_5", "H2_1", "H2_2", "H2_5", "V1_1", "V1_2")
# only these have latest Visium & per-spot annotation
# sample_ids <- c("H1_2", "H1_4", "H1_5", "H2_1", "H2_2", "H2_5", "V1_2")

# START ----------------------------------------------------------------------

image <- Read10X_Image(counts_inDir, image.name = paste0(sample_id, "_tissue_hires_image.png"))
seurat <- Load10X_Spatial(counts_inDir, filename = "filtered_feature_bc_matrix.h5", image = image, slice=sample_id)
seurat@images[[sample_id]]@scale.factors$lowres <- seurat@images[[sample_id]]@scale.factors$hires

# add per spot histological annotation
anno <- read.csv(paste0(counts_inDir, sample_id, "_Final_Consensus_Annotations.csv"))
anno <- anno %>% remove_rownames %>% column_to_rownames(var="Barcode")

seurat <- AddMetaData(seurat, anno)

p <- SpatialDimPlot(seurat, pt.size.factor = 4, group.by = "Final_Annotations")
pdf(paste0(figDir, sample_id,  "_histo_path.pdf"))
print(p)
dev.off()

# now convert to .h5ad objects
annData_file <- paste0(annDataDir, "/", sample_id, "_raw_feature_bc_matrix.h5ad")

if(!file.exists(paste0(annDataDir, "/", sample_id, "_raw_feature_bc_matrix.h5ad"))){

  # convert to h5ad object
  # intermediary object
  SaveH5Seurat(seurat, filename = paste0(annDataDir, sample_id, "_raw_feature_bc_matrix.h5Seurat"), overwrite=TRUE)

  # Creating h5Seurat file for version 3.1.5.9900
  # Adding counts for RNA
  # Adding data for RNA
  # Adding variable features for RNA
  # Adding feature-level metadata for RNA
  
  # convert
  Convert(paste0(annDataDir, sample_id, "_raw_feature_bc_matrix.h5Seurat"), dest = "h5ad", filename = paste0(annDataDir, "/", sample_id, "_raw_feature_bc_matrix.h5ad"), overwrite=TRUE)
  
  # Seurat to .h5ad ------------------------------------------------------------
  
  # bear in mind that if the object has normalised and scaled data, raw counts won't be copied over -.- therefore, remove scale.data slot
  ## counts = TRUE, data = TRUE, scale.data = FALSE,

  # Validating h5Seurat file
  # Adding data from RNA as X
  # Transfering meta.features to var
  # Adding counts from RNA as raw
  # Transfering meta.features to raw/var
  # Transfering meta.data to obs
  
  # this object can now be used in cell2Location

} else {
  message(".h5ad file already exists")
}

# SAVE LOGS -----------------------------------------------------------------

# save sessionInfo
writeLines(capture.output(sessionInfo()), paste0(resultDir, "sessionInfo.txt"))