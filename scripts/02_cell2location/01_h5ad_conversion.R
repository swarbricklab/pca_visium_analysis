# SCRIPT INFORMATION -------------------------------------------------------

# author: Eva Apostolov
# date: 2024-03-12
# last edit: 

# script name: 01_h5ad_conversion.R
# conda env: convert_data

# NOTES ---------------------------------------------------------------------

# Guidelines on data formatting (for stereoscope)

## Single Cell Count Data File --------------------------
## A .h5ad file with cells as rows and genes as columns. Cell type annotations can be read from this file as well, and should then be put in the .obs slot.
## Single Cell Annotation Data --------------------------
## tsv - Columns should be named 'bio_celltype'
## Use labels from .h5ad file. Make sure your labels are located in the .obs slot and then use the flat --label_colname KEY to indicate what the key to access these labels is (i.e., KEY).
## Gene list of top 5000 highly variable genes ----------

# TO DO ----------------------------------------------------------------------

# SETUP ----------------------------------------------------------------------

library(Seurat)
#Sys.setenv(RETICULATE_PYTHON="/home/evaapo/miniconda3/envs/convert_data/bin/python3.10")
library(SeuratDisk)

# set up directory structure
projectName="PCa_Visium"
projectDir=paste0("/share/ScratchGeneral/evaapo/projects/",projectName,"/")
repo="pca_visium_analysis"

# rerun directory - trying to organise PCa analysis
exp = "02_cell2location"
analysis = "01_h5ad_conversion"

# directory structure
resultDir=paste0(projectDir, repo, "/results/", exp, "/", analysis, "/")
rObjectDir=paste0(resultDir, "rObjects/")
figDir=paste0(resultDir, "figures/")
tabDir=paste0(resultDir, "tables/")

# ensure output directories exist
system(paste0("mkdir -p ",resultDir))
system(paste0("mkdir -p ",rObjectDir))
system(paste0("mkdir -p ",figDir))
system(paste0("mkdir -p ",tabDir))

# paths to objects
inFile = paste0("/share/ScratchGeneral/evaapo/projects/PCa/results/20230714_collate_annotation/03_combine_minor_and_malignant/rObjects/PCa_atlas_annotated_Mar24_update.rds")

# START ----------------------------------------------------------------------
# read in merged object 
cells <- readRDS(inFile)

# keep only relevant meta data
to_keep <- c("orig.ident", "barcode", "nCount_RNA", "nFeature_RNA", "percent.mito", "sample_id", "sample",
  "ERG_status_merged", "malignant_anno_merged", "celltype_major_v2", "celltype_minor_v2", "celltype_minor_v2", "celltype_subset_v2", "celltype_minor_mal", "type")

to_remove <- cells@meta.data[, !(colnames(cells@meta.data) %in% to_keep)]
meta_to_remove <- names(to_remove)
  
for(i in meta_to_remove) {
  cells[[i]] <- NULL
}

# to keep meta data names in AnnData obj convert to characters
i <- sapply(cells@meta.data, is.factor)
cells@meta.data[i] <- lapply(cells@meta.data[i], as.character)

# Seurat to .h5ad ------------------------------------------------------------

# bear in mind that if the object has normalised and scaled data, raw counts won't be copied over -.- therefore, remove scale.data slot
## counts = TRUE, data = TRUE, scale.data = FALSE,

cells2 <- DietSeurat(cells)

# convert to h5ad object
# intermediary object
SaveH5Seurat(cells2, filename = paste0(rObjectDir, "PCa_merged_filtered.h5Seurat"), overwrite=TRUE)

# Creating h5Seurat file for version 3.1.5.9900
# Adding counts for RNA
# Adding data for RNA
# Adding variable features for RNA
# Adding feature-level metadata for RNA

# convert
Convert(paste0(rObjectDir, "PCa_merged_filtered.h5Seurat"), dest = "h5ad", filename = paste0(rObjectDir, "PCa_merged_filtered.h5ad"), overwrite=TRUE)

# Validating h5Seurat file
# Adding data from RNA as X
# Transfering meta.features to var
# Adding counts from RNA as raw
# Transfering meta.features to raw/var
# Transfering meta.data to obs

# this object can now be used in cell2Location

# save sessionInfo
writeLines(capture.output(sessionInfo()), paste0(resultDir, "sessionInfo.txt"))

# 12/03/2024 - FYI, in python

## raw counts - made a new adata object - this accesses the raw counts
# adata_1 = adata.raw.to_adata()

## just extracting a pandas df - this is X, it is default -access normalised counts
# adata_2 = adata.to_df()

## to access meta data
# adata.obs[['celltype_minor_v2', 'sample_id']]
## of
# adata.obs['celltype_minor_v2']

# NOTES from 2023 --------

# # go to rObject Dir, open python and check object
# import scanpy
# import pandas as pd
# adata = scanpy.read_h5ad("PCa_merged_filtered.h5ad")
# 
# #adata
# #adata.X
# #adata.raw.X 
# #adata.to_df()
# 
# # check your raw matrix and write a csv file
# t=adata.raw.X.toarray()
# pd.DataFrame(data=t, index=adata.obs_names, columns=adata.raw.var_names).to_csv('adata_raw_x.csv')
# 
# # sparse matrices don't have a native representation, so would need to be converted to csr format
# # for a quick look, do print(adata.X) or print(adata.raw.X)
# 
# # back in R
# 
# # a super complicated way to check if the data frame contains integers only
# # library(data.table)
# # test <- fread('adata_raw_x.csv') # fread is much faster than read.csv
# # library(tibble)
# # df <- tibble::column_to_rownames(test, "V1")
# # 
# # # check for integer
# # is_whole <- function(x) {
# #   all.equal(x, as.integer(x)) # is.integer() checks for integer class, not whole number
# # }
# # 
# # ts <- data.frame(lapply(df, is_whole))
# # ts2 <- unlist(ts, use.names=FALSE)
# # unique(ts2) # all TRUE