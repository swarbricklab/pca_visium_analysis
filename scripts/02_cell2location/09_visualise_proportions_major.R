# SCRIPT INFORMATION -------------------------------------------------------

# author: Eva Apostolov
# date: 2024-03-15
# last edit: 

# script name: 08_visualise_proportions.R
# conda env: ggpubr_v2

# NOTES ---------------------------------------------------------------------

# play around with output of cell2location - figure out how to visualise cell2location as proportions

# TO DO ----------------------------------------------------------------------

# SETUP ----------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(RColorBrewer)
library(corrplot)

# set up directory structure
projectName="PCa_Visium"
projectDir=paste0("/share/ScratchGeneral/evaapo/projects/",projectName,"/")
repo="pca_visium_analysis"

# rerun directory - trying to organise PCa analysis
exp = "02_cell2location"
analysis = "09_visualise_proportions_major"

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

# START ----------------------------------------------------------------------

# read in histopath data as a df ----------

# link to histo annotations (copied manually from Dropbox) - only 6 files have it:
histo_path <- paste0(projectDir, "/", repo, "/resources/cell2location/20240318_reviewed_histo_annotation/")

histo_df <- NULL

for(file in unique(list.files(histo_path, full.name=TRUE))){
  print(file)

  temp_df <- read.csv(file)

  # Remove everything before one letter one digit pattern at the start (including the pattern) and everything after "_Histology_Reviewed.csv"
  sample_id <- sub(".*/([A-Za-z][0-9])_", "", file)
  sample_id <- sub("_Histology_Reviewed.csv", "", sample_id)

  temp_df$sample_id <- sample_id

  histo_df <- rbind(histo_df, temp_df)

}

# > table(histo_df$sample_id)
# 
# 19617-2   20033 20111-2 20130-2 20153-2 20216-1
#    1546    2941    2213    1477    1248    1718

# path to config file used in cell2location - for sample ids --------------------
configFile = paste0(projectDir, repo, "/config/sample_sheet.csv")

samples_config <- read.csv(configFile)

# don't include Batch 1 samples, they are low QC
samples_config_sub <- samples_config %>% filter(batch == 2)

df <- NULL

for (sample_id in unique(samples_config_sub$sample_id)){

  print(sample_id)

  # path to cell2location (Cell type major) - per sample
  prop_df <- read.csv(paste0("/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/02_cell2location/06a_cell_type_mapping_major/", sample_id, "/objects/q05_cell_abundance_w_sf.csv"))
  prop_df$sample_id <- sample_id

  # Extracting the type corresponding to the sample_id
  core_type <- samples_config_sub$type[samples_config_sub$sample_id == sample_id]
  prop_df$type <- core_type
  
  # prop_df <- prop_df %>% remove_rownames %>% column_to_rownames(var="X")
  names(prop_df)[names(prop_df) == "X"] <- "Barcode"

  df <- rbind(df, prop_df)

}

# Given vector of column names to exclude
exclude_cols <- c("Barcode", "sample_id", "type")

# Extract column names not found in exclude_cols
celltypes <- names(df)[!names(df) %in% exclude_cols]

df_tidy <- df %>% pivot_longer(celltypes, names_to = "celltype_major_v2", values_to = "c2l_value")

# Count how many unique barcodes per sample_id
barcode_counts <- df_tidy %>%
  group_by(sample_id) %>%
  summarise(unique_barcode_count = n_distinct(Barcode))

# > barcode_counts
# # A tibble: 8 x 2
#   sample_id unique_barcode_count
#   <chr>                    <int>
# 1 19617-2                   1546
# 2 20033                     2941
# 3 20111-2                   2213
# 4 20130-1                   2180
# 5 20130-2                   1477
# 6 20153-1                   1517
# 7 20153-2                   1248
# 8 20216-1                   1718

# combine the two dfs --------------------
# keeping only the 6 samples that have histo anno

# Merge based on both sample_id and Barcode
merged_df <- inner_join(df_tidy, histo_df, by = c("sample_id", "Barcode"))

# > dim(merged_df)
# [1] 579436      6

# plot some simple dfs --------------------

# MINOR ----------------------------------------------------------------------------

# 1 ------------
# start plotting cell2location results 

celltype_prop_df <- merged_df %>%
  dplyr::group_by(sample_id, Barcode, celltype_major_v2, type) %>%
  dplyr::summarise(c2l_cell_number = sum(c2l_value)) %>%
  dplyr::group_by(sample_id, Barcode) %>%
  dplyr::mutate(total_c2l_cell_number = sum(c2l_cell_number)) %>%
  dplyr::mutate(celltype_proportion = c2l_cell_number / total_c2l_cell_number) %>%
  dplyr::select(Barcode, sample_id, celltype_major_v2, celltype_proportion, type) %>%
  dplyr::filter(celltype_proportion >= 0.05) # only include cell types that constitute at least 5% of a spot

# > dim(merged_df)
# [1] 579436      6

# > dim(celltype_prop_df)
# [1] 41917     5

# Generate a set of colors using a ColorBrewer palette
nice_colors_func <- colorRampPalette(brewer.pal(9, "Set1"))
nice_colors <- nice_colors_func(12)

df_plot <- celltype_prop_df

# plot per core type
p <- ggplot(df_plot, aes(x = type, y = celltype_proportion, fill = celltype_major_v2)) +
  geom_bar(position = "fill", stat="identity") +
  labs(title = "Cell type major composition (Cell2Location) per Core Type",
       subtitle = " Cell types constituting <5% of a spatial spot are filtered out",
       x = "Cell type major (c2l decon)",
       y = "Proportions") +
  scale_fill_manual(values = nice_colors) +  # Apply ColorBrewer colors
  theme_minimal()

pdf(paste0(figDir, "cell2location_per_core_type.pdf"))
print(p)
dev.off()  

# plot per sample (based on above calculations)
agg <- df_plot %>%
  filter(celltype_major_v2 == "Epithelial") %>%
  group_by(sample_id) %>%
  summarise(freq = sum(celltype_proportion)) %>%
  arrange(desc(freq)) # arrange by decreasing frequency

## add missing levels and set the order of celltype_subset_v2
df_plot$sample_id <- factor(df_plot$sample_id, levels = agg$sample_id)

df_plot <- df_plot %>% group_by(sample_id) %>% arrange(desc(celltype_proportion))

p <- ggplot(df_plot, aes(x = as.factor(sample_id), y = celltype_proportion, fill = celltype_major_v2)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(title = "Cell type major composition (Cell2Location) per Sample",
       subtitle = " Cell types constituting <5% of a spatial spot are filtered out",
       x = "Cell type major (c2l decon)",
       y = "Proportions") +
  scale_fill_manual(values = nice_colors) +
  facet_wrap(~type, scales="free") +
  theme_minimal()

pdf(paste0(figDir, "cell2location_anno_per_sample.pdf"), width=10)
print(p)
dev.off()  

# plot per histological region, per sample
celltype_prop_df <- merged_df %>%
  dplyr::group_by(sample_id, Histology, Barcode, celltype_major_v2) %>%
  dplyr::summarise(c2l_cell_number = sum(c2l_value)) %>%
  dplyr::group_by(sample_id, Histology, Barcode) %>%
  dplyr::mutate(total_c2l_cell_number = sum(c2l_cell_number)) %>%
  dplyr::mutate(celltype_proportion = c2l_cell_number / total_c2l_cell_number) %>%
  dplyr::select(sample_id, Histology, Barcode, celltype_major_v2, celltype_proportion) %>%
  dplyr::filter(celltype_proportion >= 0.05) %>% 
  dplyr::filter(Histology != "Exclude")

df_plot <- celltype_prop_df

# agg <- df_plot %>%
#  filter(celltype_major_v2 == "Luminal") %>%
#  group_by(sample_id) %>%
#  summarise(freq = sum(celltype_proportion)) %>%
#  arrange(desc(freq)) # arrange by decreasing frequency
# 
# ## add missing levels and set the order of celltype_subset_v2
# df_plot$sample_id <- factor(df_plot$sample_id, levels = agg$sample_id)

p <- ggplot(df_plot, aes(x = Histology, y = celltype_proportion, fill = celltype_major_v2)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(title = "Cell2Location decon (Cell type major) per Sample",
       subtitle = " Cell types constituting <5% of a spatial spot are filtered out",
       x = "Cell type major (c2l decon)",
       y = "Proportions") +
  scale_fill_manual(values = nice_colors) +
  facet_wrap(~as.factor(sample_id), scales="free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0(figDir, "cell2location_anno_per_histo_region_per_sample.pdf"), width=12)
print(p)
dev.off()  


# now plot populations of interest ----

celltype_prop_df <- merged_df %>%
  dplyr::group_by(sample_id, Histology, Barcode, celltype_major_v2) %>%
  dplyr::summarise(c2l_cell_number = sum(c2l_value)) %>%
  dplyr::group_by(sample_id, Histology, Barcode) %>%
  dplyr::mutate(total_c2l_cell_number = sum(c2l_cell_number)) %>%
  dplyr::mutate(celltype_proportion = c2l_cell_number / total_c2l_cell_number) %>%
  dplyr::select(sample_id, Histology, Barcode, celltype_major_v2, celltype_proportion) %>%
  dplyr::filter(celltype_proportion >= 0.05) %>% # only include cell type proportions with a value greater than 10% (0.01)
  dplyr::filter(Histology != "Exclude") %>%
  # dplyr::filter(celltype_major_v2 %in% c("Luminal", "uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "pSMCs", "Pericytes", "vSMCs", "Cycling_Epi_Luminal", "Glial_cells_1", "Glial_cells_2", "Glial_cells_3", "Satellite_cells", "Myelinating_Schwann"))
  dplyr::filter(celltype_major_v2 %in% c("Epithelial", "CAFs", "SMCs", "Cycling", "PNS_glial"))

df_plot <- celltype_prop_df

nice_colors_func <- colorRampPalette(brewer.pal(9, "Set1"))
nice_colors <- nice_colors_func(16)

p <- ggplot(df_plot, aes(x = Histology, y = celltype_proportion, fill = celltype_major_v2)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(title = "Cell2Location decon (Cell type major) per Histological region, per Sample",
       subtitle = " Cell types constituting <5% of a spatial spot are filtered out; showing only cell types of interest",
       x = "Cell type major (c2l decon)",
       y = "Proportions") +
  scale_fill_manual(values = nice_colors) +
  facet_wrap(~as.factor(sample_id), scales="free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0(figDir, "cell2location_anno_per_histo_region_per_sample_of_interest_5_percent.pdf"), width=9)
print(p)
dev.off()  

# calculate spatial correlations ------------------------------
df_plot <- merged_df 
df_plot$Barcode_sample_id <- paste0(df_plot$Barcode, "_", df_plot$sample_id)


## stroma ---
wanted_celltypes <- c("Epithelial", "CAFs", "SMCs", "Cycling", "PNS_glial")
histo_feature <- c("Stroma prostatic", "Stroma extraprostatic")

df_filt <- df_plot %>% dplyr::filter(celltype_major_v2 %in% all_of(wanted_celltypes))
df_sub <- df_filt %>% filter(Histology %in% all_of(histo_feature))
df_stroma <- df_sub %>% dplyr::select(-any_of(c("Barcode", "sample_id", "type", "Histology"))) %>% pivot_wider(names_from = "celltype_major_v2", values_from = "c2l_value") %>% remove_rownames %>% column_to_rownames(var="Barcode_sample_id")


plot_df <- cor(df_stroma)
# invert default colours
pdf(paste0(figDir, "stroma_major.pdf"))
p <- corrplot(plot_df, title="STROMA - all samples combined", mar=c(0,0,1,0), tl.cex = 0.95, tl.col='black', order='AOE', cl.pos='b', addCoef.col = 'black', number.cex=0.8, col=colorRampPalette(c("dodgerblue4","white","red4"))(200), col.lim=c(-1,1))
print(p)
dev.off()


## vessel ---
wanted_celltypes <- c("Epithelial", "CAFs", "SMCs", "Cycling", "PNS_glial")
histo_feature <- c("Vessel")

df_filt <- df_plot %>% dplyr::filter(celltype_major_v2 %in% all_of(wanted_celltypes))
df_sub <- df_filt %>% filter(Histology %in% all_of(histo_feature))
df_vessel <- df_sub %>% dplyr::select(-any_of(c("Barcode", "sample_id", "type", "Histology"))) %>% pivot_wider(names_from = "celltype_major_v2", values_from = "c2l_value") %>% remove_rownames %>% column_to_rownames(var="Barcode_sample_id")


plot_df <- cor(df_vessel)
# invert default colours
pdf(paste0(figDir, "vessel_major.pdf"))
p <- corrplot(plot_df, title="VESSEL - all samples combined", mar=c(0,0,1,0), tl.cex = 0.95, tl.col='black', order='AOE', cl.pos='b', addCoef.col = 'black', number.cex=0.8, col=colorRampPalette(c("dodgerblue4","white","red4"))(200), col.lim=c(-1,1))
print(p)
dev.off()

## epithelium, malignant ---
wanted_celltypes <- c("Epithelial", "CAFs", "SMCs", "Cycling", "PNS_glial")
histo_feature <- c("GG3", "GG4", "GG4 Cribriform")

df_filt <- df_plot %>% dplyr::filter(celltype_major_v2 %in% all_of(wanted_celltypes))
df_sub <- df_filt %>% filter(Histology %in% all_of(histo_feature))
df_vessel <- df_sub %>% dplyr::select(-any_of(c("Barcode", "sample_id", "type", "Histology"))) %>% pivot_wider(names_from = "celltype_major_v2", values_from = "c2l_value") %>% remove_rownames %>% column_to_rownames(var="Barcode_sample_id")


plot_df <- cor(df_vessel)
# invert default colours
pdf(paste0(figDir, "epithelial_malignant_major.pdf"))
p <- corrplot(plot_df, title="malignant EPITHELIUM - all samples combined", mar=c(0,0,1,0), tl.cex = 0.95, tl.col='black', order='AOE', cl.pos='b', addCoef.col = 'black', number.cex=0.8, col=colorRampPalette(c("dodgerblue4","white","red4"))(200), col.lim=c(-1,1))
print(p)
dev.off()

## epithelium, benign ---
wanted_celltypes <- c("Epithelial", "CAFs", "SMCs", "Cycling", "PNS_glial")
histo_feature <- c("Epi Benign", "Epi Benign (transitional")

df_filt <- df_plot %>% dplyr::filter(celltype_major_v2 %in% all_of(wanted_celltypes))
df_sub <- df_filt %>% filter(Histology %in% all_of(histo_feature))
df_vessel <- df_sub %>% dplyr::select(-any_of(c("Barcode", "sample_id", "type", "Histology"))) %>% pivot_wider(names_from = "celltype_major_v2", values_from = "c2l_value") %>% remove_rownames %>% column_to_rownames(var="Barcode_sample_id")


plot_df <- cor(df_vessel)
# invert default colours
pdf(paste0(figDir, "epithelial_benign_major.pdf"))
p <- corrplot(plot_df, title="benign EPITHELIUM - all samples combined", mar=c(0,0,1,0), tl.cex = 0.95, tl.col='black', order='AOE', cl.pos='b', addCoef.col = 'black', number.cex=0.8, col=colorRampPalette(c("dodgerblue4","white","red4"))(200), col.lim=c(-1,1))
print(p)
dev.off()


## nerve ---
wanted_celltypes <- c("Epithelial", "CAFs", "SMCs", "Cycling", "PNS_glial")
histo_feature <- c("Nerve")

df_plot <- merged_df
df_plot$Barcode_sample_id <- paste0(df_plot$Barcode, "_", df_plot$sample_id)

df_filt <- df_plot %>% dplyr::filter(celltype_major_v2 %in% all_of(wanted_celltypes))
df_sub <- df_filt %>% filter(Histology %in% all_of(histo_feature))
df_vessel <- df_sub %>% dplyr::select(-any_of(c("Barcode", "sample_id", "type", "Histology"))) %>% pivot_wider(names_from = "celltype_major_v2", values_from = "c2l_value") %>% remove_rownames %>% column_to_rownames(var="Barcode_sample_id")


plot_df <- cor(df_vessel)
# invert default colours
pdf(paste0(figDir, "nerve_major.pdf"))
p <- corrplot(plot_df, title="NERVE - all samples combined", mar=c(0,0,1,0), tl.cex = 0.95, tl.col='black', order='AOE', cl.pos='b', addCoef.col = 'black', number.cex=0.8, col=colorRampPalette(c("dodgerblue4","white","red4"))(200), col.lim=c(-1,1))
print(p)
dev.off()

# SAVE LOGS -----------------------------------------------------------------
writeLines(capture.output(sessionInfo()), paste0(resultDir, "sessionInfo.txt"))