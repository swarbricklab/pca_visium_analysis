# SCRIPT INFORMATION -------------------------------------------------------

# author: Eva Apostolov
# date: 2024-03-15
# last edit: 2024-07-23

# script name: 08_visualise_proportions_minor.R
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
analysis = "08_visualise_proportions_minor/20240723"

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

# link to histo annotations (copied manually from Dropbox):
histo_path <- paste0(projectDir, "/", repo, "/data/20240318_reviewed_histo_annotation_FFPE/")

histo_df <- NULL

# Get the list of files and filter for those ending in C1, C2, or containing 20033
file_list <- list.files(histo_path, full.names = TRUE)
filtered_files <- file_list[grepl("C[12]_Histology_Reviewed\\.csv$|20033", file_list)]

for(file in unique(filtered_files)){
  print(file)

  temp_df <- read.csv(file)
  
  # Check if the file is one of the special cases
  if(grepl("PCa20130_C1_20272_C2|PCa20153_C1_20128_C1", file)) {
    # Use the sample_id from the file itself
    sample_id <- temp_df$sample_id
  } else {
    # Remove everything before one letter one digit pattern at the start (including the pattern) and everything after "_Histology_Reviewed.csv"
    sample_id <- sub(".*/([A-Za-z][0-9])_", "", file)
    sample_id <- sub("_Histology_Reviewed.csv", "", sample_id)
    
    # Swap _C1 for -1 and _C2 for -2
    sample_id <- sub("_C1", "-1", sample_id)
    sample_id <- sub("_C2", "-2", sample_id)
  }

  temp_df$sample_id <- sample_id

  histo_df <- rbind(histo_df, temp_df)
}


# path to config file used in cell2location - for sample ids --------------------
configFile = paste0(projectDir, repo, "/config/sample_sheet.csv")

samples_config <- read.csv(configFile)

# don't include Batch 1 samples, they are low QC
samples_config_sub <- samples_config %>% filter(batch != 1)

df <- NULL

for (section_name in unique(samples_config_sub$section_name)){

  print(section_name)

  # path to cell2location (Cell type major) - per sample
  prop_df <- read.csv(paste0("/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/02_cell2location/03a_cell_type_mapping/", section_name, "/objects/q05_cell_abundance_w_sf.csv"))
  prop_df$section_name <- section_name

  # Extracting the type corresponding to the sample_id
  core_type <- samples_config_sub$type[samples_config_sub$section_name == section_name]
  prop_df$type <- core_type
  
  # prop_df <- prop_df %>% remove_rownames %>% column_to_rownames(var="X")
  names(prop_df)[names(prop_df) == "X"] <- "Barcode"

  df <- rbind(df, prop_df)

}

# Given vector of column names to exclude
exclude_cols <- c("Barcode", "section_name", "type")

# Extract column names not found in exclude_cols
celltypes <- names(df)[!names(df) %in% exclude_cols]

df_tidy <- df %>% pivot_longer(celltypes, names_to = "celltype_minor_v2", values_to = "c2l_value")

      #  # use proportions moving forward, a bit more interpretable.. 
      #  df_tidy <- df_tidy %>%
      #    group_by(sample_id, Barcode) %>%
      #    mutate(prop_value = c2l_value / sum(c2l_value)) %>% select(-c2l_value)

# Count how many unique barcodes per sample_id
barcode_counts <- df_tidy %>%
  group_by(section_name) %>%
  summarise(unique_barcode_count = n_distinct(Barcode))

# combine the two dfs --------------------
# keeping only the 6 samples that have histo anno

# Merge based on both sample_id and Barcode
merged_df <- inner_join(df_tidy, histo_df, by = c("section_name", "Barcode"))

# exclude "Exclude" histology and blank cells from the analysis - this will be spots outside of the tissue, etc
merged_df <- merged_df %>% filter(Histology != "Exclude" & Histology != "")

# exclude "Exclude" in sample_id from the analysis
merged_df <- merged_df %>% filter(sample_id != "Exclude")

# in Histology column, swap empty spaces with an underscore
merged_df$Histology <- gsub(" ", "_", merged_df$Histology)

# replace Nerve/Ganglia with Nerve only
merged_df$Histology <- gsub("Nerve/Ganglia", "Nerve", merged_df$Histology)

# match cases/histology in GG4_Cribriform/Stroma prostate/Vessels
merged_df$Histology <- gsub("GG4_cribriform", "GG4_Cribriform", merged_df$Histology)
merged_df$Histology <- gsub("Stroma_prostate", "Stroma_prostatic", merged_df$Histology)
merged_df$Histology <- gsub("Vessels", "Vessel", merged_df$Histology)
merged_df$Histology <- gsub("Epi_Benign_\\(transitional\\)", "Epi_Benign_transitional", merged_df$Histology)

# plot some simple dfs --------------------

# MINOR ----------------------------------------------------------------------------

# 1 ------------
# plot histo annotations per core type, per sample_id 

# Generate a set of colors using a ColorBrewer palette
nice_colors_func <- colorRampPalette(brewer.pal(17, "Set1"))
nice_colors <- nice_colors_func(17)

df_plot <- merged_df %>%
  dplyr::group_by(type, Histology) %>%
  dplyr::summarise(histology_count = n_distinct(Barcode)) %>%
  dplyr::group_by(type) %>%
  dplyr::mutate(proportion = histology_count / sum(histology_count)) 

# Plot using ggplot2 with ColorBrewer colors
p <- ggplot(df_plot, aes(x = type, y = proportion, fill = Histology)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(title = "Histo Annotations per Core Type",
       x = "Histology",
       y = "Proportions") +
  scale_fill_manual(values = nice_colors) +  # Apply ColorBrewer colors
  theme_minimal()

pdf(paste0(figDir, "histo_anno_per_core_type.pdf"))
print(p)
dev.off()  

# show composition per sample
df_plot <- merged_df %>%
  dplyr::group_by(sample_id, Histology, type) %>%
  dplyr::summarise(histology_count = n_distinct(Barcode)) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(proportion = histology_count / sum(histology_count))

agg <- df_plot %>%
  filter(Histology == "Epi Benign") %>%
  group_by(sample_id) %>%
  summarise(freq = sum(proportion)) %>%
  arrange(desc(freq)) # arrange by decreasing frequency

## add missing levels and set the order of celltype_subset_v2
df_plot$sample_id <- factor(df_plot$sample_id, levels = agg$sample_id)

p <- ggplot(df_plot, aes(x = as.factor(sample_id), y = proportion, fill = Histology)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(title = "Histo Annotations per Sample",
       x = "Histology",
       y = "Proportions") +
  scale_fill_manual(values = nice_colors) +
  facet_wrap(~type, scales="free") +
  theme_minimal()

pdf(paste0(figDir, "histo_anno_per_sample.pdf"))
print(p)
dev.off()  

# 2 ------------
# start plotting cell2location results - all cell types

# proportion of each cell type within sample id?
celltype_prop_df <- merged_df %>%
  dplyr::group_by(sample_id, Barcode, celltype_minor_v2, type) %>%
  dplyr::summarise(c2l_cell_number = sum(c2l_value)) %>%
  dplyr::group_by(sample_id, Barcode) %>%
  dplyr::mutate(total_c2l_cell_number = sum(c2l_cell_number)) %>%
  dplyr::mutate(celltype_proportion = c2l_cell_number / total_c2l_cell_number) %>%
  dplyr::select(Barcode, sample_id, celltype_minor_v2, celltype_proportion, type) %>%
  dplyr::filter(celltype_proportion >= 0.05) # only include cell types that constitute at least 5% of a spot

# Generate a set of colors using a ColorBrewer palette
nice_colors_func <- colorRampPalette(brewer.pal(9, "Set1"))
nice_colors <- nice_colors_func(52)

df_plot <- celltype_prop_df

# plot per core type
p <- ggplot(df_plot, aes(x = type, y = celltype_proportion, fill = celltype_minor_v2)) +
  geom_bar(position = "fill", stat="identity") +
  labs(title = "Cell type minor composition (Cell2Location) per Core Type",
       subtitle = " Cell types constituting <5% of a spatial spot are filtered out",
       x = "Cell type minor (c2l decon)",
       y = "Proportions") +
  scale_fill_manual(values = nice_colors) +  # Apply ColorBrewer colors
  theme_minimal()

pdf(paste0(figDir, "cell2location_per_core_type.pdf"))
print(p)
dev.off()  

# plot per sample (based on above calculations)
agg <- df_plot %>%
  filter(celltype_minor_v2 == "Luminal") %>%
  group_by(sample_id) %>%
  summarise(freq = sum(celltype_proportion)) %>%
  arrange(desc(freq)) # arrange by decreasing frequency

## add missing levels and set the order of celltype_subset_v2
df_plot$sample_id <- factor(df_plot$sample_id, levels = agg$sample_id)

df_plot <- df_plot %>% group_by(sample_id) %>% arrange(desc(celltype_proportion))

p <- ggplot(df_plot, aes(x = as.factor(sample_id), y = celltype_proportion, fill = celltype_minor_v2)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(title = "Cell type minor composition (Cell2Location) per Sample",
       subtitle = " Cell types constituting <5% of a spatial spot are filtered out",
       x = "Cell type minor (c2l decon)",
       y = "Proportions") +
  scale_fill_manual(values = nice_colors) +
  facet_wrap(~type, scales="free") +
  theme_minimal()

pdf(paste0(figDir, "cell2location_anno_per_sample.pdf"), width=10)
print(p)
dev.off()  

# plot per histological region, per sample
celltype_prop_df <- merged_df %>%
  dplyr::group_by(sample_id, Histology, Barcode, celltype_minor_v2) %>%
  dplyr::summarise(c2l_cell_number = sum(c2l_value)) %>%
  dplyr::group_by(sample_id, Histology, Barcode) %>%
  dplyr::mutate(total_c2l_cell_number = sum(c2l_cell_number)) %>%
  dplyr::mutate(celltype_proportion = c2l_cell_number / total_c2l_cell_number) %>%
  dplyr::select(sample_id, Histology, Barcode, celltype_minor_v2, celltype_proportion) %>%
  dplyr::filter(celltype_proportion >= 0.05)

df_plot <- celltype_prop_df

# agg <- df_plot %>%
#  filter(celltype_minor_v2 == "Luminal") %>%
#  group_by(sample_id) %>%
#  summarise(freq = sum(celltype_proportion)) %>%
#  arrange(desc(freq)) # arrange by decreasing frequency
# 
# ## add missing levels and set the order of celltype_subset_v2
# df_plot$sample_id <- factor(df_plot$sample_id, levels = agg$sample_id)

p <- ggplot(df_plot, aes(x = Histology, y = celltype_proportion, fill = celltype_minor_v2)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(title = "Cell2Location decon (cell type minor) per Sample",
       subtitle = " Cell types constituting <5% of a spatial spot are filtered out",
       x = "Cell type minor (c2l decon)",
       y = "Proportions") +
  scale_fill_manual(values = nice_colors) +
  facet_wrap(~as.factor(sample_id), scales="free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0(figDir, "cell2location_anno_per_histo_region_per_sample.pdf"), width=12)
print(p)
dev.off()  


# now plot populations of interest ----
# Luminal, CAFs, SMCs, and Glial

# group glial cells into PNS_glial & sum the c2l values - TO DO: this needs to be redone!
merged_df_v2 <- merged_df %>%
  mutate(celltype_minor_v2 = case_when(
    grepl("Glial_", celltype_minor_v2) ~ "PNS_glial",
    grepl("Satellite", celltype_minor_v2) ~ "PNS_glial",
    grepl("Schwann", celltype_minor_v2) ~ "PNS_glial",
    TRUE ~ celltype_minor_v2
  )) %>%
  group_by(Barcode, sample_id, celltype_minor_v2, type, Histology) %>%
  summarise(c2l_value = sum(c2l_value)) %>%
  ungroup()

celltype_prop_df <- merged_df_v2 %>%
  dplyr::group_by(sample_id, Histology, Barcode, celltype_minor_v2) %>%
  dplyr::summarise(c2l_cell_number = sum(c2l_value)) %>%
  dplyr::group_by(sample_id, Histology, Barcode) %>%
  dplyr::mutate(total_c2l_cell_number = sum(c2l_cell_number)) %>%
  dplyr::mutate(celltype_proportion = c2l_cell_number / total_c2l_cell_number) %>%
  dplyr::select(sample_id, Histology, Barcode, celltype_minor_v2, celltype_proportion) %>%
  dplyr::filter(celltype_proportion >= 0.05) %>% # only include cell type proportions with a value greater than 5%
  # dplyr::filter(celltype_minor_v2 %in% c("Luminal", "uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "pSMCs", "Pericytes", "vSMCs", "Cycling_Epi_Luminal", "Glial_cells_1", "Glial_cells_2", "Glial_cells_3", "Satellite_cells", "Myelinating_Schwann"))
  dplyr::filter(celltype_minor_v2 %in% c("Luminal", "uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "pSMCs", "Pericytes", "vSMCs", "Cycling_Epi_Luminal", "PNS_glial"))

df_plot <- celltype_prop_df

nice_colors_func <- colorRampPalette(brewer.pal(9, "Set1"))
nice_colors <- nice_colors_func(16)

p <- ggplot(df_plot, aes(x = Histology, y = celltype_proportion, fill = celltype_minor_v2)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(title = "Cell2Location decon (cell type minor) per Histological region, per Sample",
       subtitle = " Cell types constituting <5% of a spatial spot are filtered out; showing only cell types of interest",
       x = "Cell type minor (c2l decon)",
       y = "Proportions") +
  scale_fill_manual(values = nice_colors) +
  facet_wrap(~as.factor(sample_id), scales="free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0(figDir, "cell2location_anno_per_histo_region_per_sample_of_interest_5_percent.pdf"), width=9)
print(p)
dev.off()  


# 3 ------------
# plots for talk

# PART 1
## show proportions of Luminal, CAFs, SMCs, and Glial per core type

celltype_prop_df <- merged_df_v2 %>%
  dplyr::group_by(Histology, Barcode, celltype_minor_v2, type) %>%
  dplyr::summarise(c2l_cell_number = sum(c2l_value)) %>%
  dplyr::group_by(type, Barcode) %>%
  dplyr::mutate(total_c2l_cell_number = sum(c2l_cell_number)) %>%
  dplyr::mutate(celltype_proportion = c2l_cell_number / total_c2l_cell_number) %>%
  dplyr::select(type, Barcode, celltype_minor_v2, celltype_proportion, Histology) %>%
  dplyr::filter(celltype_proportion >= 0.05) %>% # only include cell type proportions with a value greater than 5%
  # dplyr::filter(celltype_minor_v2 %in% c("Luminal", "uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "pSMCs", "Pericytes", "vSMCs", "Cycling_Epi_Luminal", "Glial_cells_1", "Glial_cells_2", "Glial_cells_3", "Satellite_cells", "Myelinating_Schwann"))
  # dplyr::filter(celltype_minor_v2 %in% c("Luminal", "uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "pSMCs", "Pericytes", "vSMCs", "Cycling_Epi_Luminal", "PNS_glial"))
  dplyr::filter(celltype_minor_v2 %in% c("Luminal", "Cycling_Epi_Luminal", "uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "pSMCs", "Pericytes", "vSMCs", "PNS_glial"))

df_plot <- celltype_prop_df

nice_colors_func <- colorRampPalette(brewer.pal(9, "Set1"))
nice_colors <- nice_colors_func(12)

p <- ggplot(df_plot, aes(x = type, y = celltype_proportion, fill = celltype_minor_v2)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(title = "Cell2Location decon (cell type minor) per Histological region, per Sample",
       subtitle = " Cell types constituting <5% of a spatial spot are filtered out; showing only cell types of interest",
       x = "Cell type minor (c2l decon)",
       y = "Proportions") +
  scale_fill_manual(values = nice_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0(figDir, "Epithelial_included_C2L_per_core_type_5_percent.pdf"), width=5)
print(p)
dev.off()

# PART 2
## show proportions of CAFs, SMCs, and Glial per core type

df_plot <- celltype_prop_df %>% dplyr::filter(celltype_minor_v2 %in% c("uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "PNS_glial", "pSMCs", "Pericytes", "vSMCs"))

nice_colors_func <- colorRampPalette(brewer.pal(9, "Set1"))
nice_colors <- nice_colors_func(11)

p <- ggplot(df_plot, aes(x = type, y = celltype_proportion, fill = celltype_minor_v2)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(title = "Cell2Location decon (cell type minor) per Histological region, per Sample",
       subtitle = " Cell types constituting <5% of a spatial spot are filtered out; showing only cell types of interest",
       x = "Cell type minor (c2l decon)",
       y = "Proportions") +
  scale_fill_manual(values = nice_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0(figDir, "CAFs_SMCs_Glial_C2L_per_core_type_5_percent.pdf"), width=5)
print(p)
dev.off()

## show proportions of CAFs, SMCs, and Glial per histological_region
celltype_prop_df <- merged_df_v2 %>%
  dplyr::group_by(sample_id, Histology, Barcode, celltype_minor_v2) %>%
  dplyr::summarise(c2l_cell_number = sum(c2l_value)) %>%
  dplyr::group_by(sample_id, Histology, Barcode) %>%
  dplyr::mutate(total_c2l_cell_number = sum(c2l_cell_number)) %>%
  dplyr::mutate(celltype_proportion = c2l_cell_number / total_c2l_cell_number) %>%
  dplyr::select(sample_id, Histology, Barcode, celltype_minor_v2, celltype_proportion) %>%
  dplyr::filter(celltype_proportion >= 0.05) %>% # only include cell type proportions with a value greater than 20%
  # dplyr::filter(celltype_minor_v2 %in% c("Luminal", "uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "pSMCs", "Pericytes", "vSMCs", "Cycling_Epi_Luminal", "Glial_cells_1", "Glial_cells_2", "Glial_cells_3", "Satellite_cells", "Myelinating_Schwann"))
  dplyr::filter(celltype_minor_v2 %in% c("uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "PNS_glial", "pSMCs", "Pericytes", "vSMCs"))

df_plot <- celltype_prop_df

nice_colors_func <- colorRampPalette(brewer.pal(9, "Set1"))
nice_colors <- nice_colors_func(16)

p <- ggplot(df_plot, aes(x = Histology, y = celltype_proportion, fill = celltype_minor_v2)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(title = "Cell2Location decon (cell type minor) per Histological region, per Sample",
       subtitle = " Cell types constituting <5% of a spatial spot are filtered out; showing only cell types of interest",
       x = "Cell type minor (c2l decon)",
       y = "Proportions") +
  scale_fill_manual(values = nice_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0(figDir, "CAFs_SMCs_Glial_C2L_histo_regions_5_percent.pdf"), width=9)
print(p)
dev.off()  


## show proportions of CAFs, SMCs, and Glial per histological_region, facet by core type..
celltype_prop_df <- merged_df_v2 %>%
  dplyr::group_by(sample_id, Histology, Barcode, celltype_minor_v2, type) %>%
  dplyr::summarise(c2l_cell_number = sum(c2l_value)) %>%
  dplyr::group_by(sample_id, Histology, Barcode) %>%
  dplyr::mutate(total_c2l_cell_number = sum(c2l_cell_number)) %>%
  dplyr::mutate(celltype_proportion = c2l_cell_number / total_c2l_cell_number) %>%
  dplyr::select(sample_id, Histology, Barcode, celltype_minor_v2, celltype_proportion, type) %>%
  dplyr::filter(celltype_proportion >= 0.05) %>% # only include cell type proportions with a value greater than 5%
  # dplyr::filter(celltype_minor_v2 %in% c("Luminal", "uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "pSMCs", "Pericytes", "vSMCs", "Cycling_Epi_Luminal", "Glial_cells_1", "Glial_cells_2", "Glial_cells_3", "Satellite_cells", "Myelinating_Schwann"))
  dplyr::filter(celltype_minor_v2 %in% c("uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "PNS_glial", "pSMCs", "Pericytes", "vSMCs"))

df_plot <- celltype_prop_df

nice_colors_func <- colorRampPalette(brewer.pal(9, "Set1"))
nice_colors <- nice_colors_func(16)

df_plot$Histology <- factor(df_plot$Histology, levels = c("Epi_Benign", "Epi_Benign_transitional",  "GG3", "GG4",  "GG4_Cribriform", "Inflammation", "Stroma_prostatic", "Stroma_extraprostatic",  "Nerve", "Vessel", "Adipose"))

p <- ggplot(df_plot, aes(x = type, y = celltype_proportion, fill = celltype_minor_v2)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(title = "Cell2Location decon (cell type minor) per Histological region, per Sample",
       subtitle = " Cell types constituting <5% of a spatial spot are filtered out; showing only cell types of interest",
       x = "Cell type minor (c2l decon)",
       y = "Proportions") +
  facet_wrap(~as.factor(Histology), ncol = 11) +
  scale_fill_manual(values = nice_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0(figDir, "CAFs_SMCs_Glial_C2L_histo_regions_5_percent_type_facet.pdf"), width=18)
print(p)
dev.off()  

# PART 3
## show proportions of CAFs and SMCs per core type ----------------
df_plot <- celltype_prop_df %>% dplyr::filter(celltype_minor_v2 %in% c("uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "pSMCs", "Pericytes", "vSMCs"))

nice_colors_func <- colorRampPalette(brewer.pal(9, "Set1"))
nice_colors <- nice_colors_func(11)

p <- ggplot(df_plot, aes(x = type, y = celltype_proportion, fill = celltype_minor_v2)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(title = "Cell2Location decon (cell type minor) per Histological region, per Sample",
       subtitle = " Cell types constituting <5% of a spatial spot are filtered out; showing only cell types of interest",
       x = "Cell type minor (c2l decon)",
       y = "Proportions") +
  scale_fill_manual(values = nice_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0(figDir, "CAFs_SMCs_C2L_per_core_type_5_percent.pdf"), width=5)
print(p)
dev.off()

## show proportions of CAFs and SMCs per histological_region
celltype_prop_df <- merged_df_v2 %>%
  dplyr::group_by(sample_id, Histology, Barcode, celltype_minor_v2) %>%
  dplyr::summarise(c2l_cell_number = sum(c2l_value)) %>%
  dplyr::group_by(sample_id, Histology, Barcode) %>%
  dplyr::mutate(total_c2l_cell_number = sum(c2l_cell_number)) %>%
  dplyr::mutate(celltype_proportion = c2l_cell_number / total_c2l_cell_number) %>%
  dplyr::select(sample_id, Histology, Barcode, celltype_minor_v2, celltype_proportion) %>%
  dplyr::filter(celltype_proportion >= 0.05) %>% # only include cell type proportions with a value greater than 20%
  # dplyr::filter(celltype_minor_v2 %in% c("Luminal", "uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "pSMCs", "Pericytes", "vSMCs", "Cycling_Epi_Luminal", "Glial_cells_1", "Glial_cells_2", "Glial_cells_3", "Satellite_cells", "Myelinating_Schwann"))
  dplyr::filter(celltype_minor_v2 %in% c("uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "pSMCs", "Pericytes", "vSMCs"))

df_plot <- celltype_prop_df

nice_colors_func <- colorRampPalette(brewer.pal(9, "Set1"))
nice_colors <- nice_colors_func(16)

p <- ggplot(df_plot, aes(x = Histology, y = celltype_proportion, fill = celltype_minor_v2)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(title = "Cell2Location decon (cell type minor) per Histological region, per Sample",
       subtitle = " Cell types constituting <5% of a spatial spot are filtered out; showing only cell types of interest",
       x = "Cell type minor (c2l decon)",
       y = "Proportions") +
  scale_fill_manual(values = nice_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0(figDir, "CAFs_SMCs_C2L_histo_regions_5_percent.pdf"), width=9)
print(p)
dev.off()  


## show proportions of CAFs and SMCs per histological_region, facet by core type..
celltype_prop_df <- merged_df_v2 %>%
  dplyr::group_by(sample_id, Histology, Barcode, celltype_minor_v2, type) %>%
  dplyr::summarise(c2l_cell_number = sum(c2l_value)) %>%
  dplyr::group_by(sample_id, Histology, Barcode) %>%
  dplyr::mutate(total_c2l_cell_number = sum(c2l_cell_number)) %>%
  dplyr::mutate(celltype_proportion = c2l_cell_number / total_c2l_cell_number) %>%
  dplyr::select(sample_id, Histology, Barcode, celltype_minor_v2, celltype_proportion, type) %>%
  dplyr::filter(celltype_proportion >= 0.05) %>% # only include cell type proportions with a value greater than 5%
  # dplyr::filter(celltype_minor_v2 %in% c("Luminal", "uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "pSMCs", "Pericytes", "vSMCs", "Cycling_Epi_Luminal", "Glial_cells_1", "Glial_cells_2", "Glial_cells_3", "Satellite_cells", "Myelinating_Schwann"))
  dplyr::filter(celltype_minor_v2 %in% c("uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "pSMCs", "Pericytes", "vSMCs"))

df_plot <- celltype_prop_df

nice_colors_func <- colorRampPalette(brewer.pal(9, "Set1"))
nice_colors <- nice_colors_func(16)

df_plot$Histology <- factor(df_plot$Histology, levels = c("Epi_Benign", "Epi_Benign_transitional",  "GG3", "GG4",  "GG4_Cribriform", "Inflammation", "Stroma_prostatic", "Stroma_extraprostatic",  "Nerve", "Vessel", "Adipose"))

p <- ggplot(df_plot, aes(x = type, y = celltype_proportion, fill = celltype_minor_v2)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(title = "Cell2Location decon (cell type minor) per Histological region, per Sample",
       subtitle = " Cell types constituting <5% of a spatial spot are filtered out; showing only cell types of interest",
       x = "Cell type minor (c2l decon)",
       y = "Proportions") +
  facet_wrap(~as.factor(Histology)) +
  scale_fill_manual(values = nice_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0(figDir, "CAFs_SMCs_C2L_histo_regions_5_percent_type_facet.pdf"), width=18)
print(p)
dev.off()  


## calculate proportions per sample (based on barspot proportions)
# sample_celltype_prop_df <- t %>%
#   group_by(sample_id, celltype_minor_v2, type) %>%
#   summarise(celltype_proportion = sum(celltype_proportion)) %>%
#   group_by(sample_id) %>%
#   mutate(total_proportion = sum(celltype_proportion)) %>%
#   mutate(celltype_proportion = celltype_proportion / total_proportion) %>%
#   ungroup() %>%
#   select(sample_id, celltype_minor_v2, celltype_proportion, type)

# calculate spatial correlations ------------------------------
df_plot <- merged_df 
df_plot$Barcode_sample_id <- paste0(df_plot$Barcode, "_", df_plot$sample_id)


## stroma ---
wanted_celltypes <- c("Luminal", "uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "pSMCs", "Pericytes", "vSMCs", "Cycling_Epi_Luminal", "PNS_glial")
histo_feature <- c("Stroma_prostatic", "Stroma_extraprostatic")

df_filt <- df_plot %>% dplyr::filter(celltype_minor_v2 %in% all_of(wanted_celltypes))
df_sub <- df_filt %>% filter(Histology %in% all_of(histo_feature))
df_stroma <- df_sub %>% dplyr::select(-any_of(c("Barcode", "sample_id", "type", "Histology", "section_name"))) %>% pivot_wider(names_from = "celltype_minor_v2", values_from = "c2l_value") %>% remove_rownames %>% column_to_rownames(var="Barcode_sample_id")


plot_df <- cor(df_stroma)
# invert default colours
pdf(paste0(figDir, "stroma_minor.pdf"))
p <- corrplot(plot_df, title="STROMA - all samples combined", mar=c(0,0,1,0), tl.cex = 0.95, tl.col='black', order='AOE', cl.pos='b', addCoef.col = 'black', number.cex=0.8, col=colorRampPalette(c("dodgerblue4","white","red4"))(200), col.lim=c(-1,1))
print(p)
dev.off()


## vessel ---
wanted_celltypes <- c("Luminal", "uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "pSMCs", "Pericytes", "vSMCs", "Cycling_Epi_Luminal", "PNS_glial")
histo_feature <- c("Vessel")

df_filt <- df_plot %>% dplyr::filter(celltype_minor_v2 %in% all_of(wanted_celltypes))
df_sub <- df_filt %>% filter(Histology %in% all_of(histo_feature))
df_vessel <- df_sub %>% dplyr::select(-any_of(c("Barcode", "sample_id", "type", "Histology", "section_name"))) %>% pivot_wider(names_from = "celltype_minor_v2", values_from = "c2l_value") %>% remove_rownames %>% column_to_rownames(var="Barcode_sample_id")


plot_df <- cor(df_vessel)
# invert default colours
pdf(paste0(figDir, "vessel_minor.pdf"))
p <- corrplot(plot_df, title="VESSEL - all samples combined", mar=c(0,0,1,0), tl.cex = 0.95, tl.col='black', order='AOE', cl.pos='b', addCoef.col = 'black', number.cex=0.8, col=colorRampPalette(c("dodgerblue4","white","red4"))(200), col.lim=c(-1,1))
print(p)
dev.off()

## epithelium, malignant ---
wanted_celltypes <- c("Luminal", "uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "pSMCs", "Pericytes", "vSMCs", "Cycling_Epi_Luminal", "PNS_glial")
histo_feature <- c("GG3", "GG4", "GG4_Cribriform")

df_filt <- df_plot %>% dplyr::filter(celltype_minor_v2 %in% all_of(wanted_celltypes))
df_sub <- df_filt %>% filter(Histology %in% all_of(histo_feature))
df_vessel <- df_sub %>% dplyr::select(-any_of(c("Barcode", "sample_id", "type", "Histology", "section_name"))) %>% pivot_wider(names_from = "celltype_minor_v2", values_from = "c2l_value") %>% remove_rownames %>% column_to_rownames(var="Barcode_sample_id")


plot_df <- cor(df_vessel)
# invert default colours
pdf(paste0(figDir, "epithelial_malignant_minor.pdf"))
p <- corrplot(plot_df, title="malignant EPITHELIUM - all samples combined", mar=c(0,0,1,0), tl.cex = 0.95, tl.col='black', order='AOE', cl.pos='b', addCoef.col = 'black', number.cex=0.8, col=colorRampPalette(c("dodgerblue4","white","red4"))(200), col.lim=c(-1,1))
print(p)
dev.off()

## epithelium, benign ---
wanted_celltypes <- c("Luminal", "uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "pSMCs", "Pericytes", "vSMCs", "Cycling_Epi_Luminal", "PNS_glial")
histo_feature <- c("Epi_Benign", "Epi_Benign_transitional")

df_filt <- df_plot %>% dplyr::filter(celltype_minor_v2 %in% all_of(wanted_celltypes))
df_sub <- df_filt %>% filter(Histology %in% all_of(histo_feature))
df_vessel <- df_sub %>% dplyr::select(-any_of(c("Barcode", "sample_id", "type", "Histology", "section_name"))) %>% pivot_wider(names_from = "celltype_minor_v2", values_from = "c2l_value") %>% remove_rownames %>% column_to_rownames(var="Barcode_sample_id")


plot_df <- cor(df_vessel)
# invert default colours
pdf(paste0(figDir, "epithelial_benign_minor.pdf"))
p <- corrplot(plot_df, title="benign EPITHELIUM - all samples combined", mar=c(0,0,1,0), tl.cex = 0.95, tl.col='black', order='AOE', cl.pos='b', addCoef.col = 'black', number.cex=0.8, col=colorRampPalette(c("dodgerblue4","white","red4"))(200), col.lim=c(-1,1))
print(p)
dev.off()


## nerve ---
wanted_celltypes <- c("Luminal", "uniCAFs", "NPF_like", "pnCAFs", "whCAFs", "SMC_like", "CAFs_IFN", "pSMCs", "Pericytes", "vSMCs", "Cycling_Epi_Luminal", "PNS_glial")
histo_feature <- c("Nerve")

df_plot <- merged_df_v2
df_plot$Barcode_sample_id <- paste0(df_plot$Barcode, "_", df_plot$sample_id)

df_filt <- df_plot %>% dplyr::filter(celltype_minor_v2 %in% all_of(wanted_celltypes))
df_sub <- df_filt %>% filter(Histology %in% all_of(histo_feature))
df_vessel <- df_sub %>% dplyr::select(-any_of(c("Barcode", "sample_id", "type", "Histology", "section_name"))) %>% pivot_wider(names_from = "celltype_minor_v2", values_from = "c2l_value") %>% remove_rownames %>% column_to_rownames(var="Barcode_sample_id")


plot_df <- cor(df_vessel)
# invert default colours
pdf(paste0(figDir, "nerve_minor.pdf"))
p <- corrplot(plot_df, title="NERVE - all samples combined", mar=c(0,0,1,0), tl.cex = 0.95, tl.col='black', order='AOE', cl.pos='b', addCoef.col = 'black', number.cex=0.8, col=colorRampPalette(c("dodgerblue4","white","red4"))(200), col.lim=c(-1,1))
print(p)
dev.off()

# SAVE LOGS -----------------------------------------------------------------
writeLines(capture.output(sessionInfo()), paste0(resultDir, "sessionInfo.txt"))