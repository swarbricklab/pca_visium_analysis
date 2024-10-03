# spatial co-localisation analysis, based on cell2location output
# across entire section, cell-to-cell

# conda env: atac_2022
# date: 2024-04-15
# last edit: 2024-07-22

# 01: SETUP---------------------------------------

# setup
library(ggplot2)
library(magrittr)
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(tidyr)
library(tibble)
library(pheatmap)

# for stats
library(rstatix)

# set up directory structure
projectName="PCa_Visium"
projectDir=paste0("/share/ScratchGeneral/evaapo/projects/",projectName,"/")
repo="pca_visium_analysis"

# rerun directory - trying to organise PCa analysis
exp = "04_colocalisation_inhouse_and_published"
analysis = "04_colocalisation_scanpy_scores"

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


type_cols_darker <- c(
  "cancer" = "#9065cf",
  "adj_benign" = "#E0B57C") ###E68568"

# FUNCTIONS ------------------------------------------------------------------

# add pearson correlations to heatmap
# need to work out if i need multiple testing correction

# Function to calculate Pearson correlation and p-values for a given sample_id
calculate_correlations <- function(df, sample_id, celltype) {
  # Filter the data for the given sample_id
  sample_data <- df %>%
    filter(sample_id == !!sample_id) %>% select(-sample_id, -Barcode, -type, -Histology, -section_name, -has_vessels, -nerve_enriched)

  # Select the column of interest for correlation analysis (specified celltype)
  numeric_vector <- sample_data[[celltype]]

  # Initialize an empty data frame to store results
  result_df <- data.frame(
    Celltype_Pair = character(),  # Column to identify the signature pair
    Correlation = numeric(),       # Column for correlation coefficients
    P_Value = numeric()            # Column for p-values
  )

  # Loop through each column (signature) for correlation analysis
  for (signature_col in colnames(sample_data)) {
    # Skip the column representing the specified celltype itself
    if (signature_col != celltype) {
      # Run correlation test
      cor_pval_result <- cor.test(numeric_vector, sample_data[[signature_col]], method = "pearson")

      # Extract results
      correlation_coefficient <- cor_pval_result$estimate
      p_value <- cor_pval_result$p.value

      # Append results to the data frame
      result_df <- rbind(result_df, data.frame(
        Celltype_Pair = paste(celltype, "_vs_", signature_col, sep = ""),
        Correlation = correlation_coefficient,
        P_Value = p_value
      ))
    }
  }

  return(result_df)
}


# Function to calculate Pearson correlation and p-values for a given sample_id
calculate_correlations_histo <- function(df, histology, celltype) {
  # Filter the data for the given sample_id
  sample_data <- df %>%
    filter(Histology == !!histology) %>% select(-sample_id, -Barcode, -type, -Histology, -section_name, -has_vessels, -nerve_enriched)

  # Select the column of interest for correlation analysis (specified celltype)
  numeric_vector <- sample_data[[celltype]]

  # Initialize an empty data frame to store results
  result_df <- data.frame(
    Celltype_Pair = character(),  # Column to identify the signature pair
    Correlation = numeric(),       # Column for correlation coefficients
    P_Value = numeric()            # Column for p-values
  )

  # Loop through each column (signature) for correlation analysis
  for (signature_col in colnames(sample_data)) {
    # Skip the column representing the specified celltype itself
    if (signature_col != celltype) {
      # Run correlation test
      cor_pval_result <- cor.test(numeric_vector, sample_data[[signature_col]], method = "pearson")

      # Extract results
      correlation_coefficient <- cor_pval_result$estimate
      p_value <- cor_pval_result$p.value

      # Append results to the data frame
      result_df <- rbind(result_df, data.frame(
        Celltype_Pair = paste(celltype, "_vs_", signature_col, sep = ""),
        Correlation = correlation_coefficient,
        P_Value = p_value
      ))
    }
  }

  return(result_df)
}


# Function to replace p-values with asterisks
replace_with_asterisks <- function(p_value) {
  if (p_value < 0.0001) {
    return("***")
  } else if (p_value < 0.001) {
    return("**")
  } else if (p_value < 0.01) {
    return("*")
  } else {
    return("")
  }
}

# START ----------------------------------------------------------------------

# read in histopath data from in-house data as a df ----------

# read in histopath data as a df ----------
# this is a bit convoluted to account for differences in naming conventions
# and PCa..samples containing an extra sample_id column in Loupe annotation

# link to histo annotations (copied manually from Dropbox):
histo_path <- paste0(projectDir, "/", repo, "/data/20240318_reviewed_histo_annotation_FFPE/")

histo_df1 <- NULL

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

  histo_df1 <- rbind(histo_df1, temp_df)
}

# Leave this out for now
# # read in histopath data for Mengxiao's dataset as a df --------------
# histo_path <- paste0(projectDir, "/", repo, "/resources/published_data/PMID_35948708/Count_matrices/Patient_1/Visium_with_annotation")
# 
# # list all files in subdirectories
# all_files <- list.files(histo_path, recursive = TRUE, full.names = TRUE)
# histo_path_files <- grep("Final_Consensus_Annotations", all_files, value=TRUE)
# 
# histo_df2 <- NULL
# 
# for(file in unique(histo_path_files)){
#   print(file)
# 
#   temp_df <- read.csv(file)
# 
#   # rename column (new = old)
#   temp_df <- temp_df %>% rename(Histology = Final_Annotations)
# 
#   # Remove everything before one letter one digit pattern at the start (including the pattern) and everything after "_Histology_Reviewed.csv"
#   sample_id <- sub(".*Visium_with_annotation/", "", file)
#   sample_id <- sub("/.*", "", sample_id)
# 
#   temp_df$sample_id <- sample_id
# 
#   histo_df2 <- rbind(histo_df2, temp_df)
# 
# }
# 
# # combine the two dfs
# histo_df <- rbind(histo_df1, histo_df2)

histo_df <- histo_df1

# path to config file used in cell2location - for sample ids --------------------
configFile = paste0(projectDir, repo, "/config/sample_sheet.csv")

samples_config <- read.csv(configFile)

# only include sample ids found in histo_df
# samples_config_sub <- samples_config %>% dplyr::filter(batch == 2 | batch == "mengxiao")
# samples_config_sub <- samples_config %>% dplyr::filter(sample_id %in% all_of(unique(histo_df$sample_id)))
samples_config_sub <- samples_config %>% filter(batch != 1)

# now, list all .csv files in results
# csv_files <- list.files("/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/marker_genes/per_sample/20240415", pattern = ".csv", recursive = TRUE, full.names = TRUE)
csv_files <- list.files("/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/marker_genes/per_sample/20240723", pattern = ".csv", recursive = TRUE, full.names = TRUE)
csv_files <- grep(paste(unique(samples_config_sub$section_name), collapse = "|"), csv_files, value = TRUE)

df <- NULL

for(file in csv_files){

  print(file)

  prop_df <- read.csv(file)

  section_name <- unique(prop_df[['section_name']])

  # Extracting the type corresponding to the sample_id
  core_type <- samples_config_sub$type[samples_config_sub$section_name == section_name]
  prop_df$type <- core_type
  
  names(prop_df)[names(prop_df) == "barcode_id"] <- "Barcode"

  df <- rbind(df, prop_df)

}

# Given vector of column names to exclude
exclude_cols <- c("Barcode", "section_name", "type")

# Extract column names not found in exclude_cols
celltypes <- names(df)[!names(df) %in% exclude_cols]

df_tidy <- df %>% pivot_longer(celltypes, names_to = "score_name", values_to = "score_value")

# Count how many unique barcodes per sample_id
barcode_counts <- df_tidy %>%
  group_by(section_name) %>%
  summarise(unique_barcode_count = n_distinct(Barcode))

# combine the two dfs --------------------
# keeping only the 6 samples that have histo anno

# Merge based on both sample_id and Barcode
merged_df <- inner_join(df_tidy, histo_df, by = c("section_name", "Barcode"))

spread_df <- spread(merged_df, key = score_name, value = score_value)

# > dim(spread_df)
# [1] 20143     8

# exclude "Exclude" histology and blank cells from the analysis - this will be spots outside of the tissue, etc
spread_df <- spread_df %>% filter(Histology != "Exclude" & Histology != "")

# > dim(spread_df)
# [1] 16977    8

# exclude "Exclude" in sample_id from the analysis
spread_df <- spread_df %>% filter(sample_id != "Exclude")

# > dim(spread_df)
# [1] 16976    8


# get annotation info from spread_df
# create the has_vessels and has_nerves columns based on the presence of "Vessel" and "Nerve" in the Histology column
spread_df2 <- spread_df %>%
  group_by(sample_id) %>%
  mutate(
    has_vessels = if_else(any(grepl("Vessel", Histology)), "yes", "no"),
    nerve_enriched = if_else(any(grepl("Nerve", Histology)), "yes", "no")
  ) %>%
  ungroup() %>% filter(nerve_enriched == "yes")

# ---------------------------------------------------
# plot scatter plots - per sample id ----------------
# ---------------------------------------------------

# _PNS_glial_pnCAF_scores_no_Glial_exp_genes_addDCN

cell_types_of_interest <- c('pnCAFs_score', 'PNS_glial_score', 'NPF_score')

plot_list <- list()

for(sample in unique(spread_df2$sample_id)){
  print(sample)

  df_subset <- spread_df2 %>% dplyr::filter(sample == sample_id) %>% select('Barcode', 'pnCAFs_score', 'PNS_glial_score')
  df_subset <- df_subset %>% remove_rownames() %>% column_to_rownames('Barcode')

  sp <- ggscatter(df_subset, x='pnCAFs_score', y='PNS_glial_score',
     add = "reg.line",  # Add regressin line
     add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
     conf.int = TRUE # Add confidence interval
     )

  # Add correlation coefficient
  calculate_label_positions <- function(df_subset) {
    # Calculate label positions based on data range
    label_x <- min(df_subset$pnCAFs_score) + 0.2 * diff(range(df_subset$pnCAFs_score))
    label_y <- min(df_subset$PNS_glial_score) + 0.5 * diff(range(df_subset$PNS_glial_score))
    return(list(label_x = label_x, label_y = label_y))
  }

  label_positions <- calculate_label_positions(df_subset)

  # Add correlation coefficient
  plot_list[[sample]] <- sp + stat_cor(method = "pearson", label.x = label_positions$label_x, label.y = label_positions$label_y, color="red") + labs(subtitle = sample)

}


pdf(paste0(figDir, "per_sample_scatter_plot_without_DCN_nerve_enriched.pdf"), width=14, height=10)
p <- cowplot::plot_grid(plotlist = plot_list)
print(p)
dev.off()

# ---------------------------------------------------
# test correlations, with significance --------------
# ---------------------------------------------------

# List of cell types of interest (lineage 1 comparison to all other cell types at minor level)
# interactions_df <- read.csv("/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/config/minor_mal_interactions_of_interest_v2.csv")
cell_types_of_interest <- c('pnCAFs_score', 'PNS_glial_score', 'NPF_score')

# Create an empty list to store the correlation data frames for each cell type
correlation_list <- list()

# Loop through each cell type of interest
for (cell_type in cell_types_of_interest) {
  print(cell_type)
  
  # Loop through each unique sample_id and calculate the correlations for the current cell type
  for (sample_id in unique(as.factor(spread_df2$sample_id))) {
    print(sample_id)
    
    # Calculate correlations for the current cell type and sample ID
    correlation_df <- calculate_correlations(spread_df2, sample_id, cell_type)

    # Add cell type and sample_id as columns
    correlation_df$sample_id <- sample_id
    correlation_df$cell_type <- cell_type
    
    # Store the correlation dataframe in the list
    correlation_list[[paste(cell_type, sample_id, sep = "_")]] <- correlation_df
  }
}

# Combine all correlation data frames into a single data frame
combined_cor_df <- do.call(rbind, correlation_list)

# Extract the interaction patterns from interactions_df
# interaction_patterns <- paste(interactions_df$celltype_1, "vs", interactions_df$celltype_2, sep = "_")

# Filter combined_cor_df to keep only the interactions found in interaction_patterns
# filtered_combined_cor_df <- combined_cor_df[combined_cor_df$Celltype_Pair %in% interaction_patterns, ]

# Prepare the correlation matrix
cor_matrix <- combined_cor_df %>% dplyr::select(-P_Value)

# Assuming your correlation matrix is named cor_matrix
cor_matrix <- cor_matrix %>%
  select(-cell_type) %>%
  pivot_wider(names_from = Celltype_Pair, values_from = Correlation)   %>%
  tibble::remove_rownames() %>% column_to_rownames(var = "sample_id")

# Prepare the p-value matrix
# p_value_matrix <- filtered_combined_cor_df %>% dplyr::select(Celltype_Pair, P_Value, sample_id) # %>% tibble::remove_rownames() %>% column_to_rownames(var = "sample_id")
p_value_matrix <- combined_cor_df %>% dplyr::select(Celltype_Pair, P_Value, sample_id) # %>% tibble::remove_rownames() %>% column_to_rownames(var = "sample_id")

# Assuming your correlation matrix is named cor_matrix
p_value_matrix_wide <- p_value_matrix %>%
  pivot_wider(names_from = Celltype_Pair, values_from = P_Value) %>%
  tibble::remove_rownames() %>% column_to_rownames(var = "sample_id")

# # Count the total number of tests
# m <- 30*6
# 
# # Apply Bonferroni correction to each element in the matrix
# corrected_matrix <- p_value_matrix_wide * m
# 
# # Ensure that the corrected p-values are capped at 1
# corrected_matrix[corrected_matrix > 1] <- 1

# or adjust with p.adjust(): https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust - TO DO

# Apply the function to the p-value matrix
asterisk_matrix <- apply(p_value_matrix_wide, 2, function(column) sapply(column, replace_with_asterisks))

t_cor_matrix <- t(cor_matrix)
t_asterisk_matrix <- t(asterisk_matrix)

# remove repetitive comparisons, e.g., CAFs_vs_SMCs and SMCs_vs_CAFs - only keep one of them in ----
# Get the rownames of t_cor_matrix
rownames <- rownames(t_cor_matrix)

# manually specify comparisons of interest
keep_rows <- c("pnCAFs_score_vs_NPF_score", "PNS_glial_score_vs_NPF_score", "pnCAFs_score_vs_PNS_glial_score")
filtered_rownames <- rownames[rownames %in% keep_rows]

# Now you can subset t_cor_matrix using the filtered rownames
filtered_t_cor_matrix <- t_cor_matrix[filtered_rownames, ]

# now filter t_asterisk_matrix to match filtered_t_cor_matrix
# Filter t_asterisk_matrix based on filtered_rownames
filtered_t_asterisk_matrix <- t_asterisk_matrix[filtered_rownames, ]

# get info from spread_df (should already be filtered, but just in case any samples have squeezed through, somehow)
anno_df <- spread_df %>% filter(sample_id %in% all_of(colnames(filtered_t_cor_matrix))) %>% distinct(sample_id, type) %>% remove_rownames() %>% column_to_rownames("sample_id")
# actuall add sample_id as a column again
anno_df$sample_id <- rownames(anno_df)

anno_cols <- list(
    type = c("cancer" = "#9065cf",
             "adj_benign" = "#E0B57C"),
    sample_id = c("20216-1" = "#E69F00",
                  "19617-2" = "#56B4E9",
                  "20111-2" = "#009E73",
                  "20033" = "#F0E442",
                  "20130-2" = "#0072B2",
                  "20153-2" = "#D55E00",
                  "20153-1" = "#CC79A7",   
                  "20272-2" = "#999999",   
                  "20130-1" = "#F4A82F",   
                  "20128-1" = "#A6D854"    
    )
)

# order the matrix by adjacent benign samples, followed by cancer
# desired_order = c('20111-2', '20153-2', '20130-2', '20033', '20216-1', '19617-2', '20153-1', '20272-2', '20130-1', '20128-1')

# reindex the df with the desired order of columns
# filtered_t_cor_matrix = filtered_t_cor_matrix[, desired_order]
# filtered_t_asterisk_matrix = filtered_t_asterisk_matrix[, desired_order]

# Create the heatmap using pheatmap
p <- pheatmap(
  filtered_t_cor_matrix,
  annotation = anno_df,
  annotation_colours = anno_cols,
  display_numbers = filtered_t_asterisk_matrix,
  fontsize_number = 8,  # Adjust font size for better visibility
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  col = colorRampPalette(c("dodgerblue4", "white", "red4"))(101),
 #  breaks = seq(min(filtered_combined_cor_df$Correlation), max(filtered_combined_cor_df$Correlation), length.out = 51),  # Specify the breaks for the color scale
  # breaks = seq(min(combined_cor_df$Correlation), max(combined_cor_df$Correlation), length.out = 51),  # Specify the breaks for the color scale
  breaks = seq(-0.5, 0.5, length.out = 101),
  na_col = "gray",
  border_color = NA,
  cellwidth = 14, cellheight = 14, 
  main = "Pearson correlation heatmap\nScanpy scores\n(based on cepo-derived genes)",
  ylab = "Cell type Pair",
  xlab = "Sample IDs",
  scale = "none"
)

# Print the heatmap plot
pdf(paste0(figDir, "scanpy_scores_without_DCN_nerve_enriched.pdf"), width = 18, height=18)
print(p)
dev.off()

# Create the heatmap using pheatmap
p <- pheatmap(
  filtered_t_cor_matrix,
  annotation = anno_df,
  annotation_colours = anno_cols,
  display_numbers = filtered_t_asterisk_matrix,
  fontsize_number = 8,  # Adjust font size for better visibility
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  col = colorRampPalette(c("dodgerblue4", "white", "red4"))(101),
  # breaks = seq(min(combined_cor_df$Correlation), max(combined_cor_df$Correlation), length.out = 51),  # Specify the breaks for the color scale
  breaks = seq(-0.5, 0.5, length.out = 101),
  na_col = "gray",
  border_color = NA,
  cellwidth = 14, cellheight = 14, 
  main = "Pearson correlation heatmap\nScanpy scores\n(based on cepo-derived genes)",
  ylab = "Cell type Pair",
  xlab = "Sample IDs",
  scale = "none"
)

# Print the heatmap plot
pdf(paste0(figDir, "scanpy_scores_clustered_without_DCN_nerve_enriched.pdf"), width = 18, height=18)
print(p)
dev.off()

# ---------------------------------------------------------------------------
# REPEAT ALL OF THE ABOVE, NOW USING HISTOLOGY AS A GROUPING VARIABLE -------
# RATHER THAN SAMPLE ID -----------------------------------------------------
# ---------------------------------------------------------------------------

# ---------------------------------------------------
# plot scatter plots - per histological feature -----
# ---------------------------------------------------

# reload/recreate df

spread_df <- spread(merged_df, key = score_name, value = score_value)

# exclude "Exclude" histology from the analysis - this will be spots outside of the tissue, etc
spread_df <- spread_df %>% filter(Histology != "Exclude")
spread_df <- spread_df %>% filter(sample_id != "Exclude")

# in Histology column, swap empty spaces with an underscore
spread_df$Histology <- gsub(" ", "_", spread_df$Histology)

# replace Nerve/Ganglia with Nerve only
spread_df$Histology <- gsub("Nerve/Ganglia", "Nerve", spread_df$Histology)

# match cases/histology in GG4_Cribriform/Stroma prostate/Vessels
spread_df$Histology <- gsub("GG4_cribriform", "GG4_Cribriform", spread_df$Histology)
spread_df$Histology <- gsub("Stroma_prostate", "Stroma_prostatic", spread_df$Histology)
spread_df$Histology <- gsub("Vessels", "Vessel", spread_df$Histology)
spread_df$Histology <- gsub("Epi_Benign_\\(transitional\\)", "Epi_Benign_transitional", spread_df$Histology)

# also remove unannotated spots (just one in this case)
spread_df <- spread_df[spread_df$Histology != "", ]

# get annotation info from spread_df
# create the has_vessels and has_nerves columns based on the presence of "Vessel" and "Nerve" in the Histology column
spread_df <- spread_df %>%
  group_by(sample_id) %>%
  mutate(
    has_vessels = if_else(any(grepl("Vessel", Histology)), "yes", "no"),
    nerve_enriched = if_else(any(grepl("Nerve", Histology)), "yes", "no")
  ) %>%
  ungroup()


# visualise all scores with 0.5 as threshold
p <- ggplot(spread_df, aes(x = pnCAFs_score, y = PNS_glial_score, color = nerve_enriched)) +
  geom_point() +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  theme_minimal()

pdf(paste0(figDir, "all_scores_0.5_threshold.pdf"))
print(p)
dev.off()

# with 0 as threshold
p <- ggplot(spread_df, aes(x = pnCAFs_score, y = PNS_glial_score, color = nerve_enriched)) +
  geom_point() +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  theme_minimal()

pdf(paste0(figDir, "all_scores_0_threshold.pdf"))
print(p)
dev.off()

# filter for nerve-enriched samples and filter for signature threshold above >0.5
spread_df2 <- spread_df %>% 
              filter(PNS_glial_score > 0.5 & pnCAFs_score > 0.5)
             # filter(nerve_enriched == "yes") # %>%


spread_df2 <- spread_df2 %>% filter(Histology %in% c('Epi_Benign', 'GG3', 'GG4', 'Inflammation', 'Stroma_prostatic', 'Stroma_extraprostatic', 'Vessel', 'Nerve', 'Adipose'))             

cell_types_of_interest <- c('pnCAFs_score', 'PNS_glial_score') #, 'NPF_score')

# desired_order = c('Epi_Benign', 'Epi_Benign_transitional', 'GG3', 'GG4', 'GG4_Cribriform', 'Inflammation', 'Stroma_prostatic', 'Stroma_extraprostatic', 'Vessel', 'Nerve', 'Adipose')
desired_order = c('Epi_Benign', 'GG3', 'GG4', 'Inflammation', 'Stroma_prostatic', 'Stroma_extraprostatic', 'Vessel', 'Nerve', 'Adipose')

spread_df2$Histology <- factor(spread_df2$Histology, levels = desired_order)
spread_df2 <- spread_df2[order(spread_df2$Histology), ]

# only include Histology with more than 10 spots
# spread_df2 <- spread_df2 %>%
#               group_by(Histology) %>%
#               filter(n() > 10) %>%
#               ungroup()

plot_list <- list()

for(feature in unique(spread_df2$Histology)){
  print(feature)

  df_subset <- spread_df2 %>% dplyr::filter(Histology == feature) %>% select('sample_id', 'Barcode', 'pnCAFs_score', 'PNS_glial_score')
  # this will barcodes from multiple samples so there will be duplicates, make them unique
  df_subset$Barcode <- paste(df_subset$Barcode, df_subset$sample_id)
  df_subset <- df_subset %>% remove_rownames() %>% column_to_rownames('Barcode')

  sp <- ggscatter(df_subset, x='pnCAFs_score', y='PNS_glial_score',
     add = "reg.line",  # Add regressin line
     add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
     conf.int = TRUE # Add confidence interval
     ) + 
     geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
     geom_vline(xintercept = 0.5, linetype = "dashed", color = "red")

  # Add correlation coefficient
  calculate_label_positions <- function(df_subset) {
    # Calculate label positions based on data range
    label_x <- min(df_subset$pnCAFs_score) + 0.2 * diff(range(df_subset$pnCAFs_score))
    label_y <- min(df_subset$PNS_glial_score) + 0.5 * diff(range(df_subset$PNS_glial_score))
    return(list(label_x = label_x, label_y = label_y))
  }

  label_positions <- calculate_label_positions(df_subset)

  plot_list[[feature]] <- sp + stat_cor(method = "pearson", label.x = label_positions$label_x, label.y = label_positions$label_y, color="red") + labs(subtitle = feature)

}

pdf(paste0(figDir, "feature_scatter_plot_without_DCN_nerve_enriched_filtered_at_0.5.pdf"), width=9, height=6)
p <- cowplot::plot_grid(plotlist = plot_list)
print(p)
dev.off()

# List of cell types of interest (lineage 1 comparison to all other cell types at minor level)
cell_types_of_interest <- c('pnCAFs_score', 'PNS_glial_score', 'NPF_score')

# Create an empty list to store the correlation data frames for each cell type
correlation_list <- list()

# Loop through each cell type of interest
for (cell_type in cell_types_of_interest) {
  print(cell_type)
  
  # Loop through each unique Histology and calculate the correlations for the current cell type
  for (histology in unique(as.factor(spread_df2$Histology))) {
    print(histology)
    
    # Calculate correlations for the current cell type and sample ID
    correlation_df <- calculate_correlations_histo(spread_df2, histology, cell_type)

    correlation_df$Histology <- histology
    correlation_df$cell_type <- cell_type
    
    # Store the correlation dataframe in the list
    correlation_list[[paste(cell_type, histology, sep = "_")]] <- correlation_df
  }
}

# Combine all correlation data frames into a single data frame
combined_cor_df <- do.call(rbind, correlation_list)

# add histo feature classification
spread_df <- spread_df %>%
  mutate(main_histology = case_when(
    grepl("^Epi|GG|PIN|Benign", Histology) ~ "Epithelial",
    grepl("^Transition_State", Histology) ~ "Epithelial",
    grepl("^Stroma", Histology) ~ "Stromal",
    grepl("flammation", Histology) ~ "Inflammation",
    grepl("^Vessel", Histology) ~ "Vessel",
    grepl("^Adipose|Fat", Histology) ~ "Adipose",
    grepl("^Nerve", Histology) ~ "Nerve"))

# Extract the interaction patterns from interactions_df
# interaction_patterns <- paste(interactions_df$celltype_1, "vs", interactions_df$celltype_2, sep = "_")

# Filter combined_cor_df to keep only the interactions found in interaction_patterns
# filtered_combined_cor_df <- combined_cor_df[combined_cor_df$Celltype_Pair %in% interaction_patterns, ]

# Prepare the correlation matrix
cor_matrix <- combined_cor_df %>% dplyr::select(-P_Value)

# Assuming your correlation matrix is named cor_matrix
cor_matrix <- cor_matrix %>%
  select(-cell_type) %>%
  pivot_wider(names_from = Celltype_Pair, values_from = Correlation)   %>%
  tibble::remove_rownames() %>% column_to_rownames(var = "Histology")

# Prepare the p-value matrix
p_value_matrix <- combined_cor_df %>% dplyr::select(Celltype_Pair, P_Value, Histology) # %>% tibble::remove_rownames() %>% column_to_rownames(var = "sample_id")

# Assuming your correlation matrix is named cor_matrix
p_value_matrix_wide <- p_value_matrix %>%
  pivot_wider(names_from = Celltype_Pair, values_from = P_Value) %>%
  tibble::remove_rownames() %>% column_to_rownames(var = "Histology")

# # Count the total number of tests
# m <- 30*6
# 
# # Apply Bonferroni correction to each element in the matrix
# corrected_matrix <- p_value_matrix_wide * m
# 
# # Ensure that the corrected p-values are capped at 1
# corrected_matrix[corrected_matrix > 1] <- 1

# Apply the function to the p-value matrix
asterisk_matrix <- apply(p_value_matrix_wide, 2, function(column) sapply(column, replace_with_asterisks))

t_cor_matrix <- t(cor_matrix)
t_asterisk_matrix <- t(asterisk_matrix)

# remove repetitive comparisons, e.g., CAFs_vs_SMCs and SMCs_vs_CAFs - only keep one of them in ----
# Get the rownames of t_cor_matrix
rownames <- rownames(t_cor_matrix)

# manually specify comparisons of interest
keep_rows <- c("pnCAFs_score_vs_NPF_score", "PNS_glial_score_vs_NPF_score", "pnCAFs_score_vs_PNS_glial_score")
filtered_rownames <- rownames[rownames %in% keep_rows]

# Now you can subset t_cor_matrix using the filtered rownames
filtered_t_cor_matrix <- t_cor_matrix[filtered_rownames, ]

# now filter t_asterisk_matrix to match filtered_t_cor_matrix
# Filter t_asterisk_matrix based on filtered_rownames
filtered_t_asterisk_matrix <- t_asterisk_matrix[filtered_rownames, ]

# get info from spread_df (should already be filtered, but just in case any samples have squeezed through, somehow)
anno_df <- spread_df %>% filter(Histology %in% all_of(colnames(filtered_t_cor_matrix))) %>% distinct(Histology, main_histology) %>% remove_rownames() %>% column_to_rownames("Histology")
# anno_cols <- list(type = type_cols_darker)

# order the matrix by histology features
# desired_order = c("Benign", "Benign*", 'Epi_Benign', 'Epi_Benign_transitional', "Transition_State", "PIN", 'GG1', 'GG2', 'GG3', 'GG4', 'GG4_Cribriform', 'Inflammation', 'Chronic_inflammation', 'Stroma', 'Stroma_prostatic', 'Stroma_extraprostatic', 'Vessel', 'Nerve', 'Adipose', 'Fat')
# desired_order = c('Epi_Benign', 'Epi_Benign_transitional', 'GG3', 'GG4', 'GG4_Cribriform', 'Inflammation', 'Stroma_prostatic', 'Stroma_extraprostatic', 'Vessel', 'Nerve', 'Adipose')
desired_order = c('Epi_Benign', 'GG3', 'GG4', 'GG4_Cribriform', 'Inflammation', 'Stroma_prostatic', 'Stroma_extraprostatic', 'Vessel', 'Nerve', 'Adipose')

# reindex the df with the desired order of columns
filtered_t_cor_matrix = filtered_t_cor_matrix[, desired_order]
filtered_t_asterisk_matrix = filtered_t_asterisk_matrix[, desired_order]

# Create the heatmap using pheatmap
p <- pheatmap(
  filtered_t_cor_matrix,
  annotation = anno_df,
  # annotation_colours = anno_cols,
  display_numbers = filtered_t_asterisk_matrix,
  fontsize_number = 8,  # Adjust font size for better visibility
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  col = colorRampPalette(c("dodgerblue4", "white", "red4"))(101),
  # breaks = seq(min(combined_cor_df$Correlation), max(combined_cor_df$Correlation), length.out = 51),  # Specify the breaks for the color scale
  breaks = seq(-0.5, 0.5, length.out = 101),
  na_col = "gray",
  border_color = NA,
  cellwidth = 14, cellheight = 14, 
  main = "Pearson correlation heatmap\nScanpy scores\n(based on cepo-derived genes)",
  ylab = "Cell type Pair",
  xlab = "Sample IDs",
  scale = "none"
)

# Print the heatmap plot
pdf(paste0(figDir, "scanpy_scores_histology_without_DCN_nerve_enriched.pdf"), width = 15, height=15)
print(p)
dev.off()

# Create the heatmap using pheatmap
p <- pheatmap(
  filtered_t_cor_matrix,
  annotation = anno_df,
  # annotation_colours = anno_cols,
  display_numbers = filtered_t_asterisk_matrix,
  fontsize_number = 8,  # Adjust font size for better visibility
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  col = colorRampPalette(c("dodgerblue4", "white", "red4"))(101),
  # breaks = seq(min(combined_cor_df$Correlation), max(combined_cor_df$Correlation), length.out = 51),  # Specify the breaks for the color scale
  breaks = seq(-0.5, 0.5, length.out = 101),
  na_col = "gray",
  border_color = NA,
  cellwidth = 14, cellheight = 14, 
  main = "Pearson correlation heatmap\nScanpy scores\n(based on cepo-derived genes)",
  ylab = "Cell type Pair",
  xlab = "Sample IDs",
  scale = "none"
)

# Print the heatmap plot
pdf(paste0(figDir, "scanpy_scores_histology_clustered_without_DCN_nerve_enriched.pdf"), width = 15, height=15)
print(p)
dev.off()

# SAVE LOGS -----------------------------------------------------------------

# save sessionInfo
writeLines(capture.output(sessionInfo()), paste0(resultDir, "sessionInfo.txt"))