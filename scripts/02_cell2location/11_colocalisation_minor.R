# spatial co-localisation analysis, based on cell2location output
# across entire section, cell-to-cell

# conda env: atac_2022

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
exp = "02_cell2location"
analysis = "11_colocalisation_minor"

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


# > table(histo_df$sample_id)
# 
# 19617-2   20033 20111-2 20128-1 20130-1 20130-2 20153-1 20153-2 20216-1 20272-2
#    1546    2941    2213    2186    2152    1477    1284    1248    1718    1573
# Exclude
#    1805


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

# > barcode_counts
# # A tibble: 10 × 2
#    section_name         unique_barcode_count
#    <chr>                               <int>
#  1 A1_20033                             2941
#  2 A1_20216_C1                          1718
#  3 B1_19617_C2                          1546
#  4 B1_20111_C2                          2213
#  5 C1_20130_C1                          2180
#  6 C1_20153_C1                          1517
#  7 D1_20130_C2                          1477
#  8 D1_20153_C2                          1248
#  9 PCa20130_C1_20272_C2                 4753
# 10 PCa20153_C1_20128_C1                 4247

# combine the two dfs --------------------
# keeping only the 6 samples that have histo anno

# Merge based on both sample_id and Barcode
merged_df <- inner_join(df_tidy, histo_df, by = c("section_name", "Barcode"))

spread_df <- spread(merged_df, key = celltype_minor_v2, value = c2l_value)

# > dim(spread_df)
# [1] 20143    21

# exclude "Exclude" histology and blank cells from the analysis - this will be spots outside of the tissue, etc
spread_df <- spread_df %>% filter(Histology != "Exclude" & Histology != "")

# > dim(spread_df)
# [1] 16977    21

# exclude "Exclude" in sample_id from the analysis
spread_df <- spread_df %>% filter(sample_id != "Exclude")

# > dim(spread_df)
# [1] 16976    21

# ---------------------------------------------------
# test correlations, with significance --------------
# ---------------------------------------------------

## add pearson correlation coefficients to heatmap
## need to work out if i need multiple testing correction

# Function to calculate Pearson correlation and p-values for a given sample_id
calculate_correlations <- function(df, sample_id, celltype) {
  # Filter the data for the given sample_id
  sample_data <- df %>% ungroup() %>%
    filter(sample_id == !!sample_id) %>% dplyr::select(-sample_id, -Barcode, -type, -Histology, -section_name)

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

# List of cell types of interest (lineage 1 comparison to all other cell types at minor level)
# interactions_df <- read.csv("/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/config/minor_interactions_of_interest_v3.csv")
# interactions_df <- read.csv("/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/config/minor_interactions_CAFs_all.csv")
interactions_df <- read.csv("/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/config/minor_interactions_CAFs_subset.csv")

cell_types_of_interest <- unique(interactions_df$celltype_1)

# Create an empty list to store the correlation data frames for each cell type
correlation_list <- list()

# Loop through each cell type of interest
for (cell_type in cell_types_of_interest) {
  print(cell_type)
  
  # Loop through each unique sample_id and calculate the correlations for the current cell type
  for (sample_id in unique(as.factor(spread_df$sample_id))) {
    print(sample_id)
    
    # Calculate correlations for the current cell type and sample ID
    correlation_df <- calculate_correlations(spread_df, sample_id, cell_type)
    
    # Store the correlation dataframe in the list
    correlation_list[[paste(cell_type, sample_id, sep = "_")]] <- correlation_df
  }
}

# Combine all correlation data frames into a single data frame
combined_cor_df <- do.call(rbind, correlation_list)

# Add cell type and sample_id as columns
combined_cor_df$cell_type <- gsub("_\\d+.*", "", rownames(combined_cor_df))
combined_cor_df$sample_id <- gsub("[A-Za-z_]+_", "", rownames(combined_cor_df))
# remove the decimal part in the sample_id column
combined_cor_df$sample_id <- sub("\\..*$", "", combined_cor_df$sample_id)


# Extract the interaction patterns from interactions_df
interaction_patterns <- paste(interactions_df$celltype_1, "vs", interactions_df$celltype_2, sep = "_")

# Filter combined_cor_df to keep only the interactions found in interaction_patterns
filtered_combined_cor_df <- combined_cor_df[combined_cor_df$Celltype_Pair %in% interaction_patterns, ]

# Prepare the correlation matrix
cor_matrix <- filtered_combined_cor_df %>% dplyr::select(-P_Value)

# Assuming your correlation matrix is named cor_matrix
cor_matrix <- cor_matrix %>%
  select(-cell_type) %>%
  pivot_wider(names_from = Celltype_Pair, values_from = Correlation)   %>%
  tibble::remove_rownames() %>% column_to_rownames(var = "sample_id")

# Prepare the p-value matrix
p_value_matrix <- filtered_combined_cor_df %>% dplyr::select(Celltype_Pair, P_Value, sample_id) # %>% tibble::remove_rownames() %>% column_to_rownames(var = "sample_id")

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

# Apply the function to the p-value matrix
asterisk_matrix <- apply(p_value_matrix_wide, 2, function(column) sapply(column, replace_with_asterisks))

t_cor_matrix <- t(cor_matrix)
t_asterisk_matrix <- t(asterisk_matrix)

# remove repetitive comparisons, e.g., CAFs_vs_SMCs and SMCs_vs_CAFs - only keep one of them in ----
# Get the rownames of t_cor_matrix
rownames <- rownames(t_cor_matrix)

# Function to check if a string has a reversed counterpart in a vector
has_reversed_counterpart <- function(string, vector) {
 split_string <- unlist(strsplit(string, "_vs_"))
 reversed_string <- paste(rev(split_string), collapse = "_vs_")
 any(reversed_string == vector)
}

# Identify row names with reversed counterparts
reversed_patterns <- sapply(rownames, has_reversed_counterpart, vector = rownames)

# Identify paired reverse strings
paired_reverse_strings <- rownames[reversed_patterns]

# Create a logical vector to keep only one element from each pair
keep_first_in_pair <- logical(length(paired_reverse_strings))

for (i in seq_along(paired_reverse_strings)) {

  pair <- paired_reverse_strings[i]
  index_pair <- i

  parts <- unlist(strsplit(pair, "_vs_"))
  # Reverse the parts and reassemble with "_vs_"
  reverse_pair <- paste(rev(parts), collapse = "_vs_")
  #reverse_pair <- paste(rev(strsplit(pair, "_")[[1]]), collapse = "_")
  index_reverse_pair <- which(paired_reverse_strings == reverse_pair)

  if (index_pair < index_reverse_pair) {
    keep_first_in_pair[i] <- TRUE
  } else if (index_pair == index_reverse_pair) {
    keep_first_in_pair[i] <- TRUE
  }
}

kept_pairs <- paired_reverse_strings[keep_first_in_pair]
kept_patterns <- reversed_patterns[!reversed_patterns]
kept_patterns <- names(kept_patterns)

# rows to keep 
keep_rows <- c(kept_pairs, kept_patterns)
filtered_rownames <- rownames[rownames %in% keep_rows]

# Now you can subset t_cor_matrix using the filtered rownames
filtered_t_cor_matrix <- t_cor_matrix[filtered_rownames, ]

# now filter t_asterisk_matrix to match filtered_t_cor_matrix
# Filter t_asterisk_matrix based on filtered_rownames
filtered_t_asterisk_matrix <- t_asterisk_matrix[filtered_rownames, ]

# get info from spread_df (should already be filtered, but just in case any samples have squeezed through, somehow)
spread_df2 <- spread_df %>% ungroup() %>% select(-Barcode)

anno_df <- spread_df2 %>% filter(sample_id %in% all_of(colnames(filtered_t_cor_matrix))) %>% distinct(sample_id, type) %>% remove_rownames() %>% column_to_rownames("sample_id")
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
desired_order = c('20111-2', '20153-2', '20130-2', '20033', '20216-1', '19617-2', '20153-1', '20272-2', '20130-1', '20128-1')

# reindex the df with the desired order of columns
filtered_t_cor_matrix = filtered_t_cor_matrix[, desired_order]
filtered_t_asterisk_matrix = filtered_t_asterisk_matrix[, desired_order]


# WHEN PLOTTING CAF popn's SEPARATELY --------------------------------

plot_heatmaps <- function(cor_matrix, asterisk_matrix) {
  # Function to split the matrix based on row names
  split_matrix <- function(matrix) {
    npf_like_matrix <- matrix[grep("^NPF_like", rownames(matrix)), ]
    pnCAFs_matrix <- matrix[grep("^pnCAFs", rownames(matrix)), ]
    uniCAFs_matrix <- matrix[grep("^uniCAFs", rownames(matrix)), ]
    whCAFs_matrix <- matrix[grep("^whCAFs", rownames(matrix)), ]
    matrices <- list(NPF_like = npf_like_matrix, pnCAFs = pnCAFs_matrix, uniCAFs = uniCAFs_matrix, whCAFs = whCAFs_matrix)
    return(matrices)
  }
  
  # Split the correlation matrix and asterisk matrix
  cor_matrices <- split_matrix(cor_matrix)
  asterisk_matrices <- split_matrix(asterisk_matrix)
  
  # Create heatmaps for each pair of matrices
  for (key in names(cor_matrices)) {
    cor_matrix <- cor_matrices[[key]]
    asterisk_matrix <- asterisk_matrices[[key]]
    
    p <- pheatmap(
      cor_matrix,
      annotation = anno_df,
      annotation_colours = anno_cols,
      display_numbers = asterisk_matrix,
      fontsize_number = 8,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      col = colorRampPalette(c("dodgerblue4", "white", "red4"))(50),
      na_col = "gray",
      border_color = NA,
      cellwidth = 11, cellheight = 14, 
      main = paste("Heatmap for", key),
      ylab = "Cell type Pair",
      xlab = "Sample IDs",
      scale = "none"
    )
    
    # Print the heatmap plot
    pdf(paste0(figDir, key, "_minor_cell2location.pdf"), width = 15, height=15)
    print(p)
    dev.off()

    p <- pheatmap(
      cor_matrix,
      annotation = anno_df,
      annotation_colours = anno_cols,
      display_numbers = asterisk_matrix,
      fontsize_number = 8,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      col = colorRampPalette(c("dodgerblue4", "white", "red4"))(50),
      na_col = "gray",
      border_color = NA,
      cellwidth = 11, cellheight = 14, 
      main = paste("Heatmap for", key),
      ylab = "Cell type Pair",
      xlab = "Sample IDs",
      scale = "none"
    )
    
    # Print the heatmap plot
    pdf(paste0(figDir, key, "_minor_cell2location_clustered.pdf"), width = 15, height=15)
    print(p)
    dev.off()
  }
}

# Usage: Call the function with your filtered_t_cor_matrix and filtered_t_asterisk_matrix
plot_heatmaps(filtered_t_cor_matrix, filtered_t_asterisk_matrix)

# WHEN PLOTTING ALL TOGETHER -----------------------------------------

# Create the heatmap using pheatmap
p <- pheatmap(
  filtered_t_cor_matrix,
  annotation = anno_df,
  annotation_colours = anno_cols,
  display_numbers = filtered_t_asterisk_matrix,
  fontsize_number = 8,  # Adjust font size for better visibility
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  col = colorRampPalette(c("dodgerblue4", "white", "red4"))(50),
  breaks = seq(min(filtered_combined_cor_df$Correlation), max(filtered_combined_cor_df$Correlation), length.out = 51),  # Specify the breaks for the color scale
  na_col = "gray",
  border_color = NA,
  cellwidth = 11, cellheight = 14, 
  main = "Pearson correlation heatmap",
  ylab = "Cell type Pair",
  xlab = "Sample IDs",
  scale = "none"
)

# Print the heatmap plot
pdf(paste0(figDir, "celltype_minor_cell2location_v3_new.pdf"), width = 15, height=15)
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
  col = colorRampPalette(c("dodgerblue4", "white", "red4"))(50),
  breaks = seq(min(filtered_combined_cor_df$Correlation), max(filtered_combined_cor_df$Correlation), length.out = 51),  # Specify the breaks for the color scale
  na_col = "gray",
  border_color = NA,
  cellwidth = 11, cellheight = 14, 
  main = "Pearson correlation heatmap",
  ylab = "Cell type Pair",
  xlab = "Sample IDs",
  scale = "none"
)

# Print the heatmap plot
pdf(paste0(figDir, "celltype_minor_cell2location_clustered_v3_new.pdf"), width = 15, height=15)
print(p)
dev.off()

# ---------------------------------------------------------------------------
# REPEAT ALL OF THE ABOVE, NOW USING HISTOLOGY AS A GROUPING VARIABLE -------
# RATHER THAN SAMPLE ID -----------------------------------------------------
# ---------------------------------------------------------------------------

# reload/create df

spread_df <- spread(merged_df, key = celltype_minor_v2, value = c2l_value)

# exclude "Exclude" histology from the analysis - this will be spots outside of the tissue, etc
spread_df <- spread_df %>% filter(Histology != "Exclude")

# in Histology column, swap empty spaces with an underscore
spread_df$Histology <- gsub(" ", "_", spread_df$Histology)

# replace Nerve/Ganglia with Nerve only
spread_df$Histology <- gsub("Nerve/Ganglia", "Nerve", spread_df$Histology)

# match cases/histology in GG4_Cribriform/Stroma prostate/Vessels
spread_df$Histology <- gsub("GG4_cribriform", "GG4_Cribriform", spread_df$Histology)
spread_df$Histology <- gsub("Stroma_prostate", "Stroma_prostatic", spread_df$Histology)
spread_df$Histology <- gsub("Vessels", "Vessel", spread_df$Histology)

# also remove unannotated spots (just one in this case)
spread_df <- spread_df[spread_df$Histology != "", ]

# fix transitional
spread_df$Histology <- gsub("\\(|\\)", "", spread_df$Histology)

# exclude "Exclude" in sample_id from the analysis
spread_df <- spread_df %>% filter(sample_id != "Exclude")
# > dim(spread_df)
# [1] 16976    21

# fix transitional
spread_df$Histology <- gsub("\\(|\\)", "", spread_df$Histology)

# Function to calculate Pearson correlation and p-values for a given sample_id
calculate_correlations_histo <- function(df, histology, celltype) {
  # Filter the data for the given sample_id
  sample_data <- df %>% ungroup() %>%
    filter(Histology == !!histology) %>% select(-sample_id, -Barcode, -type, -Histology, -section_name)

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

# List of cell types of interest (lineage 1 comparison to all other cell types at minor level)
# interactions_df <- read.csv("/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/config/minor_interactions_of_interest_v3.csv")
cell_types_of_interest <- unique(interactions_df$celltype_1)

# Create an empty list to store the correlation data frames for each cell type
correlation_list <- list()

# Loop through each cell type of interest
for (cell_type in cell_types_of_interest) {
  print(cell_type)
  
  # Loop through each unique Histology and calculate the correlations for the current cell type
  for (histology in unique(as.factor(spread_df$Histology))) {
    print(histology)
    
    # Calculate correlations for the current cell type and sample ID
    correlation_df <- calculate_correlations_histo(spread_df, histology, cell_type)

    correlation_df$Histology <- histology
    
    # Store the correlation dataframe in the list
    correlation_list[[paste(cell_type, histology, sep = "_")]] <- correlation_df
  }
}

# Combine all correlation data frames into a single data frame
combined_cor_df <- do.call(rbind, correlation_list)

# Add cell type and Histology as columns
# Extract the cell type from Celltype_Pair
combined_cor_df$cell_type <- gsub("_vs.*", "", combined_cor_df$Celltype_Pair)

# add histo feature classification
spread_df <- spread_df %>%
  mutate(main_histology = case_when(
    grepl("^Epi|GG", Histology) ~ "Epithelial",
    grepl("^Stroma", Histology) ~ "Stromal",
    grepl("^Inflammation", Histology) ~ "Inflammation",
    grepl("^Vessel", Histology) ~ "Vessel",
    grepl("^Adipose", Histology) ~ "Adipose",
    grepl("^Nerve", Histology) ~ "Nerve"))

# Extract the interaction patterns from interactions_df
interaction_patterns <- paste(interactions_df$celltype_1, "vs", interactions_df$celltype_2, sep = "_")

# Filter combined_cor_df to keep only the interactions found in interaction_patterns
filtered_combined_cor_df <- combined_cor_df[combined_cor_df$Celltype_Pair %in% interaction_patterns, ]

# Prepare the correlation matrix
cor_matrix <- filtered_combined_cor_df %>% dplyr::select(-P_Value)

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

# Apply the function to the p-value matrix
asterisk_matrix <- apply(p_value_matrix_wide, 2, function(column) sapply(column, replace_with_asterisks))

t_cor_matrix <- t(cor_matrix)
t_asterisk_matrix <- t(asterisk_matrix)

# remove repetitive comparisons, e.g., CAFs_vs_SMCs and SMCs_vs_CAFs - only keep one of them in ----
# Get the rownames of t_cor_matrix
rownames <- rownames(t_cor_matrix)

# Function to check if a string has a reversed counterpart in a vector
has_reversed_counterpart <- function(string, vector) {
 split_string <- unlist(strsplit(string, "_vs_"))
 reversed_string <- paste(rev(split_string), collapse = "_vs_")
 any(reversed_string == vector)
}

# Identify row names with reversed counterparts
reversed_patterns <- sapply(rownames, has_reversed_counterpart, vector = rownames)

# Identify paired reverse strings
paired_reverse_strings <- rownames[reversed_patterns]

# Create a logical vector to keep only one element from each pair
keep_first_in_pair <- logical(length(paired_reverse_strings))

for (i in seq_along(paired_reverse_strings)) {

  pair <- paired_reverse_strings[i]
  index_pair <- i

  parts <- unlist(strsplit(pair, "_vs_"))
  # Reverse the parts and reassemble with "_vs_"
  reverse_pair <- paste(rev(parts), collapse = "_vs_")
  #reverse_pair <- paste(rev(strsplit(pair, "_")[[1]]), collapse = "_")
  index_reverse_pair <- which(paired_reverse_strings == reverse_pair)

  if (index_pair < index_reverse_pair) {
    keep_first_in_pair[i] <- TRUE
  } else if (index_pair == index_reverse_pair) {
    keep_first_in_pair[i] <- TRUE
  }
}

kept_pairs <- paired_reverse_strings[keep_first_in_pair]
kept_patterns <- reversed_patterns[!reversed_patterns]
kept_patterns <- names(kept_patterns)

# rows to keep 
keep_rows <- c(kept_pairs, kept_patterns)
filtered_rownames <- rownames[rownames %in% keep_rows]

# Now you can subset t_cor_matrix using the filtered rownames
filtered_t_cor_matrix <- t_cor_matrix[filtered_rownames, ]

# now filter t_asterisk_matrix to match filtered_t_cor_matrix
# Filter t_asterisk_matrix based on filtered_rownames
filtered_t_asterisk_matrix <- t_asterisk_matrix[filtered_rownames, ]

# get info from spread_df (should already be filtered, but just in case any samples have squeezed through, somehow)
spread_df2 <- spread_df %>% ungroup() %>% select(-Barcode)

anno_df <- spread_df2 %>% filter(Histology %in% all_of(colnames(filtered_t_cor_matrix))) %>% distinct(Histology, main_histology) %>% remove_rownames() %>% column_to_rownames("Histology")

# order the matrix by histology features
desired_order = c('Epi_Benign', 'Epi_Benign_transitional', 'GG3', 'GG4', 'GG4_Cribriform', 'Inflammation', 'Stroma_prostatic', 'Stroma_extraprostatic', 'Vessel', 'Nerve', 'Adipose')

# reindex the df with the desired order of columns
filtered_t_cor_matrix = filtered_t_cor_matrix[, desired_order]
filtered_t_asterisk_matrix = filtered_t_asterisk_matrix[, desired_order]


# WHEN PLOTTING CAF popn's SEPARATELY --------------------------------

plot_heatmaps <- function(cor_matrix, asterisk_matrix) {
  # Function to split the matrix based on row names
  split_matrix <- function(matrix) {
    npf_like_matrix <- matrix[grep("^NPF_like", rownames(matrix)), ]
    pnCAFs_matrix <- matrix[grep("^pnCAFs", rownames(matrix)), ]
    uniCAFs_matrix <- matrix[grep("^uniCAFs", rownames(matrix)), ]
    whCAFs_matrix <- matrix[grep("^whCAFs", rownames(matrix)), ]
    matrices <- list(NPF_like = npf_like_matrix, pnCAFs = pnCAFs_matrix, uniCAFs = uniCAFs_matrix, whCAFs = whCAFs_matrix)
    return(matrices)
  }
  
  # Split the correlation matrix and asterisk matrix
  cor_matrices <- split_matrix(cor_matrix)
  asterisk_matrices <- split_matrix(asterisk_matrix)
  
  # Create heatmaps for each pair of matrices
  for (key in names(cor_matrices)) {
    cor_matrix <- cor_matrices[[key]]
    asterisk_matrix <- asterisk_matrices[[key]]
    
    p <- pheatmap(
      cor_matrix,
      annotation = anno_df,
      display_numbers = asterisk_matrix,
      fontsize_number = 8,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      col = colorRampPalette(c("dodgerblue4", "white", "red4"))(50),
      na_col = "gray",
      border_color = NA,
      cellwidth = 11, cellheight = 14, 
      main = paste("Heatmap for", key),
      ylab = "Cell type Pair",
      xlab = "Sample IDs",
      scale = "none"
    )
    
    # Print the heatmap plot
    pdf(paste0(figDir, key, "_minor_cell2location_histology.pdf"), width = 15, height=15)
    print(p)
    dev.off()

    p <- pheatmap(
      cor_matrix,
      annotation = anno_df,
      display_numbers = asterisk_matrix,
      fontsize_number = 8,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      col = colorRampPalette(c("dodgerblue4", "white", "red4"))(50),
      na_col = "gray",
      border_color = NA,
      cellwidth = 11, cellheight = 14, 
      main = paste("Heatmap for", key),
      ylab = "Cell type Pair",
      xlab = "Sample IDs",
      scale = "none"
    )
    
    # Print the heatmap plot
    pdf(paste0(figDir, key, "_minor_cell2location_histology_clustered.pdf"), width = 15, height=15)
    print(p)
    dev.off()
  }
}

# Usage: Call the function with your filtered_t_cor_matrix and filtered_t_asterisk_matrix
plot_heatmaps(filtered_t_cor_matrix, filtered_t_asterisk_matrix)


# WHEN PLOTTING ALL TOGETHER -----------------------------------------

# Create the heatmap using pheatmap
p <- pheatmap(
  filtered_t_cor_matrix,
  annotation = anno_df,
  # annotation_colours = anno_cols,
  display_numbers = filtered_t_asterisk_matrix,
  fontsize_number = 8,  # Adjust font size for better visibility
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  col = colorRampPalette(c("dodgerblue4", "white", "red4"))(50),
  breaks = seq(min(combined_cor_df$Correlation), max(combined_cor_df$Correlation), length.out = 51),  # Specify the breaks for the color scale
  na_col = "gray",
  border_color = NA,
  cellwidth = 11, cellheight = 14, 
  main = "Pearson correlation heatmap\nCell2Location - cell type minor",
  ylab = "Cell type Pair",
  xlab = "Sample IDs",
  scale = "none"
)

# Print the heatmap plot
pdf(paste0(figDir, "celltype_minor_cell2location_histology_v3_neww.pdf"), width = 15, height=15)
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
  col = colorRampPalette(c("dodgerblue4", "white", "red4"))(50),
  breaks = seq(min(combined_cor_df$Correlation), max(combined_cor_df$Correlation), length.out = 51),  # Specify the breaks for the color scale
  na_col = "gray",
  border_color = NA,
  cellwidth = 11, cellheight = 14, 
  main = "Pearson correlation heatmap\nCell2Location - cell type minor",
  ylab = "Cell type Pair",
  xlab = "Sample IDs",
  scale = "none"
)

# Print the heatmap plot
pdf(paste0(figDir, "celltype_minor_cell2location_clustered_histology_v3_neww.pdf"), width = 15, height=15)
print(p)
dev.off()

# SAVE LOGS -----------------------------------------------------------------

# save sessionInfo
writeLines(capture.output(sessionInfo()), paste0(resultDir, "sessionInfo.txt"))

# # full object
# cells <- readRDS("/share/ScratchGeneral/evaapo/projects/PCa/results/20230714_collate_annotation/02_update_annotation/rObjects/PCa_atlas_annotated_Oct23_update.rds")