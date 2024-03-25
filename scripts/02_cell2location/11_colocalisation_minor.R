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

# link to histo annotations (copied manually from Dropbox) - only 6 files have it:
histo_path <- paste0(projectDir, "/", repo, "/data/20240318_reviewed_histo_annotation_FFPE_v1/")

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
  prop_df <- read.csv(paste0("/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/results/02_cell2location/03a_cell_type_mapping/", sample_id, "/objects/q05_cell_abundance_w_sf.csv"))
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

spread_df <- spread(merged_df, key = celltype_major_v2, value = c2l_value)

# exclude "Exclude" histology from the analysis - this will be spots outside of the tissue, etc
spread_df <- spread_df %>% filter(Histology != "Exclude")

# ---------------------------------------------------
# test correlations, with significance --------------
# ---------------------------------------------------

## add pearson correlation coefficients to heatmap
## need to work out if i need multiple testing correction

# Function to calculate Pearson correlation and p-values for a given sample_id
calculate_correlations <- function(df, sample_id, celltype) {
  # Filter the data for the given sample_id
  sample_data <- df %>%
    filter(sample_id == !!sample_id) %>% select(-sample_id, -Barcode, -type, -Histology)

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
interactions_df <- read.csv("/share/ScratchGeneral/evaapo/projects/PCa_Visium/pca_visium_analysis/data/cell2location/minor_interactions_of_interest.csv")
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
  reverse_pair <- paste(rev(strsplit(pair, "_")[[1]]), collapse = "_")
  if (pair < reverse_pair) {
    keep_first_in_pair[i] <- TRUE
  } else if (pair == reverse_pair) {
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
anno_df <- spread_df %>% filter(sample_id %in% all_of(colnames(filtered_t_cor_matrix))) %>% distinct(sample_id, type) %>% remove_rownames() %>% column_to_rownames("sample_id")
anno_cols <- list(type = type_cols_darker)

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
  cellwidth = 14, cellheight = 14, 
  main = "Pearson correlation heatmap",
  ylab = "Cell type Pair",
  xlab = "Sample IDs",
  scale = "none"
)

# Print the heatmap plot
pdf(paste0(figDir, "celltype_minor_cell2location.pdf"), width = 15, height=15)
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
  cellwidth = 14, cellheight = 14, 
  main = "Pearson correlation heatmap",
  ylab = "Cell type Pair",
  xlab = "Sample IDs",
  scale = "none"
)

# Print the heatmap plot
pdf(paste0(figDir, "celltype_minor_cell2location_clustered.pdf"), width = 15, height=15)
print(p)
dev.off()


# SAVE LOGS -----------------------------------------------------------------

# save sessionInfo
writeLines(capture.output(sessionInfo()), paste0(resultDir, "sessionInfo.txt"))