#processing all pairwise combinations of traits and generate the matrix one trait at the time

# Load libraries
library(dplyr)
library(tibble)
library(data.table)

# Source your processing functions
source("processing_data.R")

# Define input/output
input_dir <- "preprocessed_data/"
output_dir <- "processed_data/"
distance_threshold <- 200000
significance_threshold <- 1e-5

# Read all traits
files <- list.files(input_dir, pattern = "_snp_significant.tsv$", full.names = FALSE)
traits <- gsub("_snp_significant.tsv", "", files)
traits <- sort(unique(traits))
traits <- traits[!grepl("_effect$", traits)]  # Don't compute '_effect'

# Initialize matrix if it does not exist
matrix_file <- paste0(output_dir, "significant_snp_matrix.tsv")

if (!file.exists(matrix_file)) {
  trait_matrix <- matrix(0, nrow = length(traits), ncol = length(traits))
  rownames(trait_matrix) <- traits
  colnames(trait_matrix) <- traits
  fwrite(cbind(Trait = rownames(trait_matrix), as.data.frame(trait_matrix)),
         file = matrix_file, sep = "\t", row.names = FALSE)
}

# Reload matrix each time
trait_matrix <- fread(matrix_file) %>% as.data.frame()
rownames(trait_matrix) <- trait_matrix$Trait
trait_matrix$Trait <- NULL

# ---- Function to process one trait at a time ---- #
process_single_trait <- function(base_trait) {
  
  cat("\n🔵 Processing trait:", base_trait, "\n")
  
  for (target_trait in traits) {
    if (base_trait == target_trait) next
    
    cat("🧩 ", base_trait, "->", target_trait, "\n")
    
    # File paths
    sig_trait1_file <- paste0(input_dir, base_trait, "_snp_significant.tsv")
    clean_trait2_file <- paste0(input_dir, target_trait, "_snp_clean.tsv")
    effect_clean_trait2_file <- paste0(input_dir, target_trait, "_effect_snp_clean.tsv")
    file_name <- paste0(base_trait, "_", target_trait)
    
    # 1. Merge
    merged_data <- tryCatch({
      tidy_snps(sig_trait1_file, effect_clean_trait2_file, clean_trait2_file)
    }, error = function(e) {
      message("⚠️ Skipping due to error: ", e$message)
      return(NULL)
    })
    
    if (!is.null(merged_data) && nrow(merged_data) > 0) {
      
      cat("🔹 Merged SNPs:", nrow(merged_data), "\n")
      
      # 2. Significance labeling
      merged_data <- test_significance(merged_data, significance_threshold)
      
      # 3. Distance filtering
      filtered_data <- filter_snps_by_distance(merged_data, distance_threshold)
      cat("🔹 SNPs after distance filtering:", nrow(filtered_data), "\n")
      
      # 4. Count significant filtered SNPs
      num_sig_snps <- sum(filtered_data$p_value < significance_threshold, na.rm = TRUE)
      cat("🔹 Significant & filtered SNPs:", num_sig_snps, "\n")
      
      # Save filtered file
      output_filtered_file <- paste0(output_dir, file_name, "_filtered_snps.tsv")
      fwrite(filtered_data, file = output_filtered_file, sep = "\t")
      
      # Update matrix
      trait_matrix[base_trait, target_trait] <<- num_sig_snps
      
    } else {
      cat("⚠️ No merged SNPs found for this pair.\n")
      trait_matrix[base_trait, target_trait] <<- 0
    }
    
  }
  
  # Save matrix after finishing the base_trait row
  fwrite(cbind(Trait = rownames(trait_matrix), as.data.frame(trait_matrix)),
         file = matrix_file, sep = "\t", row.names = FALSE)
  
  cat("\n✅ Finished processing:", base_trait, "\n")
}

# ---- How to use ---- #
process_single_trait("DIET_9364")

