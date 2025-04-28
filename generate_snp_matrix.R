#processing all pairwise combinations of traits and generate the matrix

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
traits <- sort(unique(traits)) # (optional)

#don't compute '_effect'
traits <- traits[!grepl("_effect$", traits)]

# Initialize asymmetric matrix
trait_matrix <- matrix(0, nrow = length(traits), ncol = length(traits))
rownames(trait_matrix) <- traits
colnames(trait_matrix) <- traits

# Loop over trait pairs
for (i in seq_along(traits)) {
  for (j in seq_along(traits)) {
    if (i != j) {
      
      trait1 <- traits[i]
      trait2 <- traits[j]
      cat("\n🧩 Processing", trait1, "->", trait2, "\n")
      
      # File paths
      sig_trait1_file <- paste0(input_dir, trait1, "_snp_significant.tsv")
      clean_trait2_file <- paste0(input_dir, trait2, "_snp_clean.tsv")
      effect_clean_trait2_file <- paste0(input_dir, trait2, "_effect_snp_clean.tsv")
      file_name <- paste0(trait1, "_", trait2)
      
      # 1. Merge
      merged_data <- tryCatch({
        tidy_snps(sig_trait1_file, effect_clean_trait2_file, clean_trait2_file)
      }, error = function(e) {
        message("⚠️ Skipping due to error: ", e$message)
        return(NULL)
      })
      
      if (!is.null(merged_data) && nrow(merged_data) > 0) {
        
        cat("🔹 Merged SNPs:", nrow(merged_data), "\n")
        
        # 2. Significance labeling (not necessary but done for consistency)
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
        trait_matrix[i, j] <- num_sig_snps
        
      } else {
        cat("⚠️ No merged SNPs found for this pair.\n")
        trait_matrix[i, j] <- 0
      }
    }
  }
}

# Save matrix
matrix_file <- paste0(output_dir, "significant_snp_matrix.tsv")
fwrite(cbind(Trait = rownames(trait_matrix), as.data.frame(trait_matrix)), matrix_file, sep = "\t", row.names = FALSE)

cat("\n✅ Asymmetric SNP matrix saved to", matrix_file, "\n")
