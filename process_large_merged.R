
# 0. Load libraries
library(dplyr)
library(tibble)
library(data.table)

# 1. Distance-based SNP filtering function
filter_snps_by_distance <- function(data, min_distance = 200000) {
  data <- data %>% arrange(p_value)
  list_of_included <- list()
  message("Filtering SNPs that are at least ", min_distance, " base pairs apart.")
  
  for (S in seq_len(nrow(data))) { 
    shall_I_include <- TRUE
    for (i in seq_along(list_of_included)) {
      if (list_of_included[[i]]$chr == data$chr[S] &&
          abs(data$pos[S] - list_of_included[[i]]$pos) < min_distance) {
        shall_I_include <- FALSE
        break
      }
    }
    if (shall_I_include) {
      list_of_included <- append(list_of_included, list(data[S, ]))
    }
  }
  return(bind_rows(list_of_included))
}

# 2. Main processing function
process_SNPs_by_chromosomes <- function(sig_trait1_file, clean_trait2_file, effect_clean_trait2_file, distance_threshold, file_name) {
  dir.create("processed_data", showWarnings = FALSE)
  
  message("Loading data files...")
  sig_trait1 <- as_tibble(fread(sig_trait1_file))
  effect_trait2 <- as_tibble(fread(effect_clean_trait2_file))
  clean_trait2 <- as_tibble(fread(clean_trait2_file))
  
  # Standardize formatting
  sig_trait1 <- sig_trait1 %>% 
    mutate(chr = as.character(chr), pos = as.numeric(pos),
           effect_allele = toupper(effect_allele), alt_allele = toupper(alt_allele))
  
  effect_trait2 <- effect_trait2 %>%
    mutate(p_value = as.numeric(p_value), chr = as.character(chr),
           pos = as.numeric(pos), beta = as.numeric(beta), SE = as.numeric(SE))
  
  clean_trait2 <- clean_trait2 %>%
    mutate(effect_allele = toupper(effect_allele), alt_allele = toupper(alt_allele),
           chr = as.character(chr), pos = as.numeric(pos))
  
  # Merge datasets
  message("Merging datasets...")
  merged_data <- sig_trait1 %>%
    inner_join(effect_trait2, by = c("chr", "pos")) %>%
    inner_join(clean_trait2, by = c("chr", "pos"), suffix = c("_sig", "_clean")) %>%
    mutate(
      final_rs_id = coalesce(rs_id.x, rs_id.y),
      effect_allele = coalesce(effect_allele_sig, effect_allele_clean),
      alt_allele = coalesce(alt_allele_sig, alt_allele_clean)
    ) %>%
    select(chr, pos, effect_allele, alt_allele, final_rs_id, p_value, beta, SE) %>%
    distinct()  # Remove potential duplicates
  
  message("Number of merged SNPs: ", nrow(merged_data))
  
  # Chromosome-wise filtering
  message("Filtering SNPs by chromosome...")
  all_filtered_snps <- list()
  for (chrom in unique(merged_data$chr)) {
    message("Processing chromosome: ", chrom)
    chr_data <- merged_data %>% filter(chr == chrom) %>% arrange(p_value)
    message("SNPs merged: ", nrow(chr_data))
    chr_filtered <- filter_snps_by_distance(chr_data, distance_threshold)
    message("SNPs filtered: ", nrow(chr_filtered))
    all_filtered_snps[[chrom]] <- chr_filtered
  }
  
  # Combine filtered SNPs across chromosomes
  combined_snps <- bind_rows(all_filtered_snps)
  message("Combined SNPs after per-chromosome filtering: ", nrow(combined_snps))
  
  # Global filtering across chromosomes
  final_filtered <- filter_snps_by_distance(combined_snps, distance_threshold)
  message("Final SNP count after global filtering: ", nrow(final_filtered))
  
  # Final significance filtering (optional before saving)
  significant_final <- final_filtered %>% filter(p_value < 1e-5)
  message("Number of significant SNPs: ", nrow(significant_final))
  
  # Save result
  output_path <- paste0("processed_data/", file_name, "_filtered_snps.tsv")
  fwrite(final_filtered, output_path, sep = "\t")
  message("Filtered SNPs saved to: ", output_path)
}

# 3. Example run
process_SNPs_by_chromosomes(
  sig_trait1_file = "preprocessed_data/MENOP_2553_snp_significant.tsv",
  clean_trait2_file = "preprocessed_data/APNEA_7362_snp_clean.tsv",
  effect_clean_trait2_file = "preprocessed_data/APNEA_7362_effect_snp_clean.tsv",
  distance_threshold = 200000,
  file_name = "MENOP1_APNEA2_Test"
)
