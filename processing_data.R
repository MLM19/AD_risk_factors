#Processing code up to creating table that relates how significant SNPs are between traits

#0. Load libraries
library(dplyr)
library(data.table)

#1. Find significant SNPs
find_significant_snps <- function(sig_data_trait1, clean_data_trait2) {
  message("Starting search for significant SNPs in another trait's data...")
  
  sig_data_trait1 <- sig_data_trait1 %>%
    mutate(chr = as.character(chr),
           pos = as.numeric(pos),
           rs_id = as.character(rs_id)) %>%
    distinct(chr, pos, .keep_all = TRUE)
  
  clean_data_trait2 <- clean_data_trait2 %>%
    mutate(chr = as.character(chr),
           pos = as.numeric(pos),
           rs_id = as.character(rs_id)) %>%
    distinct(chr, pos, effect_allele, alt_allele, .keep_all = TRUE)
  
  message("Number of unique significant SNPs from trait 1: ", nrow(sig_data_trait1))
  message("Dimension of clean data from trait 2:", paste(dim(clean_data_trait2), collapse = " x "))
  
  found_snps <- clean_data_trait2 %>%
    inner_join(sig_data_trait1 %>% select(chr, pos, rs_id), by = c("chr", "pos"))
  
  message("Initial number of SNPs found in trait 2 based on position and chromosome: ", nrow(found_snps))
  message("Columns in found_snps: ", paste(colnames(found_snps), collapse = ", "))
  
  # Decide which rs_id to keep
  final_table <- found_snps %>%
    mutate(
      final_rs_id = ifelse(rs_id.x == rs_id.y, rs_id.x, 
                           ifelse(nrow(sig_data_trait1) > nrow(clean_data_trait2), rs_id.x, rs_id.y))
    ) %>%
    select(chr, pos, effect_allele, alt_allele, final_rs_id)
  
  message("Number of SNPs after processing: ", nrow(final_table))
  message("Columns in final_table: ", paste(colnames(final_table), collapse = ", "))
  message("Search complete.")
  return(final_table)
}


#2. Test significance function
test_significance <- function(result_table, effect_data_trait2, significance_threshold = 10e-5) {
  message("Assessing significance of SNPs in trait 2...")
  
  if (nrow(result_table) == 0) {
    message("No SNPs to assess for significance.")
    return(result_table)
  }
  
  effect_data_trait2 <- effect_data_trait2 %>%
    mutate(chr = as.character(chr),
           pos = as.numeric(pos),
           rs_id = as.character(rs_id))
  
  message("Number of significant SNPs in trait 2: ", sum(effect_data_trait2$p_value < significance_threshold))
  
  merged_data <- result_table %>%
    left_join(effect_data_trait2 %>%
                select(chr, pos, p_value, beta, SE, rs_id),
              by = c("chr", "pos"))
  
  if (nrow(merged_data) == 0) {
    warning("No overlapping SNPs between traits.")
    return(merged_data)
  }
  
  merged_result <- merged_data %>%
    mutate(
      significance_level = case_when(
        p_value < significance_threshold ~ "Significant",
        p_value < 0.05 ~ "Potentially Significant",
        TRUE ~ "Not Significant"
      )
    ) %>%
    select(chr, pos, effect_allele, alt_allele, final_rs_id, p_value, beta, SE, significance_level)
  
  message("Significance assessment complete.")
  return(merged_result)
}


##Function to save the results
save_result_table <- function(result_table, file_name) {
  #Create the processed_data directory if it doesn't exist
  dir.create("processed_data", showWarnings = FALSE)
  
  # Create the full file path
  file_path <- paste0("processed_data/", file_name, "_merged_table.tsv")
  
  # Save the result table
  fwrite(result_table, file_path, sep = "\t")
  
  message("Result table saved as: ", file_path)
}

##3. Main Script ##

##1. Load the files
snp_sig_trait1_file <- "preprocessed_data/AD_7158_snp_significant.tsv"
effect_data_trait2_file <- "preprocessed_data/ALCOH_3727_effect_snp_clean.tsv"
clean_data_trait2_file <- "preprocessed_data/ALCOH_3727_snp_clean.tsv"

snp_significant <- fread(snp_sig_trait1_file, sep = "\t")
effect_data_trait2 <- fread(effect_data_trait2_file, sep = "\t")
clean_data_trait2 <- fread(clean_data_trait2_file, sep = "\t")

##2. Run the analysis
result_table <- find_significant_snps(snp_significant, clean_data_trait2)

##3. Assess significance
result_table <- test_significance(result_table, effect_data_trait2)

##4. Print the results
print(result_table)
significance_freq <- table(result_table$significance_level)
significance_df <- as.data.frame(significance_freq)
colnames(significance_df) <- c("Significance Level", "Frequency")
print(significance_df)

##5. Save result in the directory
file_name <- "AD1_ALCOH2"
save_result_table(result_table, file_name)
