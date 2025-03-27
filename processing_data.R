#processing code: automatic merging + snp enrichment, this results are not good

#0. Load libraries
library(dplyr)
library(tibble)
library(data.table)

#1. Merging: Function to find significant SNPs in trait1 and martch with trait2
tidy_snps <- function(sig_trait1_file, effect_clean_trait2_file, clean_trait2_file) {
  message("Loading data files...")
  sig_trait1 <- as_tibble(fread(sig_trait1_file))
  effect_trait2 <- as_tibble(fread(effect_clean_trait2_file))
  clean_trait2 <- as_tibble(fread(clean_trait2_file))
  
  message("Checking column names for debugging:")
  print(colnames(sig_trait1))
  print(colnames(effect_trait2))
  print(colnames(clean_trait2))
  
  message("Merging the significant SNPs from trait1 to trait2...")
  merged_data <- sig_trait1 %>%
    inner_join(effect_trait2, by = c("chr", "pos")) %>%
    inner_join(clean_trait2, by = c("chr", "pos"), suffix = c("_sig", "_clean")) %>%
    mutate(
      final_rs_id = coalesce(rs_id.x, rs_id.y, rs_id),
      effect_allele = coalesce(effect_allele_sig, effect_allele_clean),  
      alt_allele = coalesce(alt_allele_sig, alt_allele_clean)
    ) %>%
    dplyr::select(chr, pos, effect_allele, alt_allele, final_rs_id, p_value, beta, SE)
  
  message("Number of merged SNPs: ", nrow(merged_data))
  return(merged_data)
}

#2. Test significance function
test_significance <- function(data, significance_threshold = 1e-5) {
  message("Assessing significance of SNPs in trait 2...")
  data %>%
    mutate(
      significance_level = case_when(
        is.na(p_value) ~ "Unknown",
        p_value < significance_threshold ~ "Significant",
        p_value < 0.05 ~ "Potentially Significant",
        TRUE ~ "Not Significant"
      )
    )
}

##3. Filter by distance
#filter SNPs based on the minimum distance threshold
filter_snps_by_distance <- function(data, min_distance = 200000) {
  #Sorting by p_value from lowest to highest, then chr and position
  data <- data %>% arrange(p_value, chr, pos)
  
  #Initialize filtered list
  filtered_snps <- list()
  prev_pos <- -Inf #to ensure that the first SNP is also counted in the process.
  
  message("Filtering SNPs to find those that are further than: ", min_distance, " base pairs.") #marker
  
  #Go through the data and keep as signal those SNPs that are further away from minimal distance (using list to filter)
  for (i in seq_len(nrow(data))) {
    if (data$pos[i] - prev_pos >= min_distance) {
      filtered_snps <- append(filtered_snps, list(data[i, ]))
      prev_pos <- data$pos[i]
    }
  }
  
  #Convert list into dataframe 
  return(bind_rows(filtered_snps))
}

#4. Calculate enrichment: probability of SNPs being significant by chance and compare it to the observed probability
calculate_snp_enrichment <- function(filtered_result_table, all_snps_file, num_simulations = 1000, pval_threshold = 1e-5) {
  #Load all SNPs dataset as a background file
  all_snps_data <- fread(all_snps_file,  sep = "\t")
  
  #Filter out NA values for p_value calculations (there should be none)
  filtered_result_table <- filtered_result_table %>% filter(!is.na(p_value))
  all_snps_data <- all_snps_data %>% filter(!is.na(p_value))
  
  #Calculate observed probability of significance
  observed_significant_count <- sum(filtered_result_table$p_value < pval_threshold, na.rm = TRUE)
  observed_prob <- observed_significant_count / nrow(filtered_result_table)
  message("The number of observed significant SNPs is ", observed_significant_count)

  
  #Simulate random selection of SNPs to compute expected probability
  simulated_probs <- replicate(num_simulations, {         #replicate: store output of each iteration in a vector
    random_snps <- all_snps_data %>% sample_n(min(nrow(filtered_result_table), nrow(all_snps_data)), replace = FALSE)  #random sampling of SNPs to simulate a null distribution
    sum(random_snps$p_value < pval_threshold, na.rm = TRUE) / nrow(random_snps)  #compute proportion of randomly selected SNPs that meet significance threshold
  })
  
  #Compute expected probability
  mean_simulated_prob <- mean(na.omit(simulated_probs))
  
  #Calculate enrichment statistics
  enrichment_fold <- ifelse(mean_simulated_prob == 0, NA, observed_prob / mean_simulated_prob)
  z_score <- (observed_prob - mean_simulated_prob) / sd(na.omit(simulated_probs))
  p_value <- pnorm(-abs(z_score)) * 2
  
  #Convert enrichment results into a data frame
  enrichment_results_df <- data.frame(
    observed_probability = observed_prob,
    expected_probability = mean_simulated_prob,
    enrichment_fold = enrichment_fold,
    z_score = z_score,
    p_value = p_value
  )
  
  print(enrichment_results_df)
  return(enrichment_results_df)
}
  
##5. Main function
process_SNPs <- function(sig_trait1_file, clean_trait2_file, effect_clean_trait2_file, distance_threshold, file_name){
  #Make sure that the directory exists, otherwise create it 
  dir.create("processed_data", showWarnings = FALSE)
  
  #Find significant SNPs in trait 2 based on trait 1
  result_table <- tidy_snps(sig_trait1_file, effect_clean_trait2_file, clean_trait2_file)
  
  #Assess significance of found SNPs in trait 2
  result_table <- test_significance(result_table, effect_clean_trait2_file)
  
  #Filter SNPs based on proximity to the most significant ones
  filtered_result_table <- filter_snps_by_distance(result_table, distance_threshold)
  
  message("Number of SNPs after filtering: ", nrow(filtered_result_table))
  
  #Save filtered results to a dynamically named file
  filtered_output_path <- paste0("processed_data/", file_name, "_filtered_snps.tsv")
  fwrite(filtered_result_table, filtered_output_path, sep = "\t")
  
  message("Filtered SNPs saved to: ", filtered_output_path)
  
  #Perform enrichment analysis
  enrichment_results <- calculate_snp_enrichment(
    filtered_result_table,
    all_snps_file = effect_clean_trait2_file, 
    num_simulations = 10000
  )
  
  # Convert enrichment results to a dataframe for saving
  enrichment_results_df <- as.data.frame(t(enrichment_results))
  colnames(enrichment_results_df) <- c("Value")
  
  # Save enrichment results to a dynamically named file
  enrichment_output_path <- paste0("processed_data/", file_name, "_enriched_table.tsv")
  fwrite(enrichment_results_df, enrichment_output_path, sep = "\t")
  
  message("Enrichment results saved to: ", enrichment_output_path)
}


##6. Main Script Execution ##

#Define parameters and call processing function
process_SNPs(
  sig_trait1_file <- "preprocessed_data/AD_7158_snp_significant.tsv",
  clean_trait2_file <- "preprocessed_data/COVID_0403_snp_clean.tsv",
  effect_clean_trait2_file <- "preprocessed_data/COVID_0403_effect_snp_clean.tsv",
  distance_threshold <- 200000,  # Reduce from 200000 to 50000 ?
  file_name <- "AD1_COVID2_Test"
)

