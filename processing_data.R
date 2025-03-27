#processing code: automatic merging + snp enrichment, changed distance_filtering

#0. Load libraries
library(dplyr)
library(tibble)
library(data.table)

#1. Merging: Function to find significant SNPs in trait1 and match with trait2
tidy_snps <- function(sig_trait1_file, effect_clean_trait2_file, clean_trait2_file) {
  message("Loading data files...")
  sig_trait1 <- as_tibble(fread(sig_trait1_file)) #making sure that the format is data frames (tibble)
  effect_trait2 <- as_tibble(fread(effect_clean_trait2_file))
  clean_trait2 <- as_tibble(fread(clean_trait2_file))
  
  message("Checking column names for debugging:")
  print(colnames(sig_trait1))
  print(colnames(effect_trait2))
  print(colnames(clean_trait2))
  
  #transformations before merging
  sig_trait1 <- sig_trait1 %>%
    mutate(
      chr = as.character(chr),   #debugging
      pos = as.numeric(pos),
      effect_allele = toupper(effect_allele),
      alt_allele = toupper(alt_allele),
    )
  
  effect_trait2 <- effect_trait2 %>%
    mutate(
      p_value = as.numeric(p_value),
      chr = as.character(chr),   #debugging
      pos = as.numeric(pos),
      beta = as.numeric(beta),
      SE = as.numeric(SE)
    )
  
  clean_trait2 <- clean_trait2 %>%
    mutate(
      effect_allele = toupper(effect_allele),
      alt_allele = toupper(alt_allele),
      chr = as.character(chr),   #debugging
      pos = as.numeric(pos),
    )
  
  message("DEBUG the unique values in rs_id in sig_trait1, effect_trait2 and clean_trait2:")
  print(length(unique(sig_trait1$rs_id)))
  print(length(unique(effect_trait2$rs_id)))
  print(length(unique(clean_trait2$rs_id)))
  
  pre_merge <- inner_join(sig_trait1, effect_trait2, by = c("chr", "pos"))
  message("Number of SNPs match before filtering..")
  print(nrow(pre_merge)) 
  message("Number of unmatched SNPs: ")
  print(nrow(anti_join(sig_trait1, effect_trait2, by = c("chr", "pos"))))
  
  message("Merging the significant SNPs from trait1 to trait2...")
  #merging the data and selecting only the relevant columns
  merged_data <- sig_trait1 %>%
    inner_join(effect_trait2, by = c("chr", "pos")) %>%
    inner_join(clean_trait2, by = c("chr", "pos"), suffix = c("_sig", "_clean")) %>%
    mutate(
      #Prioritize given rs_id (which follows "rs123456" format) over manually created ones
      final_rs_id = case_when(
        grepl("^rs[0-9]+$", rs_id.x) ~ rs_id.x,   #If rs_id.x is in "rs123456" format, use it
        grepl("^rs[0-9]+$", rs_id.y) ~ rs_id.y,   #Else, if rs_id.y is valid, use it
        grepl("^rs[0-9]+$", rs_id) ~ rs_id,       #Else, if rs_id from clean_trait2 is valid, use it
        TRUE ~ coalesce(rs_id.x, rs_id.y, rs_id)  #If none are valid, use whatever is available
      ),
      effect_allele = coalesce(effect_allele_sig, effect_allele_clean),
      alt_allele = coalesce(alt_allele_sig, alt_allele_clean)
    ) %>%
    dplyr::select(chr, pos, effect_allele, alt_allele, final_rs_id, p_value, beta, SE)
  
  message("DEBUG")
  pre_merge <- inner_join(sig_trait1, effect_trait2, by = c("chr", "pos"))
  message("Number of SNPs match before filtering..")
  print(nrow(pre_merge)) 
  message("Number of unmatched SNPs: ")
  print(nrow(anti_join(sig_trait1, effect_trait2, by = c("chr", "pos"))))
  
  message("Number of merged SNPs: ", nrow(merged_data))
  return(merged_data)
}

#2. Test significance function #Oscar said not really needed bc we can see it by the p_value
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
  data <- data %>% arrange(p_value)
  
  #new list of included SNPs
  list_of_included <- list()
  
  message("Filtering SNPs that are at least ", min_distance, " base pairs apart.") #marker
  
  #Iterate over SNPs in data
  for (S in seq_len(nrow(data))) { 
    shall_I_include <- TRUE  # Switch variable
    
    #Check against list of previously included SNPs
    for (i in seq_along(list_of_included)) {
      if (list_of_included[[i]]$chr == data$chr[S] &&
          abs(data$pos[S] - list_of_included[[i]]$pos) < min_distance) {
        shall_I_include <- FALSE
        break
      }
    }
    
    #If SNP is at a valid distance, add to the list
    if (shall_I_include) {
      list_of_included <- append(list_of_included, list(data[S, ]))
    }
  }
  
  #Convert list to dataframe
  return(bind_rows(list_of_included))
}

#4. Calculate enrichment: probability of SNPs being significant by chance and compare it to the observed probability
calculate_snp_enrichment <- function(filtered_result_table, all_snps_file, num_simulations = 1000, pval_threshold = 1e-5) {
  message("Loading all SNPs dataset...")
  all_snps_data <- fread(all_snps_file, sep = "\t") %>%
    filter(!is.na(p_value))  # Remove NA values for consistency
  
  message("Calculating observed probability of significance...")
  observed_significant_count <- sum(filtered_result_table$p_value < pval_threshold, na.rm = TRUE)
  observed_prob <- observed_significant_count / nrow(filtered_result_table)
  message("Observed significant SNPs: ", observed_significant_count, " out of ", nrow(filtered_result_table))
  
  #If no significant SNPs found, return NA
  if (nrow(filtered_result_table) == 0 || observed_prob == 0) {
    warning("No significant SNPs found in the observed set. Enrichment cannot be calculated.")
    return(data.frame(observed_probability = observed_prob, expected_probability = NA, enrichment_fold = NA, z_score = NA, p_value = NA))
  }
  
  message("Running ", num_simulations, " simulations to estimate expected probability...")
  
  #Simulated probability distribution
  simulated_probs <- replicate(num_simulations, {
    random_snps <- all_snps_data %>% sample_n(nrow(filtered_result_table), replace = FALSE)
    sum(random_snps$p_value < pval_threshold, na.rm = TRUE) / nrow(random_snps)
  })
  
  #Expected probability under random sampling
  mean_simulated_prob <- mean(na.omit(simulated_probs))
  
  #Enrichment fold change
  enrichment_fold <- ifelse(mean_simulated_prob == 0, NA, observed_prob / mean_simulated_prob)
  
  #Statistical significance: Z-score & p-value
  z_score <- (observed_prob - mean_simulated_prob) / sd(na.omit(simulated_probs))
  p_value <- pnorm(-abs(z_score)) * 2  #Two-tailed test
  
  #Store results in a data frame
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
  
  #Convert enrichment results to a dataframe for saving
  enrichment_results_df <- as.data.frame(t(enrichment_results))
  colnames(enrichment_results_df) <- c("Value")
  
  #Save enrichment results
  enrichment_output_path <- paste0("processed_data/", file_name, "_enriched_table.tsv")
  fwrite(enrichment_results_df, enrichment_output_path, sep = "\t")
  
  message("Enrichment results saved to: ", enrichment_output_path)
}


##6. Main Script Execution ##

#Define parameters and call processing function
process_SNPs(
  sig_trait1_file <- "preprocessed_data/AD_7158_snp_significant.tsv",
  clean_trait2_file <- "preprocessed_data/DIAB_GEST_6553_snp_clean.tsv",
  effect_clean_trait2_file <- "preprocessed_data/DIAB_GEST_6553_effect_snp_clean.tsv",
  distance_threshold <- 200000,  
  file_name <- "AD1_DIAB_GEST2_Test2"
)

