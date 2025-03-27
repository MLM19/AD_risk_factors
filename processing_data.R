#processing code: automatic merging + no enrichment since validation after, changed distance_filtering

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
      chr = as.factor(chr),
      pos = as.numeric(pos),
      effect_allele = toupper(effect_allele),
      alt_allele = toupper(alt_allele),
    )
  
  effect_trait2 <- effect_trait2 %>%
    mutate(
      p_value = as.numeric(p_value),
      chr = as.factor(chr),
      pos = as.numeric(pos),
      beta = as.numeric(beta),
      SE = as.numeric(SE)
    )
  
  clean_trait2 <- clean_trait2 %>%
    mutate(
      effect_allele = toupper(effect_allele),
      alt_allele = toupper(alt_allele),
      chr = as.factor(chr),
      pos = as.numeric(pos),
    )
  
  message("Merging the significant SNPs from trait1 to trait2...")
  #merging the data and selecting only the relevant columns
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
  
  #Count significant SNPs
  num_sig_snps_filtered <- sum(filtered_result_table$p_value < 1e-5, na.rm = TRUE)
  message("Number of significant SNPs after filtering: ", num_sig_snps_filtered, " out of ", nrow(filtered_result_table), " total filtered SNPs.")
  
  #Save filtered results to a dynamically named file
  filtered_output_path <- paste0("processed_data/", file_name, "_filtered_snps.tsv")
  fwrite(filtered_result_table, filtered_output_path, sep = "\t")
  
  message("Filtered SNPs saved to: ", filtered_output_path)
  
  
  
}

##6. Main Script Execution ##

#Define parameters and call processing function
process_SNPs(
  sig_trait1_file <- "preprocessed_data/AD_7158_snp_significant.tsv",
  clean_trait2_file <- "preprocessed_data/STROK_4723_snp_clean.tsv",
  effect_clean_trait2_file <- "preprocessed_data/STROK_4723_effect_snp_clean.tsv",
  distance_threshold <- 200000,  
  file_name <- "AD1_STROK2_Test"
)

