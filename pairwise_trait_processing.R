#pairwise_trait_processing.R
#Processing all pairwise trait combinations

#1. Load libraries
library(dplyr)
library(data.table)
library(tibble)

#2. Set parameters and directories
##Define parameters: input and output directories, distance threshold, significance threshold, and traits to skip
input_dir <- "preprocessed_data/"
output_dir <- "processed_data/"
distance_threshold <- 200000
significance_threshold <- 1e-6 #In this project we have used 1e-5, 1e-6 and 1e-8
skip_traits <- c(
  #"AD_7158", "ALCOH_3727", "APNEA_7362", "ASBEST_3282", "BRE_CANCER_174", <- Here you can add traits to skip if already processed
  )

##Setup: identify input files and traits, and create the documents to fill later
dir.create(output_dir, showWarnings = FALSE)
files <- list.files(input_dir, pattern = "_snp_significant.tsv$", full.names = FALSE) #Select all the files from the input directory that match the pattern
traits <- gsub("_snp_significant.tsv", "", files) #Extract trait names from the files
traits <- sort(unique(traits)) #Make sure the trait names are unique and sorted
traits <- traits[!grepl("_effect$", traits)] #Exclude traits ending with "_effect" to avoid duplicates

##Matrix setup
#Create a matrix to store the number of significant SNPs for each merged trait pair
matrix_file <- file.path(output_dir, "significant_snp_matrix.tsv") #Create a file path for the output matrix
if (!file.exists(matrix_file)) {
  trait_matrix <- matrix(0, nrow = length(traits), ncol = length(traits))
  rownames(trait_matrix) <- traits #Set row names of the matrix to traits
  colnames(trait_matrix) <- traits #Set column names of the matrix to traits
  fwrite(cbind(Trait = traits, as.data.frame(trait_matrix)), matrix_file, sep = "\t", row.names = FALSE)
}
trait_matrix <- fread(matrix_file) %>% column_to_rownames("Trait")

##Control log file: to track the processing of trait pairs
#Initialize the file if it does not already exist
control_file <- file.path(output_dir, "pairwise_validation_control.tsv")
if (!file.exists(control_file)) {
  fwrite(data.frame(), control_file, sep = "\t")
}

#3. Distance filter; Linkage disequilibrium (LD) prunning to avoid bias
filter_snps_by_distance <- function(data, min_distance = 200000) {
  data <- data %>% arrange(p_value) #Sort by p_value to prioritize significant SNPs
  list_of_included <- list() #Initialize list of candidate variants
  for (S in seq_len(nrow(data))) { #Iterate through each SNP
    shall_I_include <- TRUE #Inicialize switch
    for (i in seq_along(list_of_included)) {
      if (list_of_included[[i]]$chr == data$chr[S] && 
          abs(data$pos[S] - list_of_included[[i]]$pos) < min_distance) { #Check if SNP is within the distance threshold in the same chromosome
        shall_I_include <- FALSE #change switch: therefore exclude SNP from candidates
        break
      }
    }
    if (shall_I_include) list_of_included <- append(list_of_included, list(data[S, ]))
  }
  bind_rows(list_of_included) #Return the filtered SNPs as a data frame
}

#4. Function to compute validation statistics
compute_validations <- function(N, K, n, k, simulations = 10000) {
    # Contingency Table (2x2) for Fisher's Exact Test:
  # ┌────────────────────────────┬────────────────────────────┐
  # │         In Trait2          │     Not in Trait2          │
  # ├────────────┬───────────────┼────────────┬───────────────┤
  # │ In Trait1  │       k       │ Not In T2  │    K - k      │
  # │            │  (overlap)    │ (only T1)  │               │
  # ├────────────┼───────────────┼────────────┼───────────────┤
  # │ Not in T1  │     n - k     │ Neither    │ N - K - n + k │
  # │ (only T2)  │               │            │               │
  # └────────────┴───────────────┴────────────┴───────────────┘
  #
  # Where:
  # N = total number of SNPs in universe
  # K = number of significant SNPs in Trait1
  # n = number of SNPs tested in Trait2 (after distance filtering)
  # k = number of overlapping significant SNPs between Trait1 and Trait2 
  p_hyper <- phyper(q = k - 1, m = K, n = N - K, k = n, lower.tail = FALSE) #Hypergeometric test: overlap significance
  fisher_matrix <- matrix(c(k, K - k, n - k, N - K - n + k), nrow = 2) #Create contingency table for Fisher's test
  p_fisher <- fisher.test(fisher_matrix, alternative = "greater")$p.value #Fisher's Exact Test: independence of overlap
  sim_k <- replicate(simulations, { #Monte Carlo simulation to empirically overlap probability
    length(intersect(sample(N, K), sample(N, n)))
  })
  p_monte <- mean(sim_k >= k)
  p_obs <- k / N  #observed proportion of overlap
  return(c(p_hyper, p_fisher, p_monte, p_obs))
}

#5. Main loop to process each trait pair
for (trait1 in traits) { #Iterate through each trait
  if (trait1 %in% skip_traits) {
    cat("# Skipping already processed trait:", trait1, "\n")
    next
  }
  
  sig_trait1_file <- file.path(input_dir, paste0(trait1, "_snp_significant.tsv"))
  sig_trait1 <- fread(sig_trait1_file)
  if (nrow(sig_trait1) == 0) { #If no significant SNPs found, skip this trait
    cat("# Skipping base trait due to no significant SNPs:", trait1, "\n")
    next
  }
  
  sig_trait1 <- sig_trait1 %>%
    mutate(chr = as.character(chr), pos = as.numeric(pos),
           effect_allele = toupper(effect_allele), alt_allele = toupper(alt_allele)) %>%
    distinct(chr, pos, .keep_all = TRUE)  #Fix the structure of the data
  
  cat("\n# Processing base trait:", trait1, "\n")
  
  for (trait2 in traits) { #Chose the second trait to compare with the first, iterate through all the possible traits
    if (trait1 == trait2) next #no self comparison
    
    ##Define file paths for trait2
    cat("##", trait1, " → ", trait2, "\n")
    clean_file <- file.path(input_dir, paste0(trait2, "_snp_clean.tsv"))
    effect_file <- file.path(input_dir, paste0(trait2, "_effect_snp_clean.tsv")) 
    output_filtered <- file.path(output_dir, paste0(trait1, "_", trait2, "_filtered_snps.tsv"))
    
    ##Load and preprocess trait2 data
    effect2 <- fread(effect_file) %>%
      mutate(chr = as.character(chr), pos = as.numeric(pos), p_value = as.numeric(p_value),
             beta = as.numeric(beta), SE = as.numeric(SE)) %>%
      distinct(chr, pos, .keep_all = TRUE)
    
    clean2 <- fread(clean_file) %>%
      mutate(chr = as.character(chr), pos = as.numeric(pos),
             effect_allele = toupper(effect_allele), alt_allele = toupper(alt_allele)) %>%
      distinct(chr, pos, .keep_all = TRUE)
    
    merged <- tryCatch({ ##MErge by chromosome and position
      sig_trait1 %>%
        inner_join(effect2, by = c("chr", "pos")) %>%
        inner_join(clean2, by = c("chr", "pos"), suffix = c("_sig", "_clean")) %>%
        mutate(
          final_rs_id = coalesce(rs_id.x, rs_id.y, rs_id),
          effect_allele = coalesce(effect_allele_sig, effect_allele_clean),
          alt_allele = coalesce(alt_allele_sig, alt_allele_clean)
        ) %>%
        select(chr, pos, final_rs_id, effect_allele, alt_allele, p_value, beta, SE)
    }, error = function(e) {
      message("⚠️ Failed to merge ", trait1, " -> ", trait2, ": ", e$message)
      return(data.frame())
    })
    
    if (nrow(merged) == 0) { #If no overlapping SNPs found, log 0 and continues with other traits skipping this one
      cat("No overlapping SNPs between", trait1, "and", trait2, "\n")
      trait_matrix[trait1, trait2] <- 0
      control_row <- data.frame(
        trait1 = trait1, trait2 = trait2,
        merged_data_size = 0, filtered_data_size = 0, significant_filtered = 0,
        p_hyper = NA, p_fisher = NA, p_monte = NA, observed_p = NA
      )
      fwrite(control_row, control_file, sep = "\t", append = TRUE)
      fwrite(cbind(Trait = rownames(trait_matrix), as.data.frame(trait_matrix)),
             matrix_file, sep = "\t", row.names = FALSE)
      next
    }
    
    ##Filter, validate and log the results of the overlap
    cat("· Merged SNPs:", nrow(merged), "\n")
    
    filtered <- filter_snps_by_distance(merged, distance_threshold) #filtering
    sig_filtered <- filtered %>% filter(p_value < significance_threshold) #select significant
    cat("· Filtered SNPs:", nrow(filtered), " | Significant:", nrow(sig_filtered), "\n")
    
    fwrite(filtered, output_filtered, sep = "\t") #save filtered files
    trait_matrix[trait1, trait2] <- nrow(sig_filtered) #update matrix iteratibly
    
    validations <- compute_validations( #perform validation of the overlap
      N = nrow(effect2),
      K = nrow(effect2 %>% filter(p_value < significance_threshold)),
      n = nrow(filtered),
      k = nrow(sig_filtered)
    )
    
    control_row <- data.frame( #save log iteratively
      trait1 = trait1, trait2 = trait2,
      merged_data_size = nrow(merged),
      filtered_data_size = nrow(filtered),
      significant_filtered = nrow(sig_filtered),
      p_hyper = validations[1],
      p_fisher = validations[2],
      p_monte = validations[3],
      observed_p = validations[4]
    )
    fwrite(control_row, control_file, sep = "\t", append = TRUE)
    
    fwrite(cbind(Trait = rownames(trait_matrix), as.data.frame(trait_matrix)),
           matrix_file, sep = "\t", row.names = FALSE)
  }
  
  cat("Done:", trait1, "\n")
}

#6. run this command to execute the script
#source("pairwise_trait_processing.R") 
