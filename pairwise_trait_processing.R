#pairwise_trait_processing.R
#Processing all pairwise trait combinations

library(dplyr)
library(data.table)
library(tibble)

# Parameters
input_dir <- "preprocessed_data/"
output_dir <- "processed_data/"
distance_threshold <- 200000
significance_threshold <- 1e-6
skip_traits <- c("AD_7158", "ALCOH_3727", "APNEA_7362", "ASBEST_3282", "BRE_CANCER_174",
                 "COR_ARTERY_2314", "COVID_0403", "DENTAL_521.1", "DEPRES_6475",
                 "DIAB_GEST_6553", "DIAB_TYP1_3791")

# Setup
dir.create(output_dir, showWarnings = FALSE)
files <- list.files(input_dir, pattern = "_snp_significant.tsv$", full.names = FALSE)
traits <- gsub("_snp_significant.tsv", "", files)
traits <- sort(unique(traits))
traits <- traits[!grepl("_effect$", traits)]

# Matrix setup
matrix_file <- file.path(output_dir, "significant_snp_matrix.tsv")
if (!file.exists(matrix_file)) {
  trait_matrix <- matrix(0, nrow = length(traits), ncol = length(traits))
  rownames(trait_matrix) <- traits
  colnames(trait_matrix) <- traits
  fwrite(cbind(Trait = traits, as.data.frame(trait_matrix)), matrix_file, sep = "\t", row.names = FALSE)
}
trait_matrix <- fread(matrix_file) %>% column_to_rownames("Trait")

# Control file
control_file <- file.path(output_dir, "pairwise_validation_control.tsv")
if (!file.exists(control_file)) {
  fwrite(data.frame(), control_file, sep = "\t")
}

# Distance filter
filter_snps_by_distance <- function(data, min_distance = 200000) {
  data <- data %>% arrange(p_value)
  list_of_included <- list()
  for (S in seq_len(nrow(data))) {
    shall_I_include <- TRUE
    for (i in seq_along(list_of_included)) {
      if (list_of_included[[i]]$chr == data$chr[S] &&
          abs(data$pos[S] - list_of_included[[i]]$pos) < min_distance) {
        shall_I_include <- FALSE
        break
      }
    }
    if (shall_I_include) list_of_included <- append(list_of_included, list(data[S, ]))
  }
  bind_rows(list_of_included)
}

# Validation stats
compute_validations <- function(N, K, n, k, simulations = 10000) {
  p_hyper <- phyper(q = k - 1, m = K, n = N - K, k = n, lower.tail = FALSE)
  fisher_matrix <- matrix(c(k, K - k, n - k, N - K - n + k), nrow = 2)
  p_fisher <- fisher.test(fisher_matrix, alternative = "greater")$p.value
  sim_k <- replicate(simulations, {
    length(intersect(sample(N, K), sample(N, n)))
  })
  p_monte <- mean(sim_k >= k)
  p_obs <- k / N
  return(c(p_hyper, p_fisher, p_monte, p_obs))
}

# Main loop
for (trait1 in traits) {
  if (trait1 %in% skip_traits) {
    cat("# Skipping already processed trait:", trait1, "\n")
    next
  }
  
  sig_trait1_file <- file.path(input_dir, paste0(trait1, "_snp_significant.tsv"))
  sig_trait1 <- fread(sig_trait1_file)
  if (nrow(sig_trait1) == 0) {
    cat("# Skipping base trait due to no significant SNPs:", trait1, "\n")
    next
  }
  
  sig_trait1 <- sig_trait1 %>%
    mutate(chr = as.character(chr), pos = as.numeric(pos),
           effect_allele = toupper(effect_allele), alt_allele = toupper(alt_allele)) %>%
    distinct(chr, pos, .keep_all = TRUE)  # deduplicate
  
  cat("\n# Processing base trait:", trait1, "\n")
  
  for (trait2 in traits) {
    if (trait1 == trait2) next
    
    cat("##", trait1, " → ", trait2, "\n")
    clean_file <- file.path(input_dir, paste0(trait2, "_snp_clean.tsv"))
    effect_file <- file.path(input_dir, paste0(trait2, "_effect_snp_clean.tsv"))
    output_filtered <- file.path(output_dir, paste0(trait1, "_", trait2, "_filtered_snps.tsv"))
    
    # Load and preprocess trait2 data
    effect2 <- fread(effect_file) %>%
      mutate(chr = as.character(chr), pos = as.numeric(pos), p_value = as.numeric(p_value),
             beta = as.numeric(beta), SE = as.numeric(SE)) %>%
      distinct(chr, pos, .keep_all = TRUE)
    
    clean2 <- fread(clean_file) %>%
      mutate(chr = as.character(chr), pos = as.numeric(pos),
             effect_allele = toupper(effect_allele), alt_allele = toupper(alt_allele)) %>%
      distinct(chr, pos, .keep_all = TRUE)
    
    merged <- tryCatch({
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
    
    if (nrow(merged) == 0) {
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
    
    cat("· Merged SNPs:", nrow(merged), "\n")
    
    filtered <- filter_snps_by_distance(merged, distance_threshold)
    sig_filtered <- filtered %>% filter(p_value < significance_threshold)
    cat("· Filtered SNPs:", nrow(filtered), " | Significant:", nrow(sig_filtered), "\n")
    
    fwrite(filtered, output_filtered, sep = "\t")
    trait_matrix[trait1, trait2] <- nrow(sig_filtered)
    
    validations <- compute_validations(
      N = nrow(effect2),
      K = nrow(effect2 %>% filter(p_value < significance_threshold)),
      n = nrow(filtered),
      k = nrow(sig_filtered)
    )
    
    control_row <- data.frame(
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


#source("pairwise_trait_processing.R")
