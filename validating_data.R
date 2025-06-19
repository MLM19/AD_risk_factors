##1. Load required libraries
library(dplyr)
library(ggplot2)

##2. Functions
### Function to compute N, K, k, and n 
compute_snp_counts <- function(trait_clean_file, trait_significant_file, ad_trait_file) {
  cat("Loading data files...\n")
  
  trait_clean <- read.delim(trait_clean_file, header = TRUE)
  trait_significant <- read.delim(trait_significant_file, header = TRUE)
  ad_trait <- read.delim(ad_trait_file, header = TRUE)
  
  if (!"rs_id" %in% colnames(trait_significant)) stop("Missing 'rs_id' in trait_significant")
  if (!"final_rs_id" %in% colnames(ad_trait)) stop("Missing 'final_rs_id' in ad_trait")
  if (!"p_value" %in% colnames(ad_trait)) stop("Missing 'p_value' in ad_trait")
  
  cat("Computing SNP counts...\n")
  
  ad_trait <- ad_trait[!is.na(ad_trait$p_value), ]
  
  N <- nrow(trait_clean)
  n <- nrow(trait_significant)
  K <- nrow(ad_trait)
  k <- sum(ad_trait$p_value < 1e-5)
  shared_rs_ids <- intersect(trait_significant$rs_id, ad_trait$final_rs_id)
  
  cat("Counts: K =", K, ", N =", N, ", k =", k, ", n =", n, "\n")
  
  return(list(
    N = N, K = K, k = k, n = n,
    trait_clean = trait_clean,
    trait_significant = trait_significant,
    ad_trait = ad_trait,
    shared_rs_ids = shared_rs_ids
  ))
}


### Main validation function with real p-values + plots
validate_snps <- function(trait_clean_file, trait_significant_file, ad_trait_file, trait_name = "Trait") {
  cat("Starting validation process for:", trait_name, "\n")
  
  counts <- compute_snp_counts(trait_clean_file, trait_significant_file, ad_trait_file)
  
  with(counts, {
    cat("Significant Shared SNPs (SS):", k, "\n")
    
    # --- Hypergeometric test (dhyper) ---
    hyper_p <- dhyper(k, n, N - n, K)
    cat(sprintf("Hypergeometric test p-value: %.4g\n", hyper_p))
    
    if (k == 0) {
      cat("Interpretation: No significant SNPs shared — possible weak or no association.\n\n")
    } else if (hyper_p < 0.05) {
      cat("Interpretation: Enrichment in shared SNPs is significant (unlikely due to chance).\n\n")
    } else {
      cat("Interpretation: No statistically significant enrichment of shared SNPs.\n\n")
    }
    
    # --- Fisher's Exact Test ---
    contingency_table <- matrix(c(n, N-n, k, K-k), nrow = 2)
    print(contingency_table)
    
    # ---   |    trait    |   AD>trait    ---
    #  SS   |      n      |      k        ---
    #  NSS  |     N-n     |     K-k       ---
    # TOTAL |      N      |      K        ---
    
    fisher_result <- fisher.test(contingency_table)
    cat(sprintf("Fisher's exact test p-value: %.4g\n", fisher_result$p.value))
    cat(sprintf("Fisher's exact test odds ratio: %.3f\n", fisher_result$estimate))
    print(fisher_result)
    
    if (fisher_result$p.value < 0.05) {
      cat("Interpretation: Fisher test supports non-random overlap between SNPs.\n")
    } else {
      cat("Interpretation: Overlap likely due to chance.\n")
    }
    
    if (fisher_result$estimate > 1) {
      cat("Odds ratio > 1: Enrichment — shared SNPs are more likely than expected.\n\n")
    } else if (fisher_result$estimate < 1) {
      cat("Odds ratio < 1: Depletion — shared SNPs are less likely than expected.\n\n")
    } else {
      cat("Odds ratio ≈ 1: No enrichment or depletion detected.\n\n")
    }
    
    # --- Monte Carlo simulation ---
    cat("Running Monte Carlo simulation...\n")
    
    #calculate observed_mean = mean of the merged values (the SNPs found in AD that are also found in the trait) = K
    observed_mean_p <- mean(ad_trait$p_value[ad_trait$final_rs_id %in% shared_rs_ids], na.rm = TRUE)
    
    #Simulate 1000 means from K random p-values in trait_clean (N total SNPs)
    simulated_means <- replicate(1000, {
      sampled <- trait_clean[sample(N, K), ]
      mean(sampled$p_value, na.rm = TRUE)
    })
    
    #Calculate empirical p-value (proportion of simulated means <= observed)
    monte_carlo_p <- (sum(simulated_means <= observed_mean_p) + 1) / (length(simulated_means) + 1) #pseudocount correction so that p-value is never 0
    
    cat(sprintf("Observed mean p-value (K = %d): %.4g\n", K, observed_mean_p))
    cat(sprintf("Monte Carlo simulation p-value: %.4g\n", monte_carlo_p))
    
    #Plot histogram
    min_simulated_p <- min(simulated_means, na.rm = TRUE)
    
    hist_plot <- ggplot(data.frame(simulated_means), aes(x = simulated_means)) +
      geom_histogram(color = "black", fill = "#69b3a2", bins = 30) +
      theme_minimal(base_size = 14) +
      theme(
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")
      ) +
      xlab("Simulated Mean P-values") +
      ylab("Frequency") +
      ggtitle(paste("Monte Carlo Simulation of Mean P-values\nfor", trait_name)) +
      geom_vline(aes(xintercept = observed_mean_p), color = "red", linetype = "dashed", linewidth = 1.2) +
      annotate("text", 
               x = observed_mean_p, y = Inf, 
               label = sprintf("Observed: %.2e", observed_mean_p),
               vjust = -0.5, hjust = 0.5, color = "red", fontface = "bold", size = 4.5) +
      annotate("text", 
               x = min_simulated_p, y = Inf,
               label = sprintf("Min Simulated: %.2e", min_simulated_p),
               vjust = -2.0, hjust = 0.5, color = "black", fontface = "italic", size = 4)
    
    plot_file <- paste0("images/montecarlo_", gsub(" ", "_", tolower(trait_name)), ".png")
    ggsave(plot_file, hist_plot, width = 8, height = 6)
    cat("Saved plot to:", plot_file, "\n")
  })
}


##3.Usage
#results <- validate_snps(
#  trait_clean_file = "preprocessed_data/HEART_FAIL_2626_effect_snp_clean.tsv",
#  trait_significant_file = "preprocessed_data/HEART_FAIL_2626_effect_snp_significant.tsv",
#  ad_trait_file = "processed_data/AD1_HEART_FAIL2_Test_filtered_snps.tsv",
#  trait_name = "Heart Failure"
#)

