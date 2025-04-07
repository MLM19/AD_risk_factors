##1. Load required libraries
library(dplyr)
library(ggplot2)

##2. Functions
###Function to compute N, K, k, and n with pseudocount adjustments
compute_snp_counts <- function(trait_clean_file, trait_significant_file, ad_trait_file) {
  cat("Loading data files...\n")
  
  trait_clean <- read.delim(trait_clean_file, header = TRUE)
  trait_significant <- read.delim(trait_significant_file, header = TRUE)
  ad_trait <- read.delim(ad_trait_file, header = TRUE)
  
  if (!"rs_id" %in% colnames(trait_significant)) stop("Missing 'rs_id' in trait_significant")
  if (!"final_rs_id" %in% colnames(ad_trait)) stop("Missing 'final_rs_id' in ad_trait")
  
  cat("Computing SNP counts...\n")
  
  N <- nrow(trait_clean) + 1e-6
  n <- max(nrow(trait_significant), 1)
  K <- max(sum(ad_trait$p_value < 1e-5), 1)
  k <- max(nrow(ad_trait), 1)
  
  return(list(
    N = N, K = K, k = k, n = n,
    trait_clean = trait_clean,
    trait_significant = trait_significant,
    ad_trait = ad_trait
  ))
}


###Main validation function with real p-values + plots
validate_snps <- function(trait_clean_file, trait_significant_file, ad_trait_file, trait_name = "Trait") {
  cat("Starting validation process for:", trait_name, "\n")
  
  counts <- compute_snp_counts(trait_clean_file, trait_significant_file, ad_trait_file)
  with(counts, {
    shared_snps <- intersect(trait_significant$rs_id, ad_trait$final_rs_id)
    SS <- length(shared_snps)
    
    cat("Significant Shared SNPs (SS):", SS, "\n")
    
    # --- Hypergeometric test (dhyper) ---
    hyper_p <- dhyper(SS, K, N - K, n)
    cat(sprintf("Hypergeometric test p-value: %.4g\n", hyper_p))
    if (SS == 0) {
      cat("Interpretation: No significant SNPs shared — possible weak or no association.\n\n")
    } else if (hyper_p < 0.05) {
      cat("Interpretation: Enrichment in shared SNPs is significant (unlikely due to chance).\n\n")
    } else {
      cat("Interpretation: No statistically significant enrichment of shared SNPs.\n\n")
    }
    
    # --- Fisher's Exact Test ---
    not_shared_trait <- n - SS
    not_shared_ad <- K - SS
    rest <- N - K - n + SS
    contingency_table <- matrix(c(SS, not_shared_trait, not_shared_ad, rest), nrow = 2)
    fisher_result <- fisher.test(contingency_table + 0.25)
    cat(sprintf("Fisher's exact test p-value: %.4g\n", fisher_result$p.value))
    cat(sprintf("Fisher's exact test odds ratio: %.3f\n", fisher_result$estimate))
    
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
    simulated_means <- replicate(1000, {
      random_rs_ids <- sample(trait_clean$rs_id, n)
      overlap <- trait_clean[trait_clean$rs_id %in% random_rs_ids & 
                               trait_clean$rs_id %in% ad_trait$final_rs_id, ]
      if (nrow(overlap) == 0) return(0)
      mean(overlap$p_value, na.rm = TRUE)
    })
    
    observed_mean_p <- if (SS == 0) 0 else {
      mean(ad_trait$p_value[ad_trait$final_rs_id %in% shared_snps], na.rm = TRUE)
    }
    
    cat(sprintf("Observed mean p-value: %.4g\n", observed_mean_p))
    
    monte_carlo_p <- sum(simulated_means <= observed_mean_p) / length(simulated_means)
    cat(sprintf("Monte Carlo simulation p-value: %.4g\n", monte_carlo_p))
    if (monte_carlo_p < 0.05) {
      cat("Interpretation: Observed mean p-value is lower than expected by chance.\n\n")
    } else {
      cat("Interpretation: Observed mean p-value could occur by chance.\n\n")
    }
    
    # --- Plot histogram ---
    library(ggplot2)
    hist_plot <- ggplot(data.frame(simulated_means), aes(x = simulated_means)) +
      geom_histogram(color = "black", fill = "#69b3a2", bins = 30) +
      theme_minimal(base_size = 14) +
      theme(
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")
      ) +
      xlab("Simulated Mean P-values") +
      ylab("Frequency") +
      ggtitle(paste("Monte Carlo Simulation of Mean P-values\nfor", trait_name))
    
    if (SS > 0) {
      hist_plot <- hist_plot + 
        geom_vline(xintercept = observed_mean_p, color = "red", linetype = "dashed", size = 1.2)
    }
    
    plot_file <- paste0("images/montecarlo_", gsub(" ", "_", tolower(trait_name)), ".png")
    ggsave(plot_file, hist_plot, width = 8, height = 6)
    cat("Saved plot to:", plot_file, "\n")
    
    return(list(
      SS = SS,
      hyper_p = hyper_p,
      fisher_p = fisher_result$p.value,
      monte_carlo_p = monte_carlo_p,
      observed_mean_p = observed_mean_p
    ))
  })
}

##3.Usage
results <- validate_snps(
  trait_clean_file = "preprocessed_data/PARKINSON_9325_effect_snp_clean.tsv",
  trait_significant_file = "preprocessed_data/PARKINSON_9325_effect_snp_significant.tsv",
  ad_trait_file = "processed_data/AD1_PARKINSON2_Test_filtered_snps.tsv",
  trait_name = "Parkinson"
)

