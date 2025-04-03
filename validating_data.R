## Version with pseudocount adjustments to prevent zero division,
## continuity correction on contingency table,
## and Regularized Fisher test with Haldane correction

## 1. Load necessary libraries
library(dplyr)  # Load dplyr for data manipulation functions (e.g., inner_join)

## 2. Functions

## Function to compute N, K, k, and n from files with pseudocount adjustments
compute_snp_counts <- function(trait_clean_file, trait_significant_file, ad_trait_file) {
  # Print a message indicating that data files are being loaded
  cat("Loading data files...\n")
  
  # Load the preprocessed clean trait data (all SNPs available)
  trait_clean <- read.delim(trait_clean_file, header = TRUE)
  # Load the trait significant data (SNPs significant for the environmental trait)
  trait_significant <- read.delim(trait_significant_file, header = TRUE)
  # Load the AD trait data (merged SNPs with AD, may include significant hits)
  ad_trait <- read.delim(ad_trait_file, header = TRUE)
  
  # Debugging: print the column names of the datasets to check the available join columns
  cat("Debugging...\n")
  print(colnames(trait_significant))  # Expecting a column "rs_id" in trait_significant
  print(colnames(ad_trait))           # Expecting a column "final_rs_id" in ad_trait
  
  # Check if the expected columns exist, stop execution if not found
  if (!"rs_id" %in% colnames(trait_significant)) {
    stop("Error: 'rs_id' not found in trait_significant dataset")
  }
  if (!"final_rs_id" %in% colnames(ad_trait)) {
    stop("Error: 'final_rs_id' not found in ad_trait dataset")
  }
  
  cat("Computing SNP counts...\n")
  
  # Compute total SNPs (N) from the trait_clean file with a small pseudocount added
  N <- nrow(trait_clean) + 1e-6
  # Compute significant SNPs (n) from the trait_significant file; ensure at least 1 to avoid zero
  n <- max(nrow(trait_significant), 1)
  # Compute significant SNPs in AD (K) using a threshold (p < 1e-5) with a minimum of 1 for stability
  K <- max(sum(ad_trait$p_value < 1e-5), 1)
  # Compute total overlapping SNPs (k) from the ad_trait file (representing merged overlapping SNPs); ensure at least 1
  k <- max(nrow(ad_trait), 1)
  
  # Print the computed values for debugging and verification
  cat("Computed values:\n")
  cat("N (Total SNPs):", N, "\n")
  cat("K (Significant SNPs in AD):", K, "\n")
  cat("k (Overlapping SNPs):", k, "\n")
  cat("n (Significant SNPs in Trait):", n, "\n")
  
  # Return the computed values as a list
  return(list(N = N, K = K, k = k, n = n))
}

## Function to create contingency table and perform validation tests
validate_snps <- function(trait_clean_file, trait_significant_file, ad_trait_file) {
  cat("Starting validation process...\n")
  
  # Compute SNP counts by calling the compute_snp_counts() function
  snp_counts <- compute_snp_counts(trait_clean_file, trait_significant_file, ad_trait_file)
  N <- snp_counts$N  # Total SNPs in the dataset
  K <- snp_counts$K  # Significant SNPs in AD (p < 1e-5)
  k <- snp_counts$k  # Total overlapping SNPs (from merged AD and trait)
  n <- snp_counts$n  # Significant SNPs in the environmental trait dataset
  
  cat("Creating contingency table...\n")
  
  # Construct the contingency table:
  # For Environmental Trait: 
  #       Statistically Significant (SS) = n (from trait_significant)
  #       Not Statistically Significant (NSS) = N - n (remaining SNPs in trait_clean)
  # For Overlap AD & Trait:
  #       Statistically Significant (SS) = K (from AD significant SNPs)
  #       Not Statistically Significant (NSS) = k - K (remaining overlapping SNPs)
  contingency_table <- matrix(
    c(n, K,         # First row: SS values (Environmental Trait = n, Overlap AD & Trait = K)
      N - n, k - K),  # Second row: NSS values (Environmental Trait = N - n, Overlap AD & Trait = k - K)
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("SS", "NSS"),  # Row names: Statistically Significant and Not Statistically Significant
      c("Environmental Trait", "Overlap AD & Trait")  # Column names: Categories for trait and overlapping data
    )
  )
  
  # Print the contingency table for verification
  cat("Contingency Table:\n")
  print(contingency_table)
  
  cat("Performing Hypergeometric Test...\n")
  # Hypergeometric Test:
  # This tests the probability of observing at least K significant SNPs in n draws,
  # given k overlapping SNPs and N total SNPs in the dataset.
  hyper_pval <- phyper(K - 1, k, N - k, n, lower.tail = FALSE)
  cat("Hypergeometric Test p-value:", hyper_pval, "\n")
  
  cat("Performing Regularized Fisher's Exact Test (Haldane Correction)...\n")
  # Fisher's Exact Test with Haldane correction:
  # Add a small constant (0.25) to each cell to regularize (avoid zeros) for small counts.
  fisher_result <- fisher.test(contingency_table + 0.25)
  cat("Fisher's Exact Test Results:\n")
  print(fisher_result)
  
  cat("Performing Monte Carlo Simulation...\n")
  # Monte Carlo Simulation:
  # Here, example p-values are generated for both overlapping and environmental trait sets.
  # In practice, you would use the actual p-values from your data.
  p_values_overlap <- runif(k, 0, 0.05)  # Generate random p-values for overlapping SNPs
  p_values_trait <- runif(n, 0, 0.05)     # Generate random p-values for environmental trait SNPs
  observed_mean_p <- mean(p_values_overlap)  # Compute the observed mean of the overlapping p-values
  # Replicate the random sampling process 1000 times to build a null distribution of mean p-values
  simulated_means <- replicate(1000, mean(sample(p_values_trait, k, replace = TRUE)))
  # Calculate the proportion of simulated means that are less than or equal to the observed mean
  monte_carlo_pval <- mean(simulated_means <= observed_mean_p)
  cat("Monte Carlo Simulation p-value:", monte_carlo_pval, "\n")
  
  # Return a list containing the contingency table and p-values from each test
  return(list(
    contingency_table = contingency_table,
    hyper_pval = hyper_pval,
    fisher_pval = fisher_result$p.value,
    mc_pval = monte_carlo_pval
  ))
}

## 3. Example usage

# Define file paths for the data files
trait_clean_file <- "preprocessed_data/HerpZo_8721_effect_snp_clean.tsv"         # Full SNP dataset for the trait
trait_significant_file <- "preprocessed_data/HerpZo_8721_effect_snp_significant.tsv"  # Significant SNPs for the trait
ad_trait_file <- "processed_data/AD1_HerpZo2_Test_filtered_snps.tsv"               # Merged overlapping SNPs (AD and trait)

# Run validation tests by calling the validate_snps function with the file paths
cat("Running validation tests...\n")
results <- validate_snps(trait_clean_file, trait_significant_file, ad_trait_file)
cat("Validation process complete.\n")
