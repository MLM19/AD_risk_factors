# validating results

#1. Libraries
library(dplyr)
library(data.table)

#2. Functions

##Load SNP data and process
load_snp_data <- function(trait_clean_file, trait_significant_file, ad_trait_file) {
  # Read data
  trait_clean_snps <- fread(trait_clean_file)
  trait_significant_snps <- fread(trait_significant_file)
  ad_trait_snps <- fread(ad_trait_file)
  
  # Extract SNPs
  total_trait_snps <- nrow(trait_clean_snps)
  significant_trait_snps <- trait_significant_snps$rs_id
  num_significant_trait_snps <- length(significant_trait_snps)
  num_non_significant_trait_snps <- total_trait_snps - num_significant_trait_snps
  
  # Extract significant SNPs for AD>trait
  significant_ad_trait_snps <- ad_trait_snps %>%
    filter(p_value < 1e-5) %>%
    pull(final_rs_id)
  
  num_significant_ad_trait_snps <- length(significant_ad_trait_snps)
  total_ad_trait_snps <- nrow(ad_trait_snps)
  
  # Count shared significant SNPs
  num_shared_significant_snps <- sum(significant_trait_snps %in% significant_ad_trait_snps)
  num_non_significant_ad_trait_snps <- total_ad_trait_snps - num_significant_ad_trait_snps
  
  # Return all computed values as a list
  return(list(
    total_trait_snps = total_trait_snps,
    num_significant_trait_snps = num_significant_trait_snps,
    num_non_significant_trait_snps = num_non_significant_trait_snps,
    num_significant_ad_trait_snps = num_significant_ad_trait_snps,
    num_non_significant_ad_trait_snps = num_non_significant_ad_trait_snps,
    num_shared_significant_snps = num_shared_significant_snps,
    total_ad_trait_snps = total_ad_trait_snps,
    p_value = ad_trait_snps$p_value,
    rs_id = ad_trait_snps$rs_id
  ))
}


##Create a contingency table 
create_contingency_table <- function(snp_data) {
  contingency_table <- matrix(
    c(snp_data$num_significant_trait_snps, snp_data$num_shared_significant_snps,  # Significant SNPs row
      snp_data$num_non_significant_trait_snps, snp_data$num_non_significant_ad_trait_snps,  # Non-significant SNPs row
      snp_data$total_trait_snps, snp_data$total_ad_trait_snps),  # Total row
    nrow = 3,
    byrow = TRUE
  )
  
  colnames(contingency_table) <- c("Trait", "AD > Trait")
  rownames(contingency_table) <- c("Significant SNPs", "Non-Significant SNPs", "Total")
  
  return(as.table(contingency_table))
}

#3. Validation functions

##Hypergeometric test
validate_hypergeometric <- function(snp_data) {
  p_value <- phyper(snp_data$num_shared_significant_snps - 1, 
                    snp_data$num_significant_trait_snps, 
                    snp_data$total_trait_snps - snp_data$num_significant_trait_snps, 
                    snp_data$total_ad_trait_snps, 
                    lower.tail = FALSE)  # Upper tail probability
  
  return(p_value)
}

##Fisher exact test
validate_fisher_exact <- function(snp_data) {
  # Constructing the 2x2 table for Fisher's Exact Test
  fisher_table <- matrix(
    c(snp_data$num_shared_significant_snps, snp_data$num_non_significant_ad_trait_snps,
      snp_data$num_non_significant_trait_snps, snp_data$num_significant_ad_trait_snps),
    nrow = 2,
    byrow = TRUE
  )
  
  # Fisher's Exact Test without constant addition
  fisher_test <- fisher.test(fisher_table)
  
  # Odds ratio and p-value
  odds_ratio <- fisher_test$estimate
  p_value <- fisher_test$p.value
  
  return(list(odds_ratio = odds_ratio, p_value = p_value))
}

##Montecarlo (p-value comparison)
validate_monte_carlo <- function(snp_data, num_iterations = 1000) {
  # Extract p-values of significant AD>Trait SNPs
  significant_ad_trait_snps <- snp_data$significant_ad_trait_snps
  significant_p_values <- snp_data$p_value[snp_data$rs_id %in% significant_ad_trait_snps]
  
  # Check if the significant_p_values is numeric
  if (!is.numeric(significant_p_values)) {
    stop("Significant p-values are not numeric.")
  }
  
  mean_significant_p_value <- mean(significant_p_values, na.rm = TRUE)
  
  # Get the total number of SNPs in the COVID dataset
  total_snp_data <- snp_data$total_trait_snps
  
  # Initialize variable to store the count of random samples with lower p-value mean
  count_lower_p_value <- 0
  
  # Monte Carlo simulation: Repeat random sampling of SNPs
  for (i in 1:num_iterations) {
    # Randomly sample 4 SNPs from the total SNP dataset
    random_sample <- sample(1:total_snp_data, size = length(significant_ad_trait_snps))
    
    # Calculate the mean p-value for the random sample
    random_sample_p_values <- snp_data$p_value[random_sample]
    
    # Check if random_sample_p_values is numeric
    if (!is.numeric(random_sample_p_values)) {
      next  # Skip this iteration if p-values are not numeric
    }
    
    mean_random_p_value <- mean(random_sample_p_values, na.rm = TRUE)
    
    # Compare the mean of the random sample with the mean of the significant SNPs
    if (!is.na(mean_random_p_value) && mean_random_p_value < mean_significant_p_value) {
      count_lower_p_value <- count_lower_p_value + 1
    }
  }
  
  # Calculate p-value as proportion of simulations where random sample p-value was smaller
  monte_carlo_p_value <- count_lower_p_value / num_iterations
  
  # Return p-value and other relevant results
  return(list(monte_carlo_p_value = monte_carlo_p_value, mean_significant_p_value = mean_significant_p_value))
}

# Usage
trait_clean_file <- "preprocessed_data/COVID_0403_effect_snp_clean.tsv"
trait_significant_file <- "preprocessed_data/COVID_0403_effect_snp_significant.tsv"
ad_trait_file <- "processed_data/AD1_COVID2_Test_filtered_snps.tsv"

# Load data
snp_data <- load_snp_data(trait_clean_file, trait_significant_file, ad_trait_file)

# Create contingency table
contingency_table <- create_contingency_table(snp_data)
print(contingency_table)

# Validate using hypergeometric test
p_value <- validate_hypergeometric(snp_data)
print(paste("Hypergeometric p-value:", p_value))

# Validate using Fisher's Exact Test
fisher_result <- validate_fisher_exact(snp_data)
print(paste("Fisher's Exact Test p-value:", fisher_result$p_value))
print(paste("Fisher's Exact Test Odds Ratio:", fisher_result$odds_ratio))

#Validate using Monte carlo result
monte_carlo_result <- validate_monte_carlo(snp_data)
print(paste("Monte Carlo p-value:", monte_carlo_result$monte_carlo_p_value))
print(paste("Mean p-value of significant SNPs:", monte_carlo_result$mean_significant_p_value))
