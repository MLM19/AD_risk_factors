# Validating results

# 1. Libraries
library(dplyr)
library(data.table)

# 2. Functions

## Load SNP data and process
load_snp_data <- function(trait_clean_file, trait_significant_file, ad_trait_file) {
  # Read data
  trait_clean_snps <- fread(trait_clean_file)
  trait_significant_snps <- fread(trait_significant_file)
  ad_trait_snps <- fread(ad_trait_file)
  
  # Debug: Print dimensions of input data
  print(paste("Trait clean SNPs:", nrow(trait_clean_snps)))
  print(paste("Trait significant SNPs:", nrow(trait_significant_snps)))
  print(paste("AD > Trait SNPs:", nrow(ad_trait_snps)))
  
  # Extract SNP counts
  total_trait_snps <- nrow(trait_clean_snps)
  num_significant_trait_snps <- nrow(trait_significant_snps)
  num_non_significant_trait_snps <- total_trait_snps - num_significant_trait_snps
  
  # Extract significant SNPs for AD > Trait
  num_significant_ad_trait_snps <- nrow(ad_trait_snps %>% filter(p_value < 1e-5))
  num_non_significant_ad_trait_snps <- nrow(ad_trait_snps) - num_significant_ad_trait_snps
  
  return(list(
    total_trait_snps = total_trait_snps,
    num_significant_trait_snps = num_significant_trait_snps,
    num_non_significant_trait_snps = num_non_significant_trait_snps,
    num_significant_ad_trait_snps = num_significant_ad_trait_snps,
    num_non_significant_ad_trait_snps = num_non_significant_ad_trait_snps,
    p_values = ad_trait_snps$p_value
  ))
}


# Create contingency table
create_contingency_table <- function(snp_data) {
  # Construct the contingency table
  contingency_table <- matrix(
    c(
      snp_data$num_significant_trait_snps, snp_data$num_significant_ad_trait_snps, # Significant SNPs row
      snp_data$num_non_significant_trait_snps, snp_data$num_non_significant_ad_trait_snps # Non-Significant SNPs row
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  colnames(contingency_table) <- c("Trait", "AD > Trait")
  rownames(contingency_table) <- c("Significant SNPs", "Non-Significant SNPs")
  
  return(as.table(contingency_table))
}


# Hypergeometric test
validate_hypergeometric <- function(snp_data) {
  # Handle case where there are no significant SNPs in AD > Trait
  if (snp_data$num_significant_ad_trait_snps == 0) {
    return(1)  # p-value of 1 indicates no enrichment
  }
  
  # Debug: Print inputs to phyper
  #print(paste("num_significant_trait_snps:", snp_data$num_significant_trait_snps))
  #print(paste("total_trait_snps:", snp_data$total_trait_snps))
  #print(paste("num_significant_ad_trait_snps:", snp_data$num_significant_ad_trait_snps))
  
  # Calculate complement (non-significant SNPs in Trait dataset)
  complement <- snp_data$total_trait_snps - snp_data$num_significant_trait_snps
  print(paste("Complement (total_trait_snps - num_significant_trait_snps):", complement))
  
  # Use minimum overlap as a proxy for shared SNPs
  num_shared_significant_snps <- min(snp_data$num_significant_trait_snps, snp_data$num_significant_ad_trait_snps)
  
  # Debug: Print calculated shared SNP count
  print(paste("num_shared_significant_snps:", num_shared_significant_snps))
  
  # Perform hypergeometric test
  p_value <- phyper(
    num_shared_significant_snps - 1, 
    snp_data$num_significant_trait_snps, 
    complement, 
    (snp_data$num_significant_ad_trait_snps + snp_data$num_non_significant_ad_trait_snps), 
    lower.tail = FALSE
  )
  
  return(p_value)
}

# Fisher exact test validation
validate_fisher_exact <- function(snp_data) {
  # Handle case where there are no significant SNPs in AD > Trait
  if (snp_data$num_significant_ad_trait_snps == 0) {
    return(list(
      odds_ratio = NA,
      p_value = 1
    ))
  }
  
  # Create a 2x2 contingency table
  fisher_table <- matrix(
    c(
      snp_data$num_significant_trait_snps, snp_data$num_significant_ad_trait_snps,
      snp_data$num_non_significant_trait_snps, snp_data$num_non_significant_ad_trait_snps
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  # Perform Fisher's Exact Test
  fisher_test <- fisher.test(fisher_table)
  
  # Extract odds ratio and p-value
  odds_ratio <- fisher_test$estimate
  p_value <- fisher_test$p.value
  
  return(list(odds_ratio = odds_ratio, p_value = p_value))
}

# Monte Carlo simulation
validate_monte_carlo <- function(snp_data, num_iterations = 1000) {
  # Handle case where there are no significant SNPs in AD > Trait
  if (snp_data$num_significant_ad_trait_snps == 0) {
    return(list(
      monte_carlo_message = "No significant SNPs detected in AD > Trait.",
      monte_carlo_p_value = NA,
      observed_mean_p_value = NA
    ))
  }
  
  # Step 1: Calculate the observed mean p-value for significant SNPs
  significant_p_values <- snp_data$p_values[1:snp_data$num_significant_ad_trait_snps]
  
  # Ensure p-values are numeric
  if (!is.numeric(significant_p_values)) {
    stop("Error: Significant p-values are not numeric.")
  }
  
  # Check for missing values in significant p-values
  if (any(is.na(significant_p_values))) {
    stop("Error: Significant p-values contain NA values.")
  }
  
  # Calculate the observed mean p-value
  observed_mean_p_value <- mean(significant_p_values, na.rm = TRUE)
  
  # Step 2: Initialize counter for random samples with lower mean
  count_lower <- 0
  
  # Step 3: Perform Monte Carlo simulation
  for (i in seq_len(num_iterations)) {
    # Randomly sample p-values from the total dataset
    random_sample <- sample(snp_data$p_values, size = snp_data$num_significant_ad_trait_snps, replace = FALSE)
    
    # Check for NA values in the random sample
    if (any(is.na(random_sample))) {
      next # Skip this iteration if there are NA values
    }
    
    # Calculate the mean of the random sample
    random_mean <- mean(random_sample, na.rm = TRUE)
    
    # Compare means and increment counter if random sample mean is smaller
    if (!is.na(random_mean) && random_mean < observed_mean_p_value) {
      count_lower <- count_lower + 1
    }
  }
  
  # Step 4: Calculate Monte Carlo p-value
  monte_carlo_p_value <- count_lower / num_iterations
  
  # If no random samples have a lower mean, report an upper bound instead of exact zero
  #reporting 0 implies total certainty
  if (monte_carlo_p_value == 0) {
    monte_carlo_message <- paste("<", format(1 / num_iterations, scientific = TRUE))
    monte_carlo_p_value <- paste("<", format(1 / num_iterations, scientific = TRUE))
  } else {
    monte_carlo_message <- format(monte_carlo_p_value, scientific = TRUE)
  }
  
  return(list(
    monte_carlo_p_value = monte_carlo_p_value,
    observed_mean_p_value = observed_mean_p_value
  ))
}


##3. Usage example
trait_clean_file <- "preprocessed_data/HerpZo_8721_effect_snp_clean.tsv"
trait_significant_file <- "preprocessed_data/HerpZo_8721_effect_snp_significant.tsv"
ad_trait_file <- "processed_data/AD1_HerpZo2_Test_filtered_snps.tsv"

# Load data
snp_data <- load_snp_data(trait_clean_file, trait_significant_file, ad_trait_file)

# Create contingency table
contingency_table <- create_contingency_table(snp_data)
print(contingency_table)

# Validate using hypergeometric test
hypergeo_p_value <- validate_hypergeometric(snp_data)
print(paste("Hypergeometric p-value:", hypergeo_p_value))

# validate using fisher exact test
fisher_result <- validate_fisher_exact(snp_data)
print(paste("Fisher's Exact Test p-value:", fisher_result$p_value))
print(paste("Fisher's Exact Test Odds Ratio:", fisher_result$odds_ratio))

# Validate using Monte Carlo simulation
monte_carlo_result <- validate_monte_carlo(snp_data, num_iterations = 1000000)
# Print results
print(paste("Monte Carlo p-value:", monte_carlo_result$monte_carlo_p_value))
print(paste("Observed Mean p-value of Significant SNPs:", monte_carlo_result$observed_mean_p_value))
