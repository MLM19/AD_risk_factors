# filter_and_plot_significant_snps.R

library(data.table)
library(dplyr)
library(qqman)

skip_traits <- c()

# Set threshold
p_value_threshold <- 1e-6

# List all *_effect_snp_clean.tsv files
effect_files <- list.files("preprocessed_data", pattern = "_effect_snp_clean.tsv$", full.names = TRUE)

# Output folder
dir.create("Images", showWarnings = FALSE, recursive = TRUE)

# Loop through all effect files
for (effect_path in effect_files) {
  # Extract the prefix (trait code)
  db_name <- sub("_effect_snp_clean.tsv$", "", basename(effect_path))
  if (db_name %in% skip_traits) {
    message("Skipping trait (already processed): ", db_name)
    next
  }
  message("Processing trait: ", db_name)
  
  # Define paths
  snp_path <- file.path("preprocessed_data", paste0(db_name, "_snp_clean.tsv"))
  sig_effect_path <- file.path("preprocessed_data", paste0(db_name, "_effect_snp_significant.tsv"))
  sig_snp_path <- file.path("preprocessed_data", paste0(db_name, "_snp_significant.tsv"))
  manhattan_path <- file.path("Images", paste0(db_name, "_manhattan.png"))
  
  # Load the data
  effect_df <- fread(effect_path)
  snp_df <- fread(snp_path)
  
  # Filter significant SNPs
  sig_effect_df <- effect_df %>% filter(p_value <= p_value_threshold)
  
  # Join with SNP table using rs_id
  sig_snp_df <- snp_df %>% filter(rs_id %in% sig_effect_df$rs_id)
  
  # Save filtered results
  fwrite(sig_effect_df, sig_effect_path, sep = "\t")
  fwrite(sig_snp_df, sig_snp_path, sep = "\t")
  
  # Manhattan plot (only from effect_df, showing all for context)
  effect_df <- effect_df %>%
    mutate(
      chr = case_when(
        chr == "X" ~ 23,
        chr == "Y" ~ 24,
        TRUE ~ suppressWarnings(as.numeric(chr))
      ),
      pos = as.numeric(pos),
      p_value = as.numeric(p_value)
    ) %>%
    filter(
      !is.na(chr) & chr %in% 1:24,
      !is.na(pos) & is.finite(pos),
      !is.na(p_value) & is.finite(p_value),
      p_value > 0 
    )
  
  if (nrow(effect_df) == 0 || all(is.na(effect_df$p_value)) || length(unique(effect_df$p_value)) == 1) {
    message("Skipping plot: p_value issue in trait: ", db_name)
    next
  }
  
  png(manhattan_path, width = 1200, height = 800, res = 150)
  manhattan(
    effect_df,
    chr = "chr",
    bp = "pos",
    p = "p_value",
    snp = "rs_id",
    logp = TRUE,
    main = paste(sub("_\\d+$", "", db_name), "Manhattan Plot"),
    suggestiveline = -log10(p_value_threshold),
    genomewideline = -log10(5e-8),
    col = c("black", "grey")
  )
  dev.off()
  
  message("Finished processing: ", db_name)
}

#source("filter_and_plot_significant_snps.R")