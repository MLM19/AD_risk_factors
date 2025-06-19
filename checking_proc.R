# --- Libraries ---
library(data.table)
library(dplyr)
library(ggplot2)

# --- Parameters ---
SIGNIFICANCE_THRESHOLD <- 1e-6   # Modifiable threshold
LD_WINDOW <- 200000              # LD pruning window (200kb)

trait_codes <- c(
  "ALCOH_3727", "APNEA_7362", "ASBEST_3282", "BRE_CANCER_174", "COR_ARTERY_2314", "COVID_0403",
  "DENTAL_521.1", "DEPRES_6475", "DIAB_GEST_6553", "DIAB_TYP1_3791", "DIAB_TYP2_8926",
  "DIET_9364", "EBS_20002", "EDU_1875", "FRECKL_6091", "HAIR_COL_1747.4", "HEART_FAIL_2626",
  "HerpZo_8721", "HHV_6906", "HIV_5531", "HSV_5524", "HYPERTEN_9346", "INFLU_2107",
  "INSOM_7286", "IRON_ANEM_280", "LACT_INT_4159", "LE_BO_DEM_1390", "LONE_1704", "MENOP_2553",
  "OBES_5785", "PARKINSON_9325", "PIÑA_104570", "PSY_STRESS_4320", "RHINO_4043", "SKIN_COL_9962",
  "SMOK_3686", "STRESS_5871", "STROK_4723", "VAS_DEM_290.16", "VITB_5741", "VITB12_5799",
  "VITD_5742"
)

# --- Load AD Data (phenotype 2) ---
cat("Loading AD data...\n")
ad_df <- fread("preprocessed_data/AD_7158_effect_snp_clean.tsv", select = c("rs_id", "chr", "pos", "p_value"))
ad_df <- na.omit(ad_df)
ad_df[, rs_id := as.character(rs_id)]

# --- LD pruning function ---
ld_prune <- function(df, LD_WINDOW = 200000) {
  #Order by significance (lowest p-value first), then by chromosome and position
  df <- df[order(df$p_value, df$chr, df$pos), ]
  
  #Initialize with the most significant SNP
  pruned <- df[1, ]
  
  for (i in 2:nrow(df)) {
    #Compare against the last SNP added to the pruned list
    last_chr <- pruned$chr[nrow(pruned)]
    last_pos <- pruned$pos[nrow(pruned)]
    
    #Only add SNP if it's on a different chromosome or far enough away
    if (df$chr[i] != last_chr || abs(df$pos[i] - last_pos) > LD_WINDOW) {
      pruned <- rbind(pruned, df[i, ])
    }
  }
  
  return(pruned)
}

# --- Trait label cleaner ---
clean_labels <- function(x) {
  x %>% sub("_\\d+(\\.\\d+)?$", "", .)
}

# --- Initialize results ---
results <- list()

# --- Main loop ---
for (trait in trait_codes) {
  cat("Processing trait:", trait, "\n")
  file_path <- paste0("preprocessed_data/", trait, "_effect_snp_clean.tsv")
  if (!file.exists(file_path)) {
    cat("  File not found.\n")
    results[[trait]] <- list(trait = trait, trimean = 0, n_snps = 0)
    next
  }
  
  # Step 1: Load trait (phenotype 1) data
  trait_df <- fread(file_path, select = c("rs_id", "chr", "pos", "p_value"))
  trait_df <- na.omit(trait_df)
  trait_df[, rs_id := as.character(rs_id)]
  
  # Step 2: Select significant variants
  sig_trait <- trait_df[p_value <= SIGNIFICANCE_THRESHOLD]
  if (nrow(sig_trait) == 0) {
    results[[trait]] <- list(trait = trait, trimean = 0, n_snps = 0)
    next
  }
  
  # Step 3: LD prune the significant SNPs
  if (nrow(sig_trait) <= 1){
    sig_trait_pruned <- sig_trait
  } else {
    sig_trait_pruned <- ld_prune(sig_trait)
  }

  # Step 4: Merge pruned phenotype 1 SNPs with full AD set (all SNPs)
  ##ensure data types to match
  sig_trait_pruned <- sig_trait_pruned %>%
    mutate(
      chr = as.factor(chr),
      pos = as.numeric(pos)
    )
  
  ad_df <- ad_df %>%
    mutate(
      chr = as.factor(chr),
      pos = as.numeric(pos),
      p_value = as.numeric(p_value)
    )
  
  ##merge by chr and pos
  merged <- sig_trait_pruned %>%
    inner_join(ad_df, by = c("chr", "pos"), suffix = c("_trait", "_ad"))
  
  ##keep the rs_id form either of the two files
  merged <- merged %>%
    mutate(
      final_rs_id = coalesce(rs_id_trait, rs_id_ad)
    )
  
  
  if (nrow(merged) == 0) {
    results[[trait]] <- list(trait = trait, trimean = 0, n_snps = 0)
    next
  }
  
  # Step 5: Tukey trimean of -log10(p_value) for AD (phenotype 2)
  logp <- -log10(merged$p_value_ad)
  q1 <- quantile(logp, 0.25)
  q2 <- quantile(logp, 0.50)
  q3 <- quantile(logp, 0.75)
  trimean <- (q1 + 2*q2 + q3) / 4
  
  # Optional: Save merged SNPs
  out_path <- paste0("processed_data/", trait, "_AD_overlap_", SIGNIFICANCE_THRESHOLD, ".tsv")
  #fwrite(merged, out_path, sep = "\t")
  
  # Step 6: Save results
  results[[trait]] <- list(
    trait = trait,
    trimean = trimean,
    n_snps = nrow(merged)
  )
  
  #print debug
  cat("  Significant SNPs:", nrow(sig_trait), "\n")
  cat("  After LD pruning:", nrow(sig_trait_pruned), "\n")
  cat("  After merge with AD:", nrow(merged), "\n")
  
  # Step 7: Threshold validation via Monte Carlo simulation
  cat("  Running trimean null distribution simulation...\n")
  set.seed(42)  # for reproducibility
  
  n_overlap <- nrow(merged)
  logp_ad_all <- -log10(ad_df$p_value)
  
  # Only sample from AD clean SNPs that have valid p-values
  ad_valid <- ad_df[!is.na(p_value)]
  n_ad_total <- nrow(ad_valid)
  
  simulated_trimeans <- replicate(1000, {
    sampled <- ad_valid[sample(n_ad_total, n_overlap), ]
    logp <- -log10(sampled$p_value)
    q1 <- quantile(logp, 0.25)
    q2 <- quantile(logp, 0.50)
    q3 <- quantile(logp, 0.75)
    (q1 + 2 * q2 + q3) / 4
  })
  
  # Step 8: Plot distribution of simulated trimeans
  hist_plot <- ggplot(data.frame(simulated_trimeans), aes(x = simulated_trimeans)) +
    geom_histogram(color = "black", fill = "#69b3a2", bins = 30) +
    geom_vline(aes(xintercept = trimean), color = "red", linetype = "dashed", linewidth = 1.2) +
    annotate("text",
             x = trimean, y = Inf,
             label = sprintf("Observed trimean: %.2f", trimean),
             vjust = -0.5, hjust = 0.5, color = "red", fontface = "bold", size = 4.2) +
    theme_bw(base_size = 14) +
    xlab("Simulated Tukey Trimeans") +
    ylab("Frequency") +
    ggtitle(paste("Null distribution of Tukey trimeans\n", clean_labels(trait), " vs AD"))
  
  # Save the plot
  plot_path <- file.path("Images/", paste0(trait, "_AD_trimean_distribution_", SIGNIFICANCE_THRESHOLD, ".png"))
  ggsave(plot_path, hist_plot, width = 7, height = 5.5)
  
  cat("  Trimean simulation plot saved to", plot_path, "\n")
}

# --- Prepare and plot results ---
result_df <- rbindlist(results, fill = TRUE)
result_df$trait <- clean_labels(result_df$trait)
result_df <- result_df %>% arrange(desc(trimean))
result_df$trait <- factor(result_df$trait, levels = rev(result_df$trait))

# --- Plot ---
plt <- ggplot(result_df, aes(x = trait, y = trimean)) +
  geom_bar(stat = "identity", fill = "darkred") +
  # Place SNP counts to the RIGHT of each bar
  geom_text(aes(label = n_snps), 
            hjust = -0.1, 
            size = 3.5) +
  coord_cartesian(clip = "off") +
  coord_flip() +
  ggtitle(paste0("Genetic association of traits with AD (p < ", SIGNIFICANCE_THRESHOLD, ")")) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.margin = margin(5.5, 50, 5.5, 5.5)  # Add space to the RIGHT for labels
  ) +
  labs(
    x = "Trait",
    y = "Tukey Trimean of -log10(p) in AD"
  ) +
  ylim(0, max(result_df$trimean, na.rm = TRUE) + 1)  # Add space on the right for SNP counts

# --- Save plot ---
#ggsave(file.path("processed_data", paste0("genetic_association_AD_", SIGNIFICANCE_THRESHOLD, ".png")),
#       plt, width = 7, height = 6)


