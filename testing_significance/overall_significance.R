# --- Libraries ---
library(data.table)
library(ggplot2)
library(dplyr)

# --- Parameters ---
SIGNIFICANCE_THRESHOLD <- 1e-8  # modifiable

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

# --- Load AD summary statistics ---
cat("DEBUG: Loading AD summary statistics...\n")
ad_df <- fread("preprocessed_data/AD_7158_effect_snp_clean.tsv", select = c("rs_id", "chr", "pos", "p_value"))
cat("DEBUG: AD SNPs loaded:", nrow(ad_df), "\n")

ad_df <- na.omit(ad_df)
cat("DEBUG: AD SNPs after NA removal:", nrow(ad_df), "\n")

ad_df[, `:=`(
  chr = as.character(chr),
  pos = as.integer(pos),
  rs_id = as.character(rs_id),
  p_value = as.numeric(p_value)
)]

ad_df <- ad_df %>% rename(rs_id_ad = rs_id)

results <- list()

# --- Loop through traits ---
for (trait in trait_codes) {
  file_path <- paste0("preprocessed_data/", trait, "_effect_snp_clean.tsv")
  cat("DEBUG: Processing trait:", trait, "\n")
  
  if (!file.exists(file_path)) {
    cat("WARNING: File not found:", file_path, "\n")
    next
  }
  
  # Load trait data
  trait_df <- fread(file_path, select = c("rs_id", "chr", "pos", "p_value"))
  cat("DEBUG: Loaded", nrow(trait_df), "rows for", trait, "\n")
  
  trait_df <- na.omit(trait_df)
  cat("DEBUG: After NA removal:", nrow(trait_df), "rows remain for", trait, "\n")
  
  trait_df[, `:=`(
    chr = as.character(chr),
    pos = as.integer(pos),
    rs_id = as.character(rs_id),
    p_value = as.numeric(p_value)
  )]
  
  trait_df <- trait_df %>% rename(rs_id_trait = rs_id)
  
  # Filter significant SNPs
  sig_trait_df <- trait_df[p_value < SIGNIFICANCE_THRESHOLD]
  cat("DEBUG: - Significant SNPs (p <", SIGNIFICANCE_THRESHOLD, "):", nrow(sig_trait_df), "\n")
  
  # Merge with AD data
  merged <- sig_trait_df %>%
    inner_join(ad_df, by = c("chr", "pos"))
  
  cat("DEBUG: - SNPs shared with AD:", nrow(merged), "\n")
  
  if (nrow(merged) > 0) {
    merged <- merged %>% mutate(final_rs_id = coalesce(rs_id_trait, rs_id_ad))
    
    log_p <- -log10(merged$p_value.y)
    q1 <- quantile(log_p, 0.25, na.rm = TRUE)
    q2 <- quantile(log_p, 0.50, na.rm = TRUE)
    q3 <- quantile(log_p, 0.75, na.rm = TRUE)
    trimean <- (q1 + 2 * q2 + q3) / 4
    cat("DEBUG: - Trimean:", round(trimean, 3), "\n")
  } else {
    trimean <- NA
    cat("DEBUG: - No shared SNPs for", trait, "\n")
  }
  
  results[[trait]] <- list(
    trait = trait,
    mean_log10_p_in_AD = trimean,
    n_SNPs_shared = nrow(merged)
  )
}

# --- Compile and sort ---
result_df <- rbindlist(results, fill = TRUE)
cat("DEBUG: Compiled result_df with", nrow(result_df), "rows\n")

#with this all traits appear in the plot
all_traits <- data.table(trait = trait_codes)
result_df <- merge(all_traits, result_df, by = "trait", all.x = TRUE)

# Clean trait labels
clean_labels <- function(x) {
  x %>% sub("_\\d+(\\.\\d+)?$", "", .)
}

result_df$clean_trait <- clean_labels(result_df$trait)

# Sorting
result_df <- result_df[order(is.na(result_df$mean_log10_p_in_AD), -result_df$mean_log10_p_in_AD), ]
result_df$clean_trait <- factor(result_df$clean_trait, levels = rev(result_df$clean_trait))

# --- Plot ---
cat("DEBUG: Generating plot...\n")
plt <- ggplot(result_df, aes(x = clean_trait, y = mean_log10_p_in_AD)) +
  geom_bar(stat = "identity", fill = "darkslateblue", na.rm = FALSE) +
  coord_flip() +
  ggtitle(paste0("AD Significance in Trait-Specific SNPs (p < ", SIGNIFICANCE_THRESHOLD, ")")) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10)
  ) +
  labs(
    x = "Trait",
    y = "Tukey Trimean of -log10(p) in AD"
  )

out_path <- file.path("processed_data", paste0("sig_overAD_", SIGNIFICANCE_THRESHOLD, ".png"))
ggsave(out_path, plt, width = 7, height = 6)
cat("DEBUG: Plot saved to", out_path, "\n")
