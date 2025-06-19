# --- Libraries ---
library(data.table)
library(dplyr)
library(ggplot2)
library(tibble)
library(stringr)
library(smacof)
library(ggrepel)

# --- Parameters ---
SIGNIFICANCE_THRESHOLD <- 1e-5
LD_WINDOW <- 200000
data_dir <- "preprocessed_data"
out_dir <- "processed_data"
trait_codes <- c(
  "AD_7158", "ALCOH_3727", "APNEA_7362", "ASBEST_3282", "BRE_CANCER_174", "COR_ARTERY_2314", "COVID_0403",
  "DENTAL_521.1", "DEPRES_6475", "DIAB_GEST_6553", "DIAB_TYP1_3791", "DIAB_TYP2_8926",
  "DIET_9364", "EBS_20002", "EDU_1875", "FRECKL_6091", "HAIR_COL_1747.4", "HEART_FAIL_2626",
  "HerpZo_8721", "HHV_6906", "HIV_5531", "HSV_5524", "HYPERTEN_9346", "INFLU_2107",
  "INSOM_7286", "IRON_ANEM_280", "LACT_INT_4159", "LE_BO_DEM_1390", "LONE_1704", "MENOP_2553",
  "OBES_5785", "PARKINSON_9325", "PIÑA_104570", "PSY_STRESS_4320", "RHINO_4043", "SKIN_COL_9962",
  "SMOK_3686", "STRESS_5871", "STROK_4723", "VAS_DEM_290.16", "VITB_5741", "VITB12_5799", "VITD_5742"
)

filename_suffix <- paste0("pval_", SIGNIFICANCE_THRESHOLD) #for the files

# --- LD pruning function ---
ld_prune <- function(df, LD_WINDOW = 200000) {
  df <- df[order(df$p_value, df$chr, df$pos), ]
  pruned <- df[1, ]
  for (i in 2:nrow(df)) {
    last_chr <- pruned$chr[nrow(pruned)]
    last_pos <- pruned$pos[nrow(pruned)]
    if (df$chr[i] != last_chr || abs(df$pos[i] - last_pos) > LD_WINDOW) {
      pruned <- rbind(pruned, df[i, ])
    }
  }
  return(pruned)
}

# --- Load and preprocess all traits ---
trait_data <- list()
cat("Loading trait data...\n")
for (trait in trait_codes) {
  file_path <- file.path(data_dir, paste0(trait, "_effect_snp_clean.tsv"))
  if (!file.exists(file_path)) next
  
  df <- fread(file_path, select = c("rs_id", "chr", "pos", "p_value")) %>%
    na.omit() %>%
    mutate(
      chr = as.factor(chr),
      pos = as.numeric(pos),
      p_value = as.numeric(p_value)
    )
  
  trait_data[[trait]] <- df
}

# --- Initialize output matrix ---
n_traits <- length(trait_codes)
logpval_trimean_matrix <- matrix(NA, n_traits, n_traits, dimnames = list(trait_codes, trait_codes))

# --- Fill matrix ---
cat("Processing trait pairs...\n")
for (trait1 in trait_codes) {
  cat(sprintf("Processing trait1: %s\n", trait1))
  
  df1 <- trait_data[[trait1]]
  if (is.null(df1)) {
    cat(sprintf("  → No data for %s, skipping...\n", trait1)) 
    next
  }
  
  sig1 <- df1 %>% filter(p_value <= SIGNIFICANCE_THRESHOLD)
  if (nrow(sig1) == 0) {
    cat(sprintf("  → No significant SNPs for %s, skipping...\n", trait1)) 
    next
  }
  
  sig1_pruned <- if (nrow(sig1) <= 1) sig1 else ld_prune(sig1)
  
  for (trait2 in trait_codes) {
    df2 <- trait_data[[trait2]]
    if (is.null(df2)) next
    
    cat(sprintf("    Comparing to trait2: %s\n", trait2)) 
    
    # Add bias correction for trait2
    df2 <- df2 %>%
      arrange(p_value) %>%
      mutate(pval_order = rank(p_value) / n())  #RANK OR ORDER?????
    
    # Merge by chr + pos
    merged <- inner_join(
      sig1_pruned %>% select(chr, pos),
      df2 %>% select(chr, pos, p_value, pval_order),
      by = c("chr", "pos")
    )
    
    if (nrow(merged) == 0) next
    
    logp <- -log10(merged$pval_order)    #compute trimean using the merged values of the order(pval)
    q1 <- quantile(logp, 0.25, na.rm = TRUE)
    q2 <- quantile(logp, 0.5, na.rm = TRUE)
    q3 <- quantile(logp, 0.75, na.rm = TRUE)
    trimean <- (q1 + 2 * q2 + q3) / 4
    
    logpval_trimean_matrix[trait1, trait2] <- trimean
  }
  
  # Save intermediate matrix
  cat(sprintf("  → Saving matrix after processing %s\n", trait1))
  fwrite(
    as.data.table(logpval_trimean_matrix, keep.rownames = "trait1"),
    file.path(out_dir, paste0("asymmetric_logpval_trimean_matrix_", filename_suffix, ".tsv")),
    sep = "\t"
  )
}

# --- Save matrix at the end---
fwrite(as.data.table(logpval_trimean_matrix, keep.rownames = "trait1"),
       file.path(out_dir, paste0("asymmetric_logpval_trimean_matrix_", filename_suffix, ".tsv")), sep = "\t")

# --- NMDS projection ---
cat("Running NMDS projection...\n")
keep_rows <- rowSums(!is.na(logpval_trimean_matrix)) > 1
keep_cols <- colSums(!is.na(logpval_trimean_matrix)) > 1
filtered_matrix <- logpval_trimean_matrix[keep_rows, keep_cols]


## --- Function to clean labels ---
clean_labels <- function(x) {
  sub("_\\d+(\\.\\d+)?$", "", x)  # Remove trailing _number or _number.decimal
}

# Fill NAs
filled_matrix <- filtered_matrix
na_mean <- mean(filled_matrix, na.rm = TRUE)
filled_matrix[is.na(filled_matrix)] <- na_mean

# Run NMDS
nmds_result <- smacofSym(filled_matrix, ndim = 2, type = "ordinal", ties = "tertiary")

# Save NMDS coordinates
nmds_coords <- as.data.frame(nmds_result$conf)
colnames(nmds_coords) <- c("Dim1", "Dim2")
nmds_coords$trait <- rownames(filled_matrix)
nmds_coords$Label <- clean_labels(nmds_coords$trait)
nmds_coords$IsAD <- nmds_coords$trait == "AD_7158"

#Save coordinates
fwrite(nmds_coords, file.path(out_dir, paste0("nmds_logpval_trimean_coords_", filename_suffix, ".tsv")), sep = "\t")

# --- Plot NMDS ---
p_nmds <- ggplot(nmds_coords, aes(x = Dim1, y = Dim2)) +
  geom_point(aes(color = IsAD), size = 3, alpha = 0.9) +
  scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey50"), guide = FALSE) +
  geom_text_repel(
    aes(label = Label),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    force = 1,
    segment.size = 0.25
  ) +
  theme_bw() +
  labs(
    title = paste("NMDS Projection (", SIGNIFICANCE_THRESHOLD, ") - Tukey Trimean log10(p-values)", sep = ""),
    x = "Dimension 1",
    y = "Dimension 2"
  )

ggsave(file.path(out_dir, paste0("nmds_logpval_trimean_projection_", filename_suffix, ".png")), p_nmds, width = 8, height = 6)