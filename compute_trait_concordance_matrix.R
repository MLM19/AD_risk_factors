library(data.table)

# --- Helper Functions ---
trimean <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  return((q[1] + 2 * q[2] + q[3]) / 4)
}

rev_allele <- function(allele) {
  return(chartr("ACGT", "TGCA", allele))
}

logSign <- function(beta1, beta2, se2, a1, a2, oa1, oa2) {
  pvalB <- (beta2 / se2)^2
  signA <- ifelse(beta1 != 0, sign(beta1), 1)
  signB <- ifelse(beta2 != 0, sign(beta2), 1)
  sign_total <- signA * signB
  v <- pvalB * sign_total
  
  if (a1 == a2 || a1 == rev_allele(a2)) {
    return(v)
  } else if (a1 == oa2 || a1 == rev_allele(oa2)) {
    return(-v)
  } else {
    return(NA)
  }
}

filter_distanced_snps <- function(df, threshold = 100000) {
  df <- df[order(chr, pos)]
  df <- df[!duplicated(floor(pos / threshold))]
  return(df)
}

# --- Load Trait Metadata ---
phenotype_table <- fread("GWASPhenotypesTable_Litte.txt")
db_table <- fread("DbTable.txt")
trait_names <- phenotype_table$short

# --- Initialize Matrices ---
n <- length(trait_names)
z_matrix <- matrix(NA, n, n, dimnames = list(trait_names, trait_names))
m_matrix <- matrix(NA, n, n, dimnames = list(trait_names, trait_names))
ne_matrix <- matrix(NA, n, n, dimnames = list(trait_names, trait_names))

# --- Main Loop ---
for (i in seq_along(trait_names)) {
  for (j in seq_along(trait_names)) {
    if (i == j) next
    
    trait1 <- trait_names[i]
    trait2 <- trait_names[j]
    
    file1 <- paste0("preprocessed_data/", trait1, "_effect_snp_clean.tsv")
    file2 <- paste0("preprocessed_data/", trait2, "_effect_snp_clean.tsv")
    
    if (!file.exists(file1) || !file.exists(file2)) next
    
    df1 <- fread(file1, colClasses = list(character = c("rs_id", "effect_allele", "alt_allele")))
    df2 <- fread(file2, colClasses = list(character = c("rs_id", "effect_allele", "alt_allele")))
    
    # Merge by position
    merged <- merge(df1, df2, by = c("chr", "pos"), suffixes = c(".1", ".2"))
    if (nrow(merged) < 10) next
    
    merged <- filter_distanced_snps(merged)
    if (nrow(merged) < 10) next
    
    log_sign_scores <- mapply(logSign,
      merged$beta.1, merged$beta.2, merged$SE.2,
      merged$effect_allele.1, merged$effect_allele.2,
      merged$alt_allele.1, merged$alt_allele.2
    )
    
    log_sign_scores <- log_sign_scores[is.finite(log_sign_scores)]
    if (length(log_sign_scores) < 10) next
    
    # Compute summary stat
    obs_stat <- trimean(log_sign_scores)
    
    # Permutation test
    reps <- 1000
    perm_stats <- replicate(reps, trimean(sample(log_sign_scores, replace = TRUE)))
    p_val <- mean(abs(perm_stats) >= abs(obs_stat))
    
    # Fill matrices
    z_matrix[i, j] <- obs_stat
    m_matrix[i, j] <- p_val
    ne_matrix[i, j] <- length(log_sign_scores)
  }
}

# --- Save Matrices ---
fwrite(as.data.table(z_matrix, keep.rownames = TRUE), "processed_data/z_score_matrix.tsv", sep = "\t")
fwrite(as.data.table(m_matrix, keep.rownames = TRUE), "processed_data/p_value_matrix.tsv", sep = "\t")
fwrite(as.data.table(ne_matrix, keep.rownames = TRUE), "processed_data/snp_count_matrix.tsv", sep = "\t")
