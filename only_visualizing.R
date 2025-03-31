#Script to visualize pre-processed data. 
# MUST HAVE: trait_effect_snp_clean.tsv --> Include columns rs_id, chr, pos, SE, beta, and p_value.

#1. Loading libraries
library(data.table) # Efficient data handling (fread, fwrite)
library(qqman)  # GWAS visualization (manhattan, qq)

#2. Function
##Visualization of data (clean data -> qqplot+manhattan plot)
plot_gwas_results <- function(cleandata, dataset_name, p_value_threshold = 1e-5) {
  if (!all(c("chr", "pos", "p_value") %in% colnames(cleandata))) { # modified chromosome and base_pair_location for chr and pos
    stop("Missing required columns for Manhattan plot!")
  }
  
  # Convert chr to numeric, keeping X as 23 and Y as 24
  if (!is.numeric(cleandata$chr)) {
    message("Converting 'chr' to numeric for Manhattan plot...")
    cleandata[, chr := ifelse(chr == "X", 23, ifelse(chr == "Y", 24, as.numeric(chr)))]
    cleandata <- cleandata[!is.na(chr)]  # Remove any remaining NA values
  }
  
  sig_threshold <- -log10(p_value_threshold)  # GWAS significance threshold (as reference to articles)
  dir.create("Images", showWarnings = FALSE, recursive = TRUE) # Make sure that the directory exists
  manhattan_filename <- paste0("Images/", dataset_name, "_manhattan_plot.png")
  png(manhattan_filename, width = 1200, height = 800, res = 150)  # high resolution PNG needed?
  
  message("Start the creation of the manhattan plot..")
  # Manhattan plot
  manhattan(
    cleandata,
    chr = "chr", # modified chromosome for chr
    bp = "pos",  # modified base_pair_location for pos
    p = "p_value",
    snp = "rs_id",
    logp = TRUE,
    main = paste(dataset_name, "Manhattan Plot"),
    suggestiveline = sig_threshold,
    genomewideline = -log10(5e-8),
    col = c("black", "grey")
  )
  #abline(h = sig_threshold, col = "red", lwd = 2, lty = 2)  # Add significance line
  
  dev.off()
  message("Saved Manhattan plot: ", manhattan_filename)
  
  # QQ plot
  qq_filename <- paste0("Images/", dataset_name, "_qq_plot.png")
  png(qq_filename, width = 1200, height = 800, res = 150)
  message("Start the creation of the qq plot..")
  qq(cleandata$p_value, main = paste(dataset_name, "QQ Plot of P-values"))
  dev.off()
  message("Saved QQ plot: ", qq_filename)
}

#3. Implementation
db_name <- "HerpZo_8721"
message("Uploading preprocessed data...")
cleandata_path <- paste0("preprocessed_data/", db_name, "_effect_snp_clean.tsv")
#Check if file exists to make sure that the data was already preprocessed. 
if (file.exists(cleandata_path)) {
  cleandata <- fread(cleandata_path)
  print(unique(cleandata$chr))
  str(cleandata)
  plot_gwas_results(cleandata, db_name, p_value_threshold = 1e-5)
} else {
  message("File not found: ", cleandata_path)
}
