# preprocessing all db + manhattan, qqplot but more fast and efficient (FIXING PROBLEMS FOUND BEFORE)

#local setwd("C:/Users/Usuario/Downloads")

#1. Load libraries
library(data.table) # Efficient data handling (fread, fwrite)
library(dplyr)  # Data manipulation (mutate, filter, distinct)
library(R.utils) # For gunzip (files ending .bgz)
library(qqman)  # GWAS visualization (manhattan, qq)
library(ggplot2) # General plotting (hist, abline, text)

#2. Defining functions for preprocessing

##Loading and validating the data function
load_database <- function(input_path, separator = "\t", column_names) {
  file_name <- basename(input_path) #Since the input_path can be a full path, we get the name of the file to print it later on
  
  #If file is compressed in .bgz, we decompress it first
  if (grepl("\\.bgz$", file_name)) {
    message("Detected .bgz file. Decompressing...")
    decompressed_path <- sub("\\.bgz$", "", input_path)  #Remove .bgz extension
    if (file.exists(decompressed_path)) { #Once the file is decompressed, changing the name is necessary since it is not reversible
        stop("Error: The decompressed file already exists. Please rename or remove the existing file.")
    }
    gunzip(input_path, destname = decompressed_path, overwrite = FALSE)  #Decompress the file
    input_path <- decompressed_path  #Update the path to the decompressed file so that the program can find the file
  }
  
  message("Loading dataset: ", file_name) #tracker to know which file is being loaded
  
  data <- fread(input_path, sep = separator, data.table = FALSE) #Read the data selecting the separator (by default is tab, but can specify space or comma)
  
  missing_cols <- setdiff(column_names, colnames(data))   #Check if the db has all the columns that needed (the ones given by column_names argument)
  if (length(missing_cols) > 0) {
    stop("Missing columns in dataset: ", paste(missing_cols, collapse = ", ")) #DEBUG:Return a message with the columns that are not found in the db if there are
  }
  
  message("Dataset loaded successfully:)")
  message("Dimensions of loaded database: ", paste(dim(data), collapse = " x "))   #DEBUG: print dimensions to check
  return(data)
}

##Cleaning and analysing the dataset function
clean_and_analyze_database <- function(data) {
  message("Cleaning and analyzing data...") #Marker
  
  message("Initial columns in the data:") #DEBUG: print the columns before renaming
  print(colnames(data))
  
  #if beta is not a column then incorporate it by converting odds ratio into a new column
  if (!"beta" %in% colnames(data) && "odds_ratio" %in% colnames(data)) {
    message("Converting odds ratio to beta.")
    data <- data %>% mutate(beta = ifelse(odds_ratio > 0, log(odds_ratio), NA))
  }
  
  message("Checking if the conversion has introduced NA: ") 
  message("The number of NA values: ", sum(is.na(data))) #DEBUG: since we already now the NA number of values, this way we can see if the transformation went right
  
  #Convert data types to have consistent formats 
  data <- data %>%
    mutate(
      effect_allele = toupper(effect_allele),
      alt_allele = toupper(alt_allele),
      p_value = as.numeric(p_value),
      chr = as.factor(chr),
      pos = as.numeric(pos),
      beta = as.numeric(beta),
      SE = as.numeric(SE)
    )
  
  #Filtering out
  data <- data %>% dplyr::select(where(~ !all(is.na(.))))  #Remove the columns that only contain NA as data
  #Remove rows with NA values of the essential fields (columns )
  message("Columns after removing NA-only columns: ", paste(colnames(data), collapse = ", ")) #DEBUG: print the columns after removing NA-only columns
  data <- data %>% filter(!is.na(rs_id) & !is.na(chr) & !is.na(pos) &  
                            !is.na(effect_allele) & !is.na(alt_allele) &
                            !is.na(beta) & !is.na(SE) & !is.na(p_value))
  
  data <- data %>% filter(nchar(effect_allele) == 1 & nchar(alt_allele) == 1) #Filter data -> Not interested in indels, keep SNPs with single nucleotide alleles
  
  message("Control check after filtering NA and indels: ") #DEBUG: control checks for the quality of data
  message("The number of NA values: ", sum(is.na(data)))
  message("length of the data: ", nrow(data))
  
  transitions <- c("A,G", "G,A", "C,T", "T,C")   #Filter out transitions
  data <- data %>%
    filter(paste(effect_allele, alt_allele, sep = ",") %in% transitions)   #Take the values of these two columns and check if they are on the list of transitions, then keep them if they are.
  
  message("Control check after filtering out transversions: ") #DEBUG: control checks for the quality of data
  message("The number of NA values: ", sum(is.na(data)))
  message("length of the data: ", nrow(data))
  
  
  message("Data cleaning and analysis complete.")
  message("Dimensions of cleaned database: ", paste(dim(data), collapse = " x "))  # State dimensions to check that the filtering is well done
  
  #Stop if the cleaned data has no rows or columns, no valid data
  if (dim(data)[1] == 0 || dim(data)[2] == 0) {
    stop("Error: No SNPs for which we are interested --> Stopping further processing.")
  }
  return(data) #provide the clean data
}

##Quality check function
quality_check <- function(data) {
  message("Performing quality check...")   #Marker: to know in which part of the process the program is
  
  message("Structure of the data:")
  message(str(data))       #DEBUG: Check the data structure
  
  message("Summary of the data:")
  message(summary(data))    #DEBUG: Check some summaries, to visually check that they make sense
  
  na_count <- colSums(is.na(data))    #DEBUG: Make sure there are no NA that can difficult the processing later on (supposed to be filtered out on the cleaning function)
  if (sum(na_count) > 0) {
    message("Warning: NA values found in the following columns:")
    message(na_count[na_count > 0])
  } else {
    message("No NA values found in the data.")
  }
  
  #Perform transformations if needed
  if (!is.factor(data$chr)) {
    message("Converting 'chr' to factor...")
    data$chr <- as.factor(data$chr)          #if the structure was not correct then transform it into the correct one for chr, pos, beta, SE and the p-value
  }
  if (!is.numeric(data$pos)) {
    message("Converting 'pos' to numeric...")
    data$pos <- as.numeric(data$pos)
  }
  if (!is.numeric(data$beta)) {
    message("Converting 'beta' to numeric...")
    data$beta <- as.numeric(data$beta)
  }
  if (!is.numeric(data$SE)) {
    message("Converting 'SE' to numeric...")
    data$SE <- as.numeric(data$SE)
  }
  if (!is.numeric(data$p_value)) {
    message("Converting 'p_value' to numeric...")
    data$p_value <- as.numeric(data$p_value)
  }
  
  message("Quality check complete.") #tracker
  return(data)
}

##Detect the database configuration function
detect_database_config <- function(data) {
  for (db_name in names(database_configs)) {   #check each of the configurations in the dictionary
    config <- database_configs[[db_name]]      #select one configuration to start comparing (that way we have the object and can use it's key values)
    if (all(config$columns %in% colnames(data))) {  #Check that all the columns of the selected configuration are in the column names of our data
      message("Detected database: ", db_name)   #We have selected a configuration that works for our database
      return(config)
    }
  }
  stop("No matching database structure found:(")  #There are no configurations that adapt to our database --> must add the new configuration manually
}

##Dictionary for column mappings and separators -> it contains the names of the columns we are to find, the separator that db configuration has and also a list of the names of the columns to do the mapping later on
#This allows the program to dynamically detect what configuration the database has, and adapt to processing it that way.
database_configs <- list(
  #Each db defines the requiered columns, the separator used in the db and a mapping the columns to standard field names
  db1 = list(columns = c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"), sep = " ",
             mapping = list(rs_id = "SNP", chr = "CHR", pos = "BP", effect_allele = "A1", alt_allele = "A2", beta = "BETA", SE = "SE", p_value = "P")),
  db2 = list(columns = c("rs_id", "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "p_value"), sep = "\t",
             mapping = list(rs_id = "rs_id", chr = "chromosome", pos = "base_pair_location", effect_allele = "effect_allele", alt_allele = "other_allele", beta = "beta", SE = "standard_error", p_value = "p_value")),
  db3 = list(columns = c("MarkerName", "CHR", "POS", "A1", "A2", "Beta", "SE", "Pval"), sep = "\t",
             mapping = list(rs_id = "MarkerName", chr = "CHR", pos = "POS", effect_allele = "A1", alt_allele = "A2", beta = "Beta", SE = "SE", p_value = "Pval")),
  db4 = list(columns = c("hm_rsid", "hm_chrom","hm_pos", "hm_other_allele", "hm_effect_allele", "hm_beta", "standard_error", "p_value"), sep = "\t",
             mapping = list(rs_id = "hm_rsid", chr = "hm_chrom", pos = "hm_pos", effect_allele = "hm_effect_allele", alt_allele = "hm_other_allele", beta = "hm_beta", SE = "standard_error", p_value = "p_value")),
  db5 = list(columns = c("rsid", "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "p_value"), sep = "\t",
             mapping = list(rs_id = "rsid", chr = "chromosome", pos = "base_pair_location", effect_allele = "effect_allele", alt_allele = "other_allele", beta = "beta", SE = "standard_error", p_value = "p_value")),
  db6 = list(columns = c("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele", "odds_ratio", "standard_error", "p_value"), sep = "\t",
             mapping = list(rs_id = "variant_id", chr = "chromosome", pos = "base_pair_location", effect_allele = "effect_allele", alt_allele = "other_allele", beta = "beta", SE = "standard_error", p_value = "p_value")),
  db7 = list(columns = c("MarkerName", "Chromosome", "Position", "Effect_allele", "Non_Effect_allele", "Beta", "SE", "Pvalue"), sep = " ",
             mapping = list(rs_id = "MarkerName", chr = "Chromosome", pos = "Position", effect_allele = "Effect_allele", alt_allele = "Non_Effect_allele", beta = "Beta", SE = "SE", p_value = "Pvalue")),
  db8 = list(columns = c("Name", "chromosome", "base_pair_location", "effect_allele", "other_allele", "odds_ratio", "standard_error", "p_value"), sep = "\t",
             mapping = list(rs_id = "Name", chr = "chromosome", pos = "base_pair_location", effect_allele = "effect_allele", alt_allele = "other_allele", beta = "beta", SE = "standard_error", p_value = "p_value")),
  db9 = list(columns = c("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "p_value"), sep = "\t",
             mapping = list(rs_id = "variant_id", chr = "chromosome", pos = "base_pair_location", effect_allele = "effect_allele", alt_allele = "other_allele", beta = "beta", SE = "standard_error", p_value = "p_value")),
  db10 = list(columns = c("chr", "pos", "ref", "alt", "beta_EUR", "se_EUR", "neglog10_pval_EUR"), sep = "\t",
              mapping = list(chr = "chr", pos = "pos", effect_allele = "ref", alt_allele = "alt", beta = "beta_EUR", SE = "se_EUR", p_value = "p_value")),
  db11 = list(columns = c("rsid", "#CHR", "POS", "REF", "ALT", "all_inv_var_meta_beta", "all_inv_var_meta_sebeta", "all_inv_var_meta_p"), sep = "\t",
              mapping = list(rs_id = "rsid", chr = "#CHR", pos = "POS", effect_allele = "REF", alt_allele = "ALT", beta = "all_inv_var_meta_beta", SE = "all_inv_var_meta_sebeta", p_value = "all_inv_var_meta_p")), 
  db12 = list(columns = c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "p_value"), sep = "\t",
              mapping = list( chr = "chromosome", pos = "base_pair_location", effect_allele = "effect_allele", alt_allele = "other_allele", beta = "beta", SE = "standard_error", p_value = "p_value")),
  db13 = list(columns = c("snpid", "chr", "bpos", "a1", "a2", "mtag_beta", "mtag_se", "mtag_pval"), sep = "\t",
             mapping = list(rs_id = "snpid", chr = "chr", pos = "bpos", effect_allele = "a1", alt_allele = "a2", beta = "mtag_beta", SE = "mtag_se", p_value = "mtag_pval"))
)

##Visualization of data (clean data -> qqplot+manhattan plot)
plot_gwas_results <- function(cleandata, dataset_name, p_value_threshold = 1e-5) {
  if (!all(c("chr", "pos", "p_value") %in% colnames(cleandata))) { # modified chromosome and base_pair_location for chr and pos
    stop("Missing required columns for Manhattan plot!")
  }
  
  #Convert chr to numeric, keeping X as 23 and Y as 24
  if (!is.numeric(cleandata$chr)) {
    message("Converting 'chr' to numeric for Manhattan plot...")
    cleandata[, chr := fifelse(chr == "X", 23,  #Using fifelse(), optimized and prevents unnecessary coercion warnings
                               fifelse(chr == "Y", 24, 
                                       suppressWarnings(as.numeric(chr))))] 
    cleandata <- cleandata[!is.na(chr) & is.finite(chr)]   # Remove any remaining NA values
  }
  
  #Ensure p_value is numeric
  cleandata[, p_value := as.numeric(p_value)]
  
  #Remove problematic p-values: filter out invalid ones
  cleandata <- cleandata[!is.na(p_value) & p_value > 0]
  
  message("Verify the conversion of the chromosome column")
  print(unique(cleandata$chr)) #DEBUG: print the unique values of the chromosome column to check if the conversion was done correctly
  
  sig_threshold <- -log10(p_value_threshold)  #GWAS significance threshold chosen 
  dir.create("Images", showWarnings = FALSE, recursive = TRUE) #Make sure that the directory exists
  manhattan_filename <- paste0("Images/", dataset_name, "_manhattan_plot.png")
  png(manhattan_filename, width = 1200, height = 800, res = 150)  #high resolution PNG needed
  
  message("Start the creation of the manhattan plot..")
  #Manhattan plot: to check the distribution of p-values across chromosomes
  manhattan(
    cleandata,
    chr = "chr", # modified chromosome for chr
    bp = "pos",  # modified base_pair_location for pos
    p = "p_value",
    snp = "rs_id",
    logp = TRUE,
    main = paste(dataset_name, "Manhattan Plot"),
    suggestiveline = sig_threshold, #line for the chosen significance threshold (red line)
    genomewideline = -log10(5e-8), #standard suggested line (blue line)
    col = c("black", "grey")
  )
  
  dev.off()
  message("Saved Manhattan plot: ", manhattan_filename)
  
  #QQ plot: to check if the p-values follow a uniform distribution
  qq_filename <- paste0("Images/", dataset_name, "_qq_plot.png")
  png(qq_filename, width = 1200, height = 800, res = 150)
  message("Start the creation of the qq plot..")
  qq(cleandata$p_value, main = paste(dataset_name, "QQ Plot of P-values"))
  dev.off()
  message("Saved QQ plot: ", qq_filename)
}

##MAIN FUNCTION##
process_gwas <- function(input_path, db_name, separator = "\t", p_value_threshold = 1e-5) { # by default, separator is going to be a tab (more common), and threshold is going to be 5e-8 (can be changed)
  message("Processing GWAS dataset...")
  
  #Loads first row of data to be able to check its configuration, with the given separator
  data_temp <- fread(input_path, sep = separator, nrows = 1, data.table = FALSE)
  if (ncol(data_temp) == 1) {        # If the data is loaded with 1 column, means that the separator is not good. We should have at least 8 columns
    stop("The separator seems incorrect, must be changed")   # Stop the process to change the separator
  }
  config <- detect_database_config(data_temp) # check for the best configuration for the db
  
  #Properly load the data, having checked the proper separator and configuration
  data <- load_database(input_path, separator = config$sep, column_names = config$columns)
  
  message("Control check of raw data: ")
  message("The number of NA values: ", sum(is.na(data)))
  message("Dimensions of loaded database: ", paste(dim(data), collapse = " x "))
  message("Columns before renaming: ")
  print(colnames(data))
  
  #Rename columns based on the mapping in database_configs
  for (new_name in names(config$mapping)) {
    old_name <- config$mapping[[new_name]]
    if (old_name %in% colnames(data)) {
      names(data)[names(data) == old_name] <- new_name
    }
  }
  
  message("Columns after renaming: ")
  print(colnames(data))
  
  #transformations for the columns
  setDT(data) #format data.table to use :=
  #if beta is not a column then incorporate it by converting odds ratio into a new column
  if (!"beta" %in% colnames(data) && "odds_ratio" %in% colnames(data)) {
    message("Converting odds ratio to beta.")
    data <- data %>% mutate(beta = ifelse(odds_ratio > 0, log(odds_ratio), NA))
  }
  
  ##Create rs_id column if it doesn't exist and required columns are present
  if (!"rs_id" %in% colnames(data) && all(c("chr", "pos", "effect_allele") %in% colnames(data))) {
    message("Creating rs_id column from chr, pos, and effect_allele.")
    data[, rs_id := paste(chr, pos, effect_allele, sep = "_")]
  }
  
  ##Compute p-value, handling potential absence of neglog10_pval_EUR
  if ("neglog10_pval_EUR" %in% colnames(data)) {
    message("Computing p_value from neglog10_pval_EUR.")
    data[, p_value := 10^(-neglog10_pval_EUR)]
    data[, neglog10_pval_EUR := NULL]  # Remove the original column
  } else if (!"p_value" %in% colnames(data)) {
    stop("No p_value column found, and neglog10_pval_EUR not present to compute it.")
  }
  
  message("Columns after transformations (if any): ")
  print(colnames(data))
  message("The number of NA values (to check if induced NAs): ", sum(is.na(data)))
  
  #Select only the essential columns
  essential_columns <- c("rs_id", "chr", "pos", "effect_allele", "alt_allele", "beta", "SE", "p_value")
  data <- data %>% dplyr::select(all_of(essential_columns))
  
  message("Columns after renaming and selection:")
  print(colnames(data))
  message("Dimensions after selection: ", paste(dim(data), collapse = " x "))
  message("The number of NA values after selection: ", sum(is.na(data)))
  
  #Use functions to clean and check the data
  cleaned_data <- clean_and_analyze_database(data) # colum names are found in the config of the db
  final_data <- quality_check(cleaned_data)
  
  #Filter by significant p-value
  significant_data <- final_data %>% filter(p_value <= p_value_threshold)
  message("Dimensions of significant data: ", paste(dim(significant_data), collapse = " x "))
  
  ##Generate output files
  #Create the names of the files
  dir.create("preprocessed_data", showWarnings = FALSE, recursive = TRUE)
  snp_table_file <- paste0("preprocessed_data/", db_name, "_snp_clean.tsv")
  effect_snp_table_file <- paste0("preprocessed_data/", db_name, "_effect_snp_clean.tsv")
  
  #Keep the data from the preprocessed data
  snp_table <- final_data %>% dplyr::select(rs_id, chr, pos, effect_allele, alt_allele)
  effect_snp_table <- final_data %>% dplyr::select(rs_id, chr, pos, SE, beta, p_value)
  
  #Save the results into files
  fwrite(snp_table, snp_table_file, sep = "\t")
  fwrite(effect_snp_table, effect_snp_table_file, sep = "\t")
  
  #Tracker of the work being well done
  message("Results saved in generated files: ")
  message("- ", snp_table_file)
  message("Dimensions: ", paste(dim(snp_table), collapse = " x "))
  message("- ", effect_snp_table_file)
  message("Dimensions: ", paste(dim(effect_snp_table), collapse = " x "))
  
  #Check if there are significant SNPs before writing the files
  if (nrow(significant_data) > 0) {
    #Create file names for significant data
    snp_significant_file <- paste0("preprocessed_data/", db_name, "_snp_significant.tsv")
    effect_snp_significant_file <- paste0("preprocessed_data/", db_name, "_effect_snp_significant.tsv")
    
    #Extract the data from significant_data
    snp_significant <- significant_data %>% dplyr::select(rs_id, chr, pos, effect_allele, alt_allele)
    effect_snp_significant <- significant_data %>% dplyr::select(rs_id, chr, pos, SE, beta, p_value)
    
    #Save significant data to files
    fwrite(snp_significant, snp_significant_file, sep = "\t")
    fwrite(effect_snp_significant, effect_snp_significant_file, sep = "\t")
    
    #Message that significant data has been saved
    message("- ", snp_significant_file)
    message("Dimensions: ", paste(dim(snp_significant), collapse = " x "))
    message("- ", effect_snp_significant_file)
    message("Dimensions: ", paste(dim(effect_snp_significant), collapse = " x "))
  } else {
    message("No significant SNPs found at p-value threshold = ", p_value_threshold)
  }
  
  message("Processing complete.")
  
  return(list(final_data = final_data, db_name = db_name)) # Return the final_data and db_name
}


#3. Example usage of the functions
input_path <- "Raw_Data/GCST90435524.h.tsv.gz"
db_name <- "HSV_5524"
processed_data <- process_gwas(input_path, db_name, separator = "\t", p_value_threshold = 1e-5)

# To visualize the data after processing the data
if (!is.null(processed_data)) {
  plot_gwas_results(processed_data$final_data, processed_data$db_name)
} else {
  message("No data to visualize.")
}

