#Cleaning data code, before merging.

#Set working directory
#setwd("/media/marialm/6E98-2F9B/TFG")   #ubuntu
#setwd("D:/TFG")                        #windows

#1. Load libraries
library(data.table) #Efficient data handling (fread, fwrite)
library(dplyr)  #Data manipulation (mutate, filter, distinct)

#2.Functions
##Loading and validating the data function
load_database <- function(input_path, separator = "\t", column_names) {
  file_name <- basename(input_path)
  message("Loading dataset: ", file_name)
  
  data <- fread(input_path, sep = separator)            #read the data selecting the separator (by default is tab, but can specify space or comma)
  
  missing_cols <- setdiff(column_names, colnames(data))   #check if the db has all the columns that needed (the ones given by column_names argument)
  if (length(missing_cols) > 0) {
    stop("Missing columns in dataset: ", paste(missing_cols, collapse = ", ")) #return a message with the columns that are not found in the db if there are
  }
  
  message("Dataset loaded successfully:)")                                         
  message("Dimensions of loaded database: ", paste(dim(data), collapse = " x "))   #print dimensions to check
  return(data)
}

##Cleaning and analysing the dataset funtion 
clean_and_analyze_database <- function(data, column_mapping, p_value_threshold = 5e-8) {  #column_mapping is a dictionary with the different configurations of the db, and the threshold can be modified aswell.
  message("Cleaning and analyzing data...")
  
  #if beta is not a column then incorporate it by converting odds ratio into a new column
  if (!"beta" %in% colnames(data) && "odds_ratio" %in% colnames(data)) {      #true beta not in data, true odds_ratio in data
    message("Converting odds ratio to beta.")
    data <- data %>% mutate(beta = ifelse(odds_ratio > 0, log(odds_ratio), NA))  #if odds ratio is bigger than 0, calculate the beta=log(odds_ratio) otherwise NA
  }
  
  #Rename the columns based on mapping
  data <- data %>% rename(!!!column_mapping)  #!!! allows dynamic renaming, must have multiple values (columns, threshold, mapping)
  
  #Convert data types
  data <- data %>%       #convert the data types into the ones we are most interested
    mutate(
      effect_allele = toupper(effect_allele), #ensure capital letters for alleles so that it can be analysed
      alt_allele = toupper(alt_allele),
      p_value = as.numeric(p_value), 
      chr = as.factor(chr),
      pos = as.numeric(pos),   #nt position must be numeric
      beta = as.numeric(beta),
      SE = as.numeric(SE)
    )
  
  #Filtering out 
  data <- data %>% select(where(~ !all(is.na(.))))  #remove the columns that only contain NA as data
  data <- na.omit(data)               #Remove rows with NA values
  #data <- data %>% filter(p_value <= p_value_threshold) #Filter the data with p-value smaller than given p-value threshold
  data <- data %>% filter(nchar(effect_allele) == 1 & nchar(alt_allele) == 1) #Not interested in indels, only one allele must be affected
  
  transitions <- c("A,G", "G,A", "C,T", "T,C")   #Filter out transversions, by defining transitions
  data <- data %>% 
    filter(paste(effect_allele, alt_allele, sep = ",") %in% transitions)   #Take the values of these two columns and check if they are on the list of transitions, then keep them if they are.
  
  message("Data cleaning and analysis complete.")
  message("Dimensions of cleaned database: ", paste(dim(data), collapse = " x "))  #State dimensions to check that the filtering is well done
  
  #Stop if the cleaned data has no rows or columns
  if (dim(data)[1] == 0 || dim(data)[2] == 0) {           
    stop("Error: No SNPs for which we are interested --> Stopping further processing.")
  }
  
  return(data) #provide the clean data
}

##Quality check function
quality_check <- function(data) {
  message("Performing quality check...")   #as a marker to know in which part of the process the program is
  
  message("Structure of the data:")
  print(str(data))       #Check the data structure
  
  message("Summary of the data:")
  print(summary(data))    #Check some summaries, to visually check that they make sense
  
  na_count <- colSums(is.na(data))    #Make sure there are no NA that can difficult the processing later on (supposed to be filtered out on the cleaning function)
  if (sum(na_count) > 0) {
    message("Warning: NA values found in the following columns:")
    print(na_count[na_count > 0])
  } else {
    message("No NA values found in the data.")
  }
  
  # Perform transformations if needed
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

#Dictionary for column mappings and separators -> it contains the names of the columns we are to find, the separator that db configuration has and also a list of the names of the columns to do the mapping later on
#This allows the program to dynamically detect what configuration the database has, and adapt to processing it that way. 
database_configs <- list(
  db1 = list(columns = c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"), sep = " ",
             mapping = list(rs_id = "SNP", chr = "CHR", pos = "BP", effect_allele = "A1", alt_allele = "A2", beta = "BETA", SE = "SE", p_value = "P")),
  db2 = list(columns = c("rs_id", "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "p_value"), sep = "\t", 
             mapping = list(rs_id = "rs_id", chr = "chromosome", pos = "base_pair_location", effect_allele = "effect_allele", alt_allele = "other_allele", beta = "beta", SE = "standard_error", p_value = "p_value")),
  db3 = list(columns = c("MarkerName", "CHR", "POS", "A1", "A2", "Beta", "SE", "Pval"), sep = "\t", 
             mapping = list(rs_id = "MarkerName", chr = "CHR", pos = "POS", effect_allele = "A1", alt_allele = "A2", beta = "Beta", SE = "SE", p_value = "Pval")),
  db4 = list(columns = c("hm_rsid", "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "p_value"), sep = "\t", 
             mapping = list(rs_id = "hm_rsid", chr = "chromosome", pos = "base_pair_location", effect_allele = "effect_allele", alt_allele = "other_allele", beta = "beta", SE = "standard_error", p_value = "p_value")),
  db5 = list(columns = c("rsid", "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "p_value"), sep = "\t",
             mapping = list(rs_id = "rsid", chr = "chromosome", pos = "base_pair_location", effect_allele = "effect_allele", alt_allele = "other_allele", beta = "beta", SE = "standard_error", p_value = "p_value")),
  db6 = list(columns = c("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele", "odds_ratio", "standard_error", "p_value"), sep = "\t",
             mapping = list(rs_id = "variant_id",chr = "chromosome",pos = "base_pair_location",effect_allele = "effect_allele",alt_allele = "other_allele",beta = "beta",SE = "standard_error",p_value = "p_value")),
  db7 = list(columns = c("MarkerName", "Chromosome", "Position", "Effect_allele", "Non_Effect_allele", "Beta", "SE", "Pvalue"), sep = " ", 
             mapping = list(rs_id = "MarkerName",chr = "Chromosome",pos = "Position",effect_allele = "Effect_allele",alt_allele = "Non_Effect_allele",beta = "Beta",SE = "SE",p_value = "Pvalue")),
  db8 = list(columns = c("Name", "chromosome", "base_pair_location", "effect_allele", "other_allele", "odds_ratio", "standard_error", "p_value"), sep = "\t",
             mapping = list(rs_id = "Name",chr = "chromosome",pos = "base_pair_location",effect_allele = "effect_allele",alt_allele = "other_allele",beta = "beta",SE = "standard_error",p_value = "p_value")), 
  db9 = list(columns = c("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "p_value"), sep = "\t",
             mapping = list(rs_id = "variant_id",chr = "chromosome",pos = "base_pair_location",effect_allele = "effect_allele",alt_allele = "other_allele",beta = "beta",SE = "standard_error",p_value = "p_value"))
)

##Detect the database configuration function
detect_database_config <- function(data) {
  for (db_name in names(database_configs)) {   #check each of the configurations in the dictionary 
    config <- database_configs[[db_name]]      #select one configuration to start comparing (that way we have the object and can use it's key values)
    if (all(config$columns %in% colnames(data))){  #Check that all the columns of the selected configuration are in the column names of our data
      message("Detected database: ", db_name)   #We have selected a configuration that works for our database
      return(config)
    }
  }
  stop("No matching database structure found:(")  #There are no configurations that adapt to our database --> must add the new configuration manually
}


##MAIN FUNCTION
process_gwas <- function(input_path, db_name, separator = "\t", p_value_threshold = 5e-8) { #by default, separator is going to be a tab (more common), and threshold is going to be 5e-8 (can be changed)
  message("Processing GWAS dataset...")
  
  #Loads first row of data to be able to check its configuration, with the given separator
  data <- fread(input_path, sep = separator, nrows = 1) 
  if (ncol(data)==1) {        #If the data is loaded with 1 column, means that the separator is not good. We should have at least 8 columns
    stop("The separator seems incorrect, must be changed")   #Stop the process to change the separator
  }
  config <- detect_database_config(data) #check for the best configuration for the db
  
  #Properly load the data, having checked the proper separator and configuration
  data <- load_database(input_path, separator = separator, column_names = config$columns)
  
  ##Use functions to clean and check the data
  cleaned_data <- clean_and_analyze_database(data, config$mapping, p_value_threshold) #colum names are found in the config of the db
  final_data <- quality_check(cleaned_data)
  
  ##Generate output files
  #Create the names of the files
  snp_table_file <- paste0("preprocessed_data/", db_name, "_snp_clean.tsv")
  effect_snp_table_file <- paste0("preprocessed_data/", db_name, "_effect_snp_clean.tsv")
  
  #Keep the data from the preprocessed data
  snp_table <- final_data %>% select(rs_id, chr, pos, effect_allele, alt_allele)
  effect_snp_table <- final_data %>% select(rs_id, chr, pos, SE, beta, p_value)
  
  #Save the results into files
  fwrite(snp_table, snp_table_file, sep = "\t")
  fwrite(effect_snp_table, effect_snp_table_file, sep = "\t")
  
  #Tracker of the work being well done
  message("Results saved in generated files: ")
  message("- ", snp_table_file)
  message("Dimensions: ", paste(dim(snp_table), collapse = " x "))
  message("- ", effect_snp_table_file)
  message("Dimensions: ", paste(dim(effect_snp_table), collapse = " x "))
  
  return(list(
    cleaned_data = final_data,
    snp_table = snp_table,
    effect_snp_table = effect_snp_table
  ))
}

#3. Implementation
input_path <- "Raw_Data/GCST90435742.h.tsv.gz"
db_name <- "VITD_5742"
processed_data <- process_gwas(input_path, db_name, separator = "\t", p_value_threshold = 10e-5)
