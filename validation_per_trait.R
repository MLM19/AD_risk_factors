# Required Libraries
library(dplyr)
library(ggplot2)
source("validating_data.R")

# Trait formats (short/long)
trait_formats <- list(
  AD_7158 = "long", ALCOH_3727 = "short", APNEA_7362 = "short", ASBEST_3282 = "short",
  BRE_CANCER_174 = "short", COR_ARTERY_2314 = "long", COVID_0403 = "short", DENTAL_521.1 = "short",
  DEPRES_6475 = "long", DIAB_GEST_6553 = "short", DIAB_TYP1_3791 = "short", DIAB_TYP2_8926 = "long",
  DIET_9364 = "short", EBS_20002 = "short", EDU_1875 = "short", FRECKL_6091 = "short",
  HAIR_COL_1747.4 = "short", HEART_FAIL_2626 = "short", HerpZo_8721 = "short", HHV_6906 = "short",
  HIV_5531 = "short", HSV_5524 = "short", HYPERTEN_9346 = "short", INFLU_2107 = "short",
  INSOM_7286 = "short", IRON_ANEM_280 = "short", LACT_INT_4159 = "short", LE_BO_DEM_1390 = "short",
  LONE_1704 = "short", MENOP_2553 = "long", OBES_5785 = "short", PARKINSON_9325 = "short",
  PIÑA_104570 = "short", PSY_STRESS_4320 = "short", RHINO_4043 = "short", SKIN_COL_9962 = "short",
  SMOK_3686 = "short", STRESS_5871 = "short", STROK_4723 = "short", VAS_DEM_290.16 = "short",
  VITB_5741 = "short", VITB12_5799 = "short", VITD_5742 = "short"
)

# Traits list
traits <- names(trait_formats)

# Function to validate one trait against all others
validate_trait_vs_all <- function(input_trait, input_format) {
  all_traits <- names(trait_formats)
  other_traits <- setdiff(all_traits, input_trait)
  
  for (target_trait in other_traits) {
    target_format <- trait_formats[[target_trait]]
    
    trait_clean_file <- paste0("preprocessed_data/", input_trait, "_effect_snp_clean.tsv")
    trait_significant_file <- paste0("preprocessed_data/", input_trait, "_effect_snp_significant.tsv")
    
    if (input_format == "short" && target_format == "short") {
      ad_trait_file <- paste0("processed_data/", input_trait, "_", target_trait, "_filtered_snps.tsv")
    } else if (input_format == "long" && target_format == "long") {
      ad_trait_file <- paste0("processed_data/", input_trait, "_", target_trait, "_filtered_snps.tsv")
    } else if (input_format == "short" && target_format == "long") {
      ad_trait_file <- paste0("processed_data/", input_trait, "_", target_trait, "_filtered_snps.tsv")
    } else {
      ad_trait_file <- paste0("processed_data/", input_trait, "_", target_trait, "_filtered_snps.tsv")
    }
    
    cat("\n===============================\n")
    cat("Validating", input_trait, "vs", target_trait, "\n")
    cat("===============================\n")
    
    tryCatch({
      validate_snps(
        trait_clean_file = trait_clean_file,
        trait_significant_file = trait_significant_file,
        ad_trait_file = ad_trait_file,
        trait_name = paste(input_trait, "vs", target_trait)
      )
    }, error = function(e) {
      cat("Skipped due to error:", conditionMessage(e), "\n")
    })
  }
}

validate_trait_vs_all(input_trait = "ASBEST_3282", input_format = "short")
