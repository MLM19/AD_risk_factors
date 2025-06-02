# fix_and_rescale_matrix.py

import pandas as pd

# === CONFIGURATION ===
input_file = "processed_data/significant_snp_matrix.tsv"
output_file = "rescaled_significant_snp_matrix.tsv"
significance_threshold = 1e-6  # Not directly used but could be made dynamic

# === STEP 1: Load and clean matrix ===
matrix = pd.read_csv(input_file, sep="\t", index_col=0)

#if needed: Remove any row or column with "_effect"
#matrix = matrix[~matrix.index.str.contains("_effect")]
#matrix = matrix.loc[:, ~matrix.columns.str.contains("_effect")]

# === STEP 2: Normalize using max SNP hits per trait ===
maxhit_trait = {
    "AD_7158": 5836, "ALCOH_3727": 25, "APNEA_7362": 3225, "ASBEST_3282": 71,
    "BRE_CANCER_174": 2026, "COR_ARTERY_2314": 17168, "COVID_0403": 1896,
    "DENTAL_521.1": 38, "DEPRES_6475": 7522, "DIAB_GEST_6553": 9,
    "DIAB_TYP1_3791": 18866, "DIAB_TYP2_8926": 20092, "DIET_9364": 16312,
    "EBS_20002": 11, "EDU_1875": 15945, "FRECKL_6091": 192,
    "HAIR_COL_1747.4": 28487, "HEART_FAIL_2626": 1607, "HerpZo_8721": 17,
    "HHV_6906": 102, "HIV_5531": 10, "HSV_5524": 10, "HYPERTEN_9346": 27580,
    "INFLU_2107": 61, "INSOM_7286": 1113, "IRON_ANEM_280": 828,
    "LACT_INT_4159": 160, "LE_BO_DEM_1390": 167, "LONE_1704": 1139,
    "MENOP_2553": 24234, "OBES_5785": 508, "PARKINSON_9325": 3343,
    "PIÑA_104570": 11985, "PSY_STRESS_4320": 74, "RHINO_4043": 2,
    "SKIN_COL_9962": 292, "SMOK_3686": 23, "STRESS_5871": 6,
    "STROK_4723": 852, "VAS_DEM_290.16": 49, "VITB_5741": 122,
    "VITB12_5799": 52, "VITD_5742": 6
}


# Replace zero-hit traits with small epsilon to avoid divide-by-zero
rescaled_matrix = matrix.div(
    [maxhit_trait.get(trait, 1e-6) if maxhit_trait.get(trait, 0) > 0 else 1e-6
     for trait in matrix.index],
    axis=0
)

# === STEP 3: Remove rows or columns that are all zeros ===
row_sums = rescaled_matrix.sum(axis=1)
col_sums = rescaled_matrix.sum(axis=0)

rows_to_keep = row_sums[row_sums != 0].index
cols_to_keep = col_sums[col_sums != 0].index

removed_rows = set(rescaled_matrix.index) - set(rows_to_keep)
removed_cols = set(rescaled_matrix.columns) - set(cols_to_keep)

if removed_rows:
    print("Removed rows (all zero):", sorted(removed_rows))
if removed_cols:
    print("Removed columns (all zero):", sorted(removed_cols))

rescaled_matrix = rescaled_matrix.loc[rows_to_keep, cols_to_keep]

# === STEP 4: Save cleaned and rescaled matrix ===
rescaled_matrix.to_csv(output_file, sep="\t")
print(f"✅ Cleaned and rescaled matrix saved to: {output_file}")
