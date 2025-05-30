import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# --- Parameters ---
SIGNIFICANCE_THRESHOLD = 1e-5  #modifiable

trait_codes = [
    "ALCOH_3727", "APNEA_7362", "ASBEST_3282", "BRE_CANCER_174", "COR_ARTERY_2314", "COVID_0403",
    "DENTAL_521.1", "DEPRES_6475", "DIAB_GEST_6553", "DIAB_TYP1_3791", "DIAB_TYP2_8926",
    "DIET_9364", "EBS_20002", "EDU_1875", "FRECKL_6091", "HAIR_COL_1747.4", "HEART_FAIL_2626",
    "HerpZo_8721", "HHV_6906", "HIV_5531", "HSV_5524", "HYPERTEN_9346", "INFLU_2107",
    "INSOM_7286", "IRON_ANEM_280", "LACT_INT_4159", "LE_BO_DEM_1390", "LONE_1704", "MENOP_2553",
    "OBES_5785", "PARKINSON_9325", "PIÑA_104570", "PSY_STRESS_4320", "RHINO_4043", "SKIN_COL_9962",
    "SMOK_3686", "STRESS_5871", "STROK_4723", "VAS_DEM_290.16", "VITB_5741", "VITB12_5799",
    "VITD_5742"
]

print("Loading AD summary statistics...")
ad_df = pd.read_csv("preprocessed_data/AD_7158_effect_snp_clean.tsv", sep="\t", usecols=["rs_id", "p_value"])
ad_df = ad_df.dropna()
ad_df["rs_id"] = ad_df["rs_id"].astype(str)

results = []

# --- Loop through traits ---
for trait in trait_codes:
    file_path = f"preprocessed_data/{trait}_effect_snp_clean.tsv"
    print(f"Processing trait: {trait}")

    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        continue

    # Load full trait data
    trait_df = pd.read_csv(file_path, sep="\t", usecols=["rs_id", "p_value"])
    trait_df["rs_id"] = trait_df["rs_id"].astype(str)
    trait_df = trait_df.dropna()

    # Apply significance threshold
    sig_trait_df = trait_df[trait_df["p_value"] < SIGNIFICANCE_THRESHOLD]
    print(f"- Significant SNPs at p < {SIGNIFICANCE_THRESHOLD}: {len(sig_trait_df)}")

    # Join with AD p-values
    merged = pd.merge(sig_trait_df, ad_df, on="rs_id", how="inner")

    if not merged.empty:
        merged["log_p"] = -np.log10(merged["p_value_y"])  # AD p-values
        q1 = merged["log_p"].quantile(0.25)
        q2 = merged["log_p"].quantile(0.50)
        q3 = merged["log_p"].quantile(0.75)
        trimean = (q1 + 2*q2 + q3) / 4
    else:
        trimean = np.nan


    results.append({
        "trait": trait,
        "mean_log10_p_in_AD": trimean,
        "n_SNPs_shared": len(merged)
    })

# Convert to DataFrame and sort
result_df = pd.DataFrame(results).dropna().sort_values(by="mean_log10_p_in_AD", ascending=False)

# --- Plot ---
print("Generating plot...")
plt.figure(figsize=(8, 14))
bars = plt.barh(result_df["trait"], result_df["mean_log10_p_in_AD"], color="darkslateblue")
plt.xlabel("Tukey Trimean of -log10(p) in AD")
plt.title(f"Tukey Trimean AD Significance of Trait-Specific SNPs (p < {SIGNIFICANCE_THRESHOLD})")
plt.gca().invert_yaxis()  # Highest bar on top
plt.tight_layout()
plt.show()
