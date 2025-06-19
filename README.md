
![Logo](https://eu01.edcwb.com/buscador/img/centros/logogrande/54808-a088d273bd364fc39011c52eaac0ee03.png)


# *TFG_AD_risk_factors*: Uncovering Environmental Risk Factors through GWAS


This project identifies environmental risk factors for Alzheimer’s Disease (AD) by integrating GWAS summary statistics from multiple traits. We merge and standardize these datasets, apply LD‑pruning and multiple p‑value thresholds, then quantify shared genetic signals via multiscale statistical tests (Hypergeometric, Fisher’s Exact, Monte Carlo) and a custom Tukey‑trimean statistic.

## Tech Stack

![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54) 



## Workflow and Code Overview

### 1. Preprocessing GWAS Data

- **Clean & Harmonize** raw GWAS files into uniform tables with columns:  
  `rs_id`, `chr`, `pos`, `effect_allele`, `alt_allele`, `beta`, `SE`, `p_value`.
- **LD‑Prune** to keep SNPs ≥ 200 kb apart.
- **Threshold Filtering** at user‑defined p‑value cutoffs (1e‑5, 1e‑6, 1e‑8).
- **Visualization**: Manhattan & QQ plots per trait.

**Scripts:**

- `preprocessing_data.R`  
  - _Inputs_: raw file path, trait code, separator, p‑value threshold  
  - _Outputs_:  
    - `preprocessed_data/{trait}_snp_clean.tsv`  
    - `preprocessed_data/{trait}_effect_snp_clean.tsv`  
    - `preprocessed_data/{trait}_snp_significant.tsv`  
    - `preprocessed_data/{trait}_effect_snp_significant.tsv`  
    - `Images/{trait}_manhattan_plot.png`  
    - `Images/{trait}_qq_plot.png`

- `filter_and_plot_significant_snps.R`  
  - _Inputs_: p‑value threshold, `*_effect_snp_clean.tsv`  
  - _Outputs_: updated significant files + Manhattan plots  

- `liftover_script.py`  
  - _Inputs_: UCSC chain file + raw trait text  
  - _Outputs_: `Raw_Data/{trait}.tsv.gz` (GRCh38)

### 2. Significance Assessment

- Compare trait‑specific SNPs’ p‑values in AD across thresholds.
- Compute **Tukey trimean** of –log₁₀(p) and plot horizontal barplots.

**Script:** `overall_significance.R`  
- _Inputs_: threshold, trait codes, `*_effect_snp_clean.tsv`  
- _Outputs_: `processed_data/sig_overAD_{threshold}.png`

### 3. Data Processing
#### 3.1 Overlap Matrix (Non‑scaled)

- Merge significant SNPs for each trait pair, enforce LD‑pruning, count overlaps.
- Validate with Hypergeometric, Fisher’s Exact, and Monte Carlo.

**Script:** `pairwise_trait_processing.R`  
- _Inputs_: preprocessed significant/clean files, `distance_threshold`, `p_value_threshold`, skip list  
- _Outputs_:  
  - `processed_data/significant_snp_matrix.tsv`  
  - `processed_data/pairwise_validation_control.tsv`  
  - `processed_data/{trait1}_{trait2}_filtered_snps.tsv`

#### 3.2 Scaled Matrix

- Divide each overlap count by total significant SNPs in trait₁.

**Script:** `fix_and_rescale_matrix.py`  
- _Inputs_: `significant_snp_matrix.tsv`  
- _Outputs_: `processed_data/rescaled_significant_snp_matrix.tsv`

### 4. Exploratory Analysis
- **`checking_proc.R`**  
  - Prune overlapping SNPs (≥ 200 kb), compute trimean of AD p‑values, Monte Carlo validation, and plot distributions.  
  - _Outputs_: trimean tables + distribution plots  

- **`checking_proc_part2.R`**  
  - Build asymmetric trimean matrix, run NMDS, visualize trait similarity in AD p‑value patterns.  
  - _Outputs_: NMDS coord files + projection plots  

### 5. Validation 
- **`validating_data.R`**: core functions for N, K, n, k calculations, enrichment tests, Monte Carlo histogram.  
- **`validation_per_trait.R`**: wrapper for single‑trait validation (long/short formats).

_Outputs_: console summaries + `Images/montecarlo_{trait}.png`

### 6. Dimensionality Reduction & Clustering
- **`nmds_and_clustering.R`**: reads (re)scaled matrix, runs NMDS, then Mclust & DBSCAN.  
- **Modular scripts**: `NMDS_script.R` (distance prep + NMDS) & `clustering.R` (cluster assignments & plots).

_Outputs_: NMDS coords, projection images, cluster maps.

### 7. Custom Statistic
- **`compute_trait_concordance_matrix.R`**: computes signed concordance (logSign) matrices.  
- **`directionalconditioning.R`**: visualizes Tukey Mean Z² conditioning.  
- **`1e-7visualization.R` / `1e-8visualization.R`**: threshold‑specific NMDS & clustering.

_Outputs_: concordance matrices + NMDS/clustering visualizations.


## Overall Project structure
```graphql
AD_risk_factors/
├── Raw_Data/             #Placeholder for raw GWAS summary statistics (download link provided)
├── preprocessed_data/    #Contains preprocessed trait files (download link provided)
├── processed_data/       #Contains processed trait data (download link provided), Overlap matrices, validation outputs
├── Images/               #Visualizations (Manhattan and QQ plots, NMDS, clusters)
├── testing_significance/ #Barplots and scripts for threshold comparison
└── overall_significance.R   #Script for assessing significance of AD over other traits
├── preprocessing_data.R  #R script for data preprocessing
├── filter_and_plot_significant_snps.R   #R script for data preprocessing but all pairwise combinations
├── liftover_script.py
├── pairwise_trait_processing.R    #R script for processing and validation of all pairwise combinations
├── fix_and_rescale_matrix.py #Python script reescaling the matrix
├── checking_proc.R
├── checking_proc_part2.R
├── validating_data.R
├── validation_per_trait.R
├── nmds_and_clustering.R #Generates NMDS projections and applies Mclust/DBSCAN clustering
├── NMDS_script.R
├── clustering.R
├── compute_trait_concordance_matrix.R
├── directionalconditioning.R
├── 1e-7visualization.R
├── 1e-8visualization.R
├── .gitattributes        #Git LFS configuration (not configurated)
├── .gitignore            #Git ignore configuration
├── LICENSE               #MIT License
└── README.md             #Project overview and workflow 
```

## Directory structure and data localization

In order for the scripts to work, the following directories are needed:

| **Directory**          | **Intended content of the directory** | **Comments** |
|------------------------|------------|-------------|
| **Raw_Data/**          | GWAS summary statistics files for each of the investigated traits. | The full dataset is too large for GitHub. Thus, it contains a placeholder for GWAS summary statistics files. **[Download Here](https://drive.google.com/drive/folders/1K9kSSyzRsOEMwewDZqyXIZHwqqEOPC9P?usp=drive_link)**. |
| **preprocessed_data/** | Contains four files per trait generated by the `preprocessing_data.R` script: <br> - `trait_effect_snp_clean.tsv` <br> - `trait_effect_snp_significant.tsv` <br> - `trait_snp_clean.tsv` <br> - `trait_snp_significant.tsv` <br><br> **Effect SNP files** (`*_effect_snp_*`): Include columns `rs_id`, `chr`, `pos`, `SE`, `beta`, and `p_value`. <br> **SNP files** (`*_snp_*`): Include columns `rs_id`, `chr`, `pos`, `effect_allele`, and `alt_allele`. | Placeholder for preprocessed trait files.  **[Download Here](https://drive.google.com/drive/folders/118cy-vei9lbwYPbTSJ0FD5uAj3ukUy68?usp=drive_link)**. |
| **processed_data/**   | Contains variants that overlap AD-trait, that are merged and filtered, one file per trait, generated by the `processing_data` script.  | Placeholder for processed trait data. **[Download Here](https://drive.google.com/drive/folders/17qK7DhVPpLyZ83b2vtuxb-o5ip26JMMA?usp=drive_link)**.|
| **Images/**           | Contains visualization plots of the different scripts used in the project.  | This directory is a placeholder for those figures. **[Download Here](https://drive.google.com/drive/folders/1Rj-Wmyp2CZI_iRq-nX3Ykjcv3vhR8nXM?usp=drive_link)**|
| **testing_significance/**           | Contains the script and the plots for the assessment of the p-value thresholds. | In this project, these are the values analysed as potential thresholds: 1e-5,1e-6,1e-8.|

## Run Locally -Quickstart-

1. Clone the repository

```bash
  git clone https://github.com/MLM19/AD_risk_factors.git
  cd AD_risk_factors
```

2. Download the Data
Download the GWAS summary statistics on the traits that are going to be rearched. Place those files into the Raw_Data/ directory. 
The data used in this reasearch can be found in the placeholder as a google drive link. 

3. Install dependencies
- **R:** The main analysis scripts are written in R.
- **Dependencies:** Install required of the following R packages, in order to run the scripts:
  ```R
  install.packages(c(
  "data.table","dplyr","ggplot2","qqman",
  "smacof","mclust","dbscan","R.utils",
  "ggrepel", "tibble"
  ))
  
  ```
- **Python:** Some analysis and transformations are performed in Python.
- **Dependencies:** Install and load the following libraries for Python in order to run the scripts:
  ```python
  pip install pandas numpy matplotlib
  
  ```

4. Run the scripts
- Preprocessing
  (If build conversion is needed)
  ```bash
  python liftover_script.py
  ```
  ```R
  source("preprocessing_data.R")
  source("filter_and_plot_significant_snps.R")
  ```

- Assess significance
  ```R
  source("overall_significance.R")
  ```
- Build overlap matrix & validation
  ```R
  source("pairwise_trait_processing.R")
  ```
  ```bash
  python fix_and_rescale_matrix.py
  ```
- Exploratory NMDS and clustering
  ```R
  source("checking_proc.R")
  source("checking_proc_part2.R")
  source("nmds_and_clustering.R")
  ```
- Compute custom stats and visualize it
  ```R
  source("compute_trait_concordance_matrix.R")
  source("directionalconditioning.R")
  source("1e-7visualization.R")
  source("1e-8visualization.R")
  ```


## Author

- [@MLM19](https://www.github.com/MLM19)

## Notes
- GWAS data uses **GRCh38** as the reference genome.
- Any compatible GWAS summary statistics can be used with this project, provided the format matches the configurations or is added manually into the configuration section of the `preprocessing_data.R`.

  
## License

[MIT](https://choosealicense.com/licenses/mit/)


