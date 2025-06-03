
![Logo](https://eu01.edcwb.com/buscador/img/centros/logogrande/54808-a088d273bd364fc39011c52eaac0ee03.png)


# *TFG_AD_risk_factors*: Uncovering Environmental Risk Factors through GWAS


This project aims to identify environmental risk factors for Alzheimer’s Disease (AD) using Genome-Wide Association Study (GWAS) data. By merging and filtering GWAS summary statistics from multiple traits, we investigate genetic overlaps between AD and environmental risk factors. The analysis is performed using R and Python, and various statistical tests (e.g., Hypergeometric, Fisher's Exact, Monte Carlo simulations) to validate the results.


## Tech Stack

![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54) 

## Overview

### 1. Preprocessing GWAS data
- **Cleaning data**: Ensure all the data contains the necessary columns, convert data types, filter out noise, quality check, and visualize.
- **Configuration dictionary**: Manually add any configuration (in a different format for a database that is not already in the dictionary) to the dictionary to allow configuration. Must contain names of the columns and the name to which it should mutate, so that all the databases contain the same columns. 
- **Filtering by threshold**: Once the clean data is obtained in this step, filtering by threshold can be done with two scripts: `preprocessing_data.R` and `filter_and_plot_significant_snps.R`
- **Plot variants distribution by trait**: To visualize the data in a Manhattan plot, it must be filtered by significance. 

### 2. Deciding on a significant threshold
- **Dynamically modifiable threshold**: Compared 1e-5, 1e-6, or 1e-8 significant threshold for the p values on the variants for each trait. (Any threshold can be specified, but those are the ones I checked)
- **Significant variants of a selected phenotype over AD**: Select the significant variants of a phenotype of our choice, find the level of significance of those variants in Alzheimer's disease, we perform a trimean.
- **Plotting the significance**: Create a plot to assess which traits have more significant variants for the selected threshold. 

### 3. Data Processing
- **Merging Significant SNPs**: Identifying significant SNPs in one trait and checking their presence in another trait.
- **Filtering by Distance**: Ensuring selected SNPs are at least 200kb apart to minimize LD bias.
- **Significance Testing**: Evaluating SNP significance based on a certain threshold for p-values.
- **Creates matrix and control document**: Create a matrix with significant, filtered, merged hits of each pairwise combination and a control file where each step is registered; merged size, filtered size, significant hits, and then the validation results.
- **Perform validation tests**: Perform 3 statistical tests to ensure that the obtained results are not given by random chance; Hypergeometric test, Fisher's exact test, and Montecarlo simulation.  

### 4. Scripts
- **`preprocessing_data.R`**: Manually modify configuration database and add those formats that are not included. Clean data by filtering and removing NA's, filtering noise out. Generates Manhattan plots and QQ plots. 
- **`significanceplt.py`**: Manually provide: p-value threshold and list of trait's codes. Generates an image for each threshold in which the trait that contains the variants with most significance in AD are on top. Statistic for the analysis is trimean. 
- **`filter_and_plot_significant_snps.R`**: Automated script to SNP merging and filtering, producing processed GWAS datasets, generates plots and summary statistics for GWAS comparisons. The threshold is modifiable, and traits can be skipped if re-running the script. Creates Manhattan plots to see the distribution of variants in each trait.
- **`pairwise_trait_processing.R`**: Automated script for processing each pairwise combination of traits and evaluating the number of filtered merged significant variants. Threshold and distance are modifiable, and traits can be skipped when re-running the script. Creates a matrix of all the significant hits (`processed_data/significant_snp_matrix.tsv`) and a control document that saves the validation results(`processed_data/pairwise_validation_control.tsv`). 
- **`fix_and_rescale_matrix.py`**: Standardize the matrix by dividing the number of significant hits found in the merged of two traits by the number of significant hits of trait1 (the number of significant SNPs found in trait1 when filtering for the threshold). In this script, I manually created the dictionary with the code name of the trait and the value. Generates another file `rescaled_significant_snp_matrix.tsv`
- **`nmds_and_clustering.R`**: Once the matrix is re-scaled, and the distance approximation is done, then this produces the Non-metric multidimensional scaling (NMDS), which is a projection of the data that allows us to see to how similar the traits are amongst themselves. Then we perform two types of clustering: mClust and DBSCAN. 

### 3. Results and Outputs
- **Merged SNP Lists**: Filtered SNPs for each trait pair. Can be found in the `processed_data/` directory. 
- **Statistical Reports**: Summary tables detailing the processing results at each step. 
- **Visual Representations**: Significance assessment plot, Manhattan plots, monte carlo simulation plots, NMDS projections, and clustering plots. 

## Overall Project structure
```graphql
AD_risk_factors/
├── Images/               #Visualizations (e.g., Manhattan and QQ plots)
├── Raw_Data/             #Placeholder for raw GWAS summary statistics (download link provided)
├── preprocessed_data/    #Contains preprocessed trait files (download link provided)
├── processed_data/       #Contains processed trait data (download link provided)
├── testing_significance/ #Contains significanceplt.py and 3 plots based on significance analysis based on trimeans of the data
└── significanceplt.py    #Python script for assessing significance of AD over other traits
├── .gitattributes        #Git LFS configuration (not configurated)
├── .gitignore            #Git ignore configuration
├── LICENSE               #MIT License
├── README.md             #This README file 
├── preprocessing_data.R  #R script for data preprocessing 
├── filter_and_plot_significant_snps.R   #R script for data preprocessing but all pairwise combinations
├── pairwise_trait_processing.R    #R script for processing and validation of all pairwise combinations
├── fix_and_rescale_matrix.py #Python script reescaling the matrix
└── nmds_and_clustering.R #R script for visualizing the data and interpreting
```

## Workflow
1. Data Acquisition:
Download the necessary GWAS summary statistics files from the provided Google Drive links and place them in the appropriate directories.
Or else, find other raw data to perform these tests. Most recommended places: GWAS Catalog (https://www.ebi.ac.uk/gwas/downloads/summary-statistics) and UK biobank (https://www.ukbiobank.ac.uk/).

2. Preprocessing:
Run the `preprocessing_data.R` script to clean and standardize the GWAS data. This script will pre-process each trait individually.
If you already have `*_effect_snp_clean.tsv` and want to assess for another threshold, you can use `filter_and_plot_significant_snps.R` to preprocess all the traits that can be found in the preprocessed file.
Both scripts generate visual outputs (Manhattan plots) to help interpret the results. It also produces `*_snp_significant.tsv`. Thresholds can be dynamically changed.

- 2.1. Assessing significance: 
  Run the `significanceplt.py` script to assess a defined threshold for p-value variants. 
  Must define the threshold that you want to assess and the traits that are going to be used. Must have   `*_effect_snp_clean.tsv` and `*_snp_clean.tsv` to run the script.

3. Processing:
Execute the `pairwise_trait_processing.R` script to merge significant SNPs and apply distance filtering. Then find the significant, merged, and filtered variants for each pairwise combination to create an asymmetric matrix. It also computes the validations. 
Must have `*_snp_significant.tsv` and `*_snp_clean.tsv` in the `preprocessing/` directory to run the script

- 3.1. Re-escaling matrix:
  Run the `fix_and_rescale_matrix.py` script to rescale the matrix before the visualization. 

4. Visualization:
Run the `nmds_and_clustering.R` script to perform the NMDS projection and the clustering algorithms on the data. 

## Directory structure and data localization

In order for the scripts to work, the following directories are needed:

| **Directory**          | **Intended content of the directory** | **Comments** |
|------------------------|------------|-------------|
| **Raw_Data/**          | GWAS summary statistics files for each of the investigated traits. | The full dataset is too large for GitHub. Thus, it contains a placeholder for GWAS summary statistics files. **[Download Here](https://drive.google.com/drive/folders/1K9kSSyzRsOEMwewDZqyXIZHwqqEOPC9P?usp=drive_link)**. |
| **preprocessed_data/** | Contains four files per trait generated by the `preprocessing_data.R` script: <br> - `trait_effect_snp_clean.tsv` <br> - `trait_effect_snp_significant.tsv` <br> - `trait_snp_clean.tsv` <br> - `trait_snp_significant.tsv` <br><br> **Effect SNP files** (`*_effect_snp_*`): Include columns `rs_id`, `chr`, `pos`, `SE`, `beta`, and `p_value`. <br> **SNP files** (`*_snp_*`): Include columns `rs_id`, `chr`, `pos`, `effect_allele`, and `alt_allele`. | Placeholder for preprocessed trait files.  **[Download Here](https://drive.google.com/drive/folders/118cy-vei9lbwYPbTSJ0FD5uAj3ukUy68?usp=drive_link)**. |
| **processed_data/**   | Contains variants that overlap AD-trait, that are merged and filtered, one file per trait, generated by the `processing_data` script.  | Placeholder for processed trait data. **[Download Here](https://drive.google.com/drive/folders/17qK7DhVPpLyZ83b2vtuxb-o5ip26JMMA?usp=drive_link)**.|
| **Images/**           | Contains visualization plots of variant distributions for each trait. Includes one **Manhattan plot** per trait.  | This directory is a placeholder for those figures. **[Download Here](https://drive.google.com/drive/folders/1Rj-Wmyp2CZI_iRq-nX3Ykjcv3vhR8nXM?usp=drive_link)**|
| **testing_significance/**           | Contains the script and the plots for the assessment of the p-value thresholds. | In this project, these are the values analysed as potential thresholds: 1e-5,1e-6,1e-8.|

## Prerequisites

- **R:** The main analysis scripts are written in R.
- **Dependencies:** Install required of the following R packages, in order to run the scripts:
  ```r
  #validating.R
  install.packages("data.table")
  install.packages("dplyr")
  install.packages("ggplot2")
  install.packages("qqman")
  install.packages("R.utils")

  #filter_and_plot_significant_snps.R
  install.packages("dplyr")
  install.packages("data.table")
  install.packages("qqman")

  #nmds_and_clustering.R
  install.packages("smacof")
  install.packages("ggplot2")
  install.packages("ggrepel")
  install.packages("mclust")
  install.packages("dbscan")
  install.packages("readr")

  #pairwise_trait_processing.R
  install.packages("dplyr")
  install.packages("tibble")
  install.packages("data.table")
  ```
  
- **Python:** Some analysis and transformations are performed in Python.
- **Dependencies:** Install and load the following libraries for Python in order to run the scripts:
  ```python
  #significanceplt.py
  pip install pandas numpy matplotlib

  #fix_and_rescale_matrix.py
  pip install pandas
  
  ```


## Run Locally

1. Clone the project

```bash
  git clone https://github.com/MLM19/AD_risk_factors.git
```

2. Go to the project directory

```bash
  cd AD_risk_factors
```

3. Download the Data
Download the GWAS summary statistics on the traits that are going to be rearched. Place those files into the Raw_Data/ directory. 
The data used in this reasearch can be found in the placeholder as a google drive link. 

4. Install dependencies
(See the prerequisites section) For R scripts, install necessary packages.

5. Run the scripts
From the main directory of the project run the scripts.
The script are prepared to work with the given directories however they can be adjusted by changing the paths of the data in the scripts.
To run the R scripts, you can write `source(name_of_the_script.R)`.

## Authors

- [@MLM19](https://www.github.com/MLM19)

## Notes
- GWAS data uses **GRCh38** as the reference genome.
- Any compatible GWAS summary statistics can be used with this project, provided the format matches the configurations or is added manually into the configuration section of the `preprocessing_data.R`.
## License

[MIT](https://choosealicense.com/licenses/mit/)


