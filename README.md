
![Logo](https://eu01.edcwb.com/buscador/img/centros/logogrande/54808-a088d273bd364fc39011c52eaac0ee03.png)


# *TFG_AD_risk_factors*: Uncovering Environmental Risk Factors through GWAS


This project aims to identify environmental risk factors for Alzheimer’s Disease (AD) using Genome-Wide Association Study (GWAS) data. By merging and filtering GWAS summary statistics from multiple traits, we investigate genetic overlaps between AD and environmental risk factors. The analysis is performed using R and Python, and various statistical tests (e.g., Hypergeometric, Fisher's Exact, Monte Carlo simulations) to validate the results.


## Tech Stack

![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54) 

## Overview

### 1. Preprocessing GWAS data
- **Cleaning data**:
- **Plot variants distribution by trait**:

### 2. Deciding singnificant threshold
- **Dynamically modifiable threshold**: Compared 1e-5, 1e-6 or 1e-8 signifcant threshold for the p values on the variants for each trait. (Any threshold can be specified but those are the ones I checked)
- **Significant variants of a selected phenotype over AD**: Select the significant variants of a phenotype of our choice, find the level of significance of those variants in the Alzheimer's disease, we perform a trimean.
- **Plotting the significance**: Create a plot to assess which traits have more significant variants for the selected threshold. 

### 3. Data Processing
- **Merging Significant SNPs**: Identifying significant SNPs in one trait and checking their presence in another trait.
- **Filtering by Distance**: Ensuring selected SNPs are at least 200kb apart to minimize LD bias.
- **Significance Testing**: Evaluating SNP significance based on p-values.

### 4. Main Scripts
- **`filter_and_plot_significant_snps.R`**: Automatized script to SNP merging and filtering, producing processed GWAS datasets, Generates plots and summary statistics for GWAS comparisons. The threshold is modifiable and traits can be skipped if re-running the script. Creates Manhattan plots to see the distribution of variants in each trait. 



### 3. Results and Outputs
- **Merged SNP Lists**: Filtered SNPs for each trait pair.
- **Statistical Reports**: Summary tables detailing significant associations.
- **Visual Representations**: Manhattan plots, bar charts, and networks illustrating SNP associations.

## Overall Project structure
```graphql
AD_risk_factors/
├── Images/               #Visualizations (e.g., Manhattan and QQ plots)
├── Raw_Data/             #Placeholder for raw GWAS summary statistics (download link provided)
├── preprocessed_data/    #Contains preprocessed trait files (download link provided)
├── processed_data/       #Contains processed trait data (download link provided)
├── testing_significance/ #Contains significanceplt.py and 3 plots based on significance analysis based on trimeans of the data
├── .gitattributes        #Git LFS configuration (not configurated)
├── .gitignore            #Git ignore configuration
├── LICENSE               #MIT License
├── README.md             #This README file 
├── only_visualizing.R    #R script for visualization 
├── preprocessing_data.R  #R script for data preprocessing 
├── processing_data.R     #R script for data processing
├── validating_data.R     #R script for validation
├── process_large_merged.R #R script for data processing of pairwise combination of traits with large number of merged SNPs
└── validation_per_trait.R #R script for validating one trait at the time
```

## Workflow
1. Data Acquisition:
Download the necessary GWAS summary statistics files from the provided Google Drive links and place them in the appropriate directories.

2. Preprocessing:
Run the `preprocessing_data.R` script to clean and standardize the GWAS data.
The only_visualizing.R script generates visual outputs (e.g., Manhattan and QQ plots) to help interpret the results, but those plots can also be generated as part of the `preprocessing_data.R` script.

2.1. Assessing significance: 
Run the `significanceplt.py` script to assess a defined threshold for p value variants. 
Must define the threshold that you want to assess, the traits that are going to be used.

4. Processing:
Execute the `processing_data.R` script to merge significant SNPs and apply distance filtering.

5. Validation:
Use the `validating_data.R` script to construct contingency tables and perform statistical tests, ensuring the observed associations are not due to random chance.

## Directory structure and data localization

In order for the scripts to work, the following directories are needed:

| **Directory**          | **Intended content of the directory** | **Comments** |
|------------------------|------------|-------------|
| **Raw_Data/**          | GWAS summary statistics files for each of the investigated traits. | The full dataset is too large for GitHub. Thus, it contains a placeholder for GWAS summary statistics files. **[Download Here](https://drive.google.com/drive/folders/1K9kSSyzRsOEMwewDZqyXIZHwqqEOPC9P?usp=drive_link)**. |
| **preprocessed_data/** | Contains four files per trait generated by the `preprocessing_data.R` script: <br> - `trait_effect_snp_clean.tsv` <br> - `trait_effect_snp_significant.tsv` <br> - `trait_snp_clean.tsv` <br> - `trait_snp_significant.tsv` <br><br> **Effect SNP files** (`*_effect_snp_*`): Include columns `rs_id`, `chr`, `pos`, `SE`, `beta`, and `p_value`. <br> **SNP files** (`*_snp_*`): Include columns `rs_id`, `chr`, `pos`, `effect_allele`, and `alt_allele`. | Placeholder for preprocessed trait files.  **[Download Here](https://drive.google.com/drive/folders/118cy-vei9lbwYPbTSJ0FD5uAj3ukUy68?usp=drive_link)**. |
| **processed_data/**   | Contains variants that overlap AD-trait, that are merged and filtered, one file per trait, generated by the `processing_data` script.  | Placeholder for processed trait data. **[Download Here](https://drive.google.com/drive/folders/17qK7DhVPpLyZ83b2vtuxb-o5ip26JMMA?usp=drive_link)**.|
| **Images/**           | Contains visualization plots of variant distributions for each trait. Includes one **Manhattan plot** per trait.  | This directory is a placeholder for those figures. **[Download Here](https://drive.google.com/drive/folders/1Rj-Wmyp2CZI_iRq-nX3Ykjcv3vhR8nXM?usp=drive_link)**|
| **testing_significance/**           | Contains the script and the plots for the assesing of the p value threshold | In this project, these are the values analysed as potential threshold: 1e-5,1e-6,1e-8|

## Prerequisites

- **R:** The main analysis scripts is written in R.
- **Dependencies:** Install required of the following R packages, in order to run the scripts:
  ```r
  #filter_and_plot_significant_snps.R
  install.packages("dplyr")
  install.packages("data.table")
  install.packages("qqman")

  #only_visualizing.R
  install.packages("data.table")
  install.packages("qqman")

  #processing_data.R
  install.packages("dplyr")
  install.packages("tibble")
  install.packages("data.table")

  #validating.R
  install.packages("dplyr")
  install.packages("ggplot2")
  ```
- **Python:** Some analysis and transformations are performed in Python.
- **Dependencies:** Install and load the following libraries for Python in order to run the scripts:
  ```python
  #significanceplt.py
  pip install pandas numpy matplotlib
  
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


