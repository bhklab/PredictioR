<p align="center">
  <img width="300" src="vignettes/SignatureSets_Logo.jpg">
</p>


# PredictioR: An R Package for Biomarker Discovery in Immuno-Oncology Therapy Response

## Overview

**PredictioR** is an R package for biomarker discovery and predictive modeling in **immuno-oncology (IO) therapy response**. It supports **pan-cancer**, **cancer-specific**, and **treatment-specific** analyses, and integrates clinical covariates like **age**, **sex**, and **tumor type** into its modeling workflows.

The package includes:
- Signature scoring methods for curated IO gene signatures  
- IO response prediction algorithms  
- Functions for clinical association analysis  
- Built-in support for curated datasets and integration with the [SignatureSets](https://github.com/bhklab/SignatureSets) package

---

**Data Resources**
- **IO Datasets**: Clinical and molecular profiles of IO-treated cohorts, available at: [ORCESTRA](https://www.orcestra.ca/clinical_icb)  
- **IO Signatures**: IO gene signatures available from the companion repository: [SignatureSets GitHub repository](https://github.com/bhklab/SignatureSets)


## Installation
                                                                 
**Dependencies**  
Requires **R 4.4.1 or higher**

Install required packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

dependencies <- c(
  "MultiAssayExperiment", "survival", "survcomp",
  "GSVA", "meta", "ggplot2", "ggrepel"
)

for (pkg in dependencies) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE)
  }
}

```

Install PredictioR from GitHub

``` bash

devtools::install_github("bhklab/PredictioR")
library(PredictioR) 

```

For source-level exploration:

``` bash

git clone https://github.com/bhklab/PredictioR
cd PredictioR

```
---

## Documentation and usage examples:

More details about function usage and computational methods are provided in the package documentation and [vignettes](https://github.com/bhklab/PredictioR/blob/main/vignettes/PredictioR.Rmd), or via the web application at [predictio.ca](https://predictio.ca/).

---

## Repository Structure

```plaintext

PredictioR/
├── 📁 R/            – Core package functions
├── 📁 data/         – Selected and curated IO signatures and datasets
├── 📁 man/          – Function documentation (.Rd files)
├── 📁 vignettes/    – Workflows and usage examples
├── 📄 DESCRIPTION   – Package metadata
└── 📄 README.md     – Overview and setup instructions

```
---

## Citation
                                                                  
If you use PredictioR or its datasets in your work, please cite the following papers:                                                                  
- [Bareche, Y., Kelly, D., Abbas-Aghababazadeh, F. et al., Annals of Oncology 2022](https://pubmed.ncbi.nlm.nih.gov/36055464/)
                                                                      
- [Mammoliti, Anthony, et al., Nature communications, 2021](https://pubmed.ncbi.nlm.nih.gov/34608132/).

