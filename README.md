<p align="center">
  <img width="300" src="vignettes/SignatureSets_Logo.jpg">
</p>


# PredictioR: An R Package for Biomarker Discovery in Immuno-Oncology Therapy Response

## 📖 Introduction
    
`PredictioR` is an R package designed for comprehensive biomarker discovery in Immuno-Oncology (IO) therapy response. It supports pan-cancer, cancer-specific, and treatment-specific analyses. The package also enables assessment of key clinical variables such as age and sex. Multiple algorithms are included for IO response prediction, along with methods for computing IO signature scores.

**Resources:**
- 🧬 **IO Data**: Clinical and molecular profiles – available at [ORCESTRA](https://www.orcestra.ca/clinical_icb)  
- 🧾 **IO Signatures**: Curated IO gene signatures – available from [SignatureSets GitHub repository](https://github.com/bhklab/SignatureSets)


## 🔧 Setup
                                                                 
The latest version of `PredictioR` is available at the [PredictioR GitHub repository](https://github.com/bhklab/PredictioR). The package is not yet on CRAN or Bioconductor.

✅ Requirements
- R version 4.4.1 or higher

📦 Install Dependencies

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

dependencies <- c('MultiAssayExperiment', 'survival', 'survcomp', 'GSVA', 'meta', 'ggplot2', 'ggrepel')
for (pkg in dependencies) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE)
}

```

📥 Install PredictioR from GitHub

``` bash

devtools::install_github("bhklab/PredictioR")
library(PredictioR) 

```

💻 Clone the Repository (optional)

``` bash

git clone https://github.com/bhklab/PredictioR
cd PredictioR

```

More details about function usage and computational methods are provided in the package documentation and [vignettes](https://github.com/bhklab/PredictioR/blob/main/vignettes/PredictioR.Rmd), or via the web application at [predictio.ca](https://predictio.ca/).

## 📁 Repository Structure

```plaintext

PredictioR/
├── 📁 R/            – Core package functions
├── 📁 data/         – Selected and curated IO signatures and datasets
├── 📁 man/          – Function documentation (.Rd files)
├── 📁 vignettes/    – Workflows and usage examples
├── 📄 DESCRIPTION   – Package metadata
└── 📄 README.md     – Overview and setup instructions

```

## 📝 Citation
                                                                  
If you use PredictioR or its datasets in your work, please cite the following papers:                                                                  
- [Bareche, Y., Kelly, D., Abbas-Aghababazadeh, F. et al., Annals of Oncology 2022](https://pubmed.ncbi.nlm.nih.gov/36055464/).
                                                                      
- [Mammoliti, Anthony, et al., Nature communications, 2021](https://pubmed.ncbi.nlm.nih.gov/34608132/)

## Version Information

📌 This is **version: v1.0**, corresponding to the release used in the paper. 
