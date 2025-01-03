<p align="center">
  <img width="300" src="vignettes/SignatureSets_Logo.jpg">
</p>


# PredictioR: An R Package for Biomarker Discovery in Immuno-Oncology Therapy Response

## Introduction
    
The `PredictioR` is an R package to perform comprehensive analyses for biomarker discovery in Immuno-Oncology (IO) response. It supports pan-cancer, cancer-specific, and treatment-specific settings. Additionally, it allows assessment of the importance of key clinical variables such as age and sex. The package includes multiple algorithms for predicting IO response and provides implementations for computing IO signature scores.

The curated IO data and signature resources can be accessed and downloaded from the following locations:

- IO Data: Clinical and molecular profiles are publicly available at [ORCESTRA](https://www.orcestra.ca/clinical_icb).

- IO Signatures: Curated IO signatures can be downloaded from [SignatureSets GitHub repository](https://github.com/bhklab/SignatureSets)."


## Setup
                                                                 
The latest version of PredictioR repository can be found on the [PredictioR GitHub repository](https://github.com/bhklab/PredictioR). The package is not yet on CRAN or Bioconductor. 

It is essential that you have R 4.4.1 or above already installed on your computer or server. PredictioR utilizes many other R packages that are currently available from CRAN, Bioconductor and GitHub. Before installing PredictioR, please install all dependencies by executing the following command in R console:

The dependencies includes MultiAssayExperiment, survival, survcomp, GSVA, meta, ggplot2 and ggrepel.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

depens<-c( 'MultiAssayExperiment', 'survival', 'survcomp', 'GSVA', 'meta', 'ggplot2', 'ggrepel')
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))  BiocManager::install(depen,update = FALSE)
}

```

You can install it from GitHub:

``` bash

devtools::install_github("bhklab/PredictioR")
library(PredictioR) 

```

To set up the repository, please download this folder locally:

``` bash

git clone https://github.com/bhklab/PredictioR
cd PredictioR

```

## Citation 
                                                                  
If the data or computational functions from the PredictioR package are used in your publication, please cite the following paper(s):                                                                  
- [Bareche, Y., Kelly, D., Abbas-Aghababazadeh, F. et al., Annals of Oncology 2022](https://pubmed.ncbi.nlm.nih.gov/36055464/).
                                                                      
- [Mammoliti, Anthony, et al., Nature communications, 2021](https://pubmed.ncbi.nlm.nih.gov/34608132/)


