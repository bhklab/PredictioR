## load libraries

library(MultiAssayExperiment)
library(qs)
library(data.table)

############################################################
## create directory
############################################################
## create a list with name of cohorts as discovery and validation cohorts:

datasets.dir <- "C:/PredictIO_CodeOcean_gene/PredictIO_package/results/datasets"

## if datasets directory exist, remove
if( dir.exists( datasets.dir ) ) {
  unlink( datasets.dir , recursive = TRUE )
  }

dir.create( datasets.dir )

################################################################################
## create expression data
################################################################################

createSummarizedExperiments <- function(dat_file){
  
  dat_icb <- readRDS(dat_file)
  
  ## slot names for expression data
  dat_type <- c("expr", "expr_gene_tpm")
  
  dat_type_ix <- which(dat_type %in% names(dat_icb))
  expr <- assay(dat_icb[[dat_type[dat_type_ix]]])
  clin <- as.data.frame(colData(dat_icb[[dat_type[dat_type_ix]]]))
  annot <- as.data.frame(rowData(dat_icb[[dat_type[dat_type_ix]]])) 

  ## limit to protein-coding genes and remove duplicated genes using gene names
  annot <- annot[annot$gene_type == "protein_coding", , drop=FALSE]
  
  ## remove ENSG*Par_Y
  remove_Par_Y <- grep("PAR_Y",rownames(annot))
  if(length(remove_Par_Y) > 0){
    annot <- annot[-remove_Par_Y, ]
  }
  
  ## keep the smallest ENSG for each gene symbol
  annot <- annot[order(rownames(annot)), , drop=FALSE]
  annot <- annot[!duplicated(annot$gene_name), , drop=FALSE]
  
  ## subset gene expression data
  expr <- expr[rownames(annot), , drop=FALSE]
  ## change feature names from ENSG ids to gene symbols  
  rownames(expr) <- rownames(annot) <- annot$gene_name

  ## create SummarizedExperiments
  eset <- SummarizedExperiment(assay= list("gene_expression"=expr), 
                               rowData=annot, 
                               colData=clin)

  return(eset)
  
}

################################################################################
## load all data: get expression data and prepare the SE objects
################################################################################

dir <- "C:/PredictIO_CodeOcean_gene/PredictIO_package/data/data_ICB"

## list of datasets
datasets <- paste(dir, list.files(dir), sep="/")

## extract study names
studies <- substr(list.files(dir), 5, nchar(list.files(dir))-4)

dat <- lapply(1:length(datasets), function(k) {
  cat("Processing data:", studies[k], "\n")
  createSummarizedExperiments(datasets[k])
  
})

names(dat) <- studies

qsave(dat, file=file.path(datasets.dir, "ICB_data.qs"))

#####################################################
## create data as qs R file (each cohort + cancer type)

dir <- "C:/PredictIO_CodeOcean_gene/PredictIO_package/results/datasets/qs/"

dat <- qread(file.path(datasets.dir, "ICB_data.qs"))
study <- names(dat)

for(k in 1:length(study)){
  
  cancer_type <- unique(colData(dat[[k]])$cancer_type)
  if( length(cancer_type) == 1){
    
    dat_icb <- dat[[k]]
    qsave(dat_icb, file=paste(dir, paste("ICB_", 
                                         paste(paste(study[k], cancer_type, sep="__"), sep=""), 
                                         ".qs", sep=""), sep="/")) 
    
  }
  
  if( length(cancer_type) > 1 ){
    
    for(i in 1:length(cancer_type)){
      
      dat_icb <- dat[[k]]
      sub_clin <- colData(dat_icb)[colData(dat_icb)$cancer_type == cancer_type[i], ]
      dat_icb <- dat_icb[, rownames(sub_clin)]
      qsave(dat_icb, file=paste(dir, paste("ICB_", 
                                           paste(paste(study[k], cancer_type[i], sep="__"), sep=""), 
                                           ".qs", sep=""), sep="/")) 
    }
  }
}

########################################################################
## merge all cohorts as a list as a qs file
dir <- "C:/PredictIO_CodeOcean_gene/PredictIO_package/results/datasets/qs"
ICB_data <- lapply(1:length(list.files(dir)), function(i){
  
  print(i)
  qread( paste(dir, list.files(dir)[i], sep="/"))
  
})

names(ICB_data) <- substr(list.files(dir), 5, nchar(list.files(dir)) - 3)

qsave(ICB_data, file= "C:/PredictIO_CodeOcean_gene/PredictIO_package/results/datasets/ICB_data_updated.qs")

#####################################################
## create data as rda file (each cohort + cancer type)

dir <- "C:/PredictIO_CodeOcean_gene/PredictIO_package/results/datasets/rda/"

dat <- qread(file.path(datasets.dir, "ICB_data.qs"))
study <- names(dat)

for(k in 1:length(study)){
  
  cancer_type <- unique(colData(dat[[k]])$cancer_type)
  if( length(cancer_type) == 1){
    
    dat_icb <- dat[[k]]
    save(dat_icb, file=paste(dir, paste("ICB_", 
                                         paste(paste(study[k], cancer_type, sep="__"), sep=""), 
                                         ".rda", sep=""), sep="/")) 
    
  }
  
  if( length(cancer_type) > 1 ){
    
    for(i in 1:length(cancer_type)){
      
      dat_icb <- dat[[k]]
      sub_clin <- colData(dat_icb)[colData(dat_icb)$cancer_type == cancer_type[i], ]
      dat_icb <- dat_icb[, rownames(sub_clin)]
      save(dat_icb, file=paste(dir, paste("ICB_", 
                                           paste(paste(study[k], cancer_type[i], sep="__"), sep=""), 
                                           ".rda", sep=""), sep="/")) 
    }
  }
}

########################################################################
## merge all cohorts as a list rda file
dir <- "C:/PredictIO_CodeOcean_gene/PredictIO_package/results/datasets/rda"
ICB_data <- lapply(1:1, function(i){ # length(list.files(dir))
  
  print(i)
  load( paste(dir, list.files(dir)[i], sep="/"))
  
})

names(ICB_data) <- substr(list.files(dir), 5, nchar(list.files(dir)) - 4)

save(ICB_data, file= "C:/PredictIO_CodeOcean_gene/PredictIO_package/results/datasets/ICB_data.rda")

#######################################################################
## NOTE: no duplicated samples across studies
