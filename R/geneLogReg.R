#' Fit Logistic Regression Model for Genes: Continuous Expression Variable
#' @description
#' Fits a logistic regression model with continuous expression data for each gene.
#' 
#' @param dat.icb A MultiAssayExperiment (MAE) object, SummarizedExperiment (SE) object, or a data frame or matrix of gene expression data.
#' @param clin If dat.icb is a data frame or matrix, then it contains clinical data (as data frame or matrix). By default, it is NULL.
#' @param time.censor Possible censoring in months.
#' @param missing.perc A cutoff to remove genes with zero expression across samples.
#' @param const.int A pseudocount is added to the TPM (Transcripts Per Million) values before performing a log transformation.
#' @param n.cutoff Minimum number of samples included in the association analysis.
#' @param feature A vector of character strings for selected features.
#' @param study Name of study.
#' @param n0.cutoff Minimum number of non-responders in the analysis.
#' @param n1.cutoff Minimum number of responders in the analysis.
#' @param cancer.type Name of the cancer type for the given study.
#' @param treatment Name of the treatment for the given study. 
#'
#' @return A subset of results using an object of class logistic regression representing the fit. 
#' Outcome: Immunotherapy response outcome i.e., R (responder) and NR (non-responder).
#' Gene: Name of selected genes.
#' Study: Name of study.
#' Coef: Estimate of treatment effect i.e., log odds ratio.
#' SE: Standard error of treatment estimate.
#' N: Number of samples.
#' Pval: Estimated p-value.
#' Cancer_type: A character shows the cancer type.
#' Treatment: A character shows the treatment type.
#' @export
#'
#' @examples
#' Assess the association between selected genes and response (R vs NR) in immunotherapy.
#' geneLogReg(dat.icb = ICB_small_Mariathasan,
#'            missing.perc = 0.5,
#'            const.int = 0.001,
#'            n.cutoff = 15,
#'            feature = c('CXCL9', 'CXCL10', 'TIGIT', 
#'                          'CD83', 'STAT1', 'CXCL11',
#'                          'CXCL13', 'CD8A', 'CTLA4'),
#'            study = 'ICB_Mariathasan', 
#'            n0.cutoff = 10,
#'            n1.cutoff = 10,
#'            cancer.type = 'Bladder',
#'            treatment = 'PD-1/PD-L1')
#' 
geneLogReg <- function(dat.icb, clin = NULL, missing.perc, const.int=0.001, n.cutoff, feature, study,
                       n0.cutoff, n1.cutoff, cancer.type, treatment){
  
  if( !class(dat.icb) %in% c("SummarizedExperiment", "MultiAssayExperiment", "matrix", "data.frame") ){
    
    stop(message("function requires SummarizedExperiment, MultiAssayExperiment, data.frame, or matrix class of data"))
    
  }
  
  if( class(dat.icb) == "MultiAssayExperiment"){
    
    dat <- createSE(dat.icb)
    dat_expr <- assay(dat)
    dat_clin <- colData(dat)
    
  }
  
  if( class(dat.icb) == "SummarizedExperiment"){
    
    dat_expr <- assay(dat.icb)
    dat_clin <- colData(dat.icb)
    
  }
  
  if( sum(nrow(clin)) > 0  ){
    
    dat_expr <- dat.icb
    dat_clin <- clin
    
  }
  
  #cancer_type <- names( table( dat_clin$cancer_type )[ table( dat_clin$cancer_type ) >= n.cutoff ] )
  
  #message(paste(study, cancer_type, sep="/"))
  
  #data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type]
  data <- dat_expr
  remove <- rem(data, missing.perc, const.int)
  
  if( length(remove) ){
    data <- data[-remove,]
  }
  
  data <- as.matrix( data[ rownames(data) %in% feature , ] )
  
  if( nrow(data) & length(dat_clin$response) >= n.cutoff &
      sum(dat_clin$response == "NR", na.rm = TRUE) >= n1.cutoff &
      sum(dat_clin$response == "R", na.rm = TRUE) >= n0.cutoff ){
    
    res <- lapply(1:nrow(data), function(k){
      
      g <- as.numeric( scale( data[k , ] ) )
      names( g ) <- colnames( data )
      
      x <- ifelse( dat_clin$response %in% "R" , 0 ,
                   ifelse( dat_clin$response %in% "NR" , 1 , NA ) )
      
      fit <- glm( x ~ g , family=binomial( link="logit" ) )
      
      data.frame( Outcome = "R vs NR",
                  Gene = rownames(data)[k],
                  Study = study,
                  Coef = round( summary(fit)$coefficients[ "g" , "Estimate"  ] , 3 ),
                  SE = round( summary(fit)$coefficients[ "g" , "Std. Error" ] , 3 ),
                  N = length(x[!is.na(x)]),
                  Pval = summary(fit)$coefficients[ "g" , "Pr(>|z|)" ],
                  Cancer_type = cancer.type,
                  Treatment = treatment)
      
    })
    
    res <- do.call(rbind, res)
    rownames(res) <- NULL
    #  res$FDR <- p.adjust(res$Pval, method = "BH")
    
  }else{  res <- data.frame( Outcome = "R vs NR",
                             Gene = NA,
                             Study = study,
                             Coef = NA,
                             SE = NA,
                             N = NA,
                             Pval = NA,
                             Cancer_type = NA,
                             Treatment = NA)
  
  message("lack of number of samples and/or genes with known immunotherapy survival outcome")
  
  }
  
  return(res)
  
}

