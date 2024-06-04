#' Fit Logistic Regression Model for Gene Signature: Continuous Variable
#' @description
#' Fits a logistic regression model for continuous signature data.
#' 
#' @param dat.icb A MultiAssayExperiment (MAE) object, SummarizedExperiment (SE) object, or a data frame or matrix of gene expression data.
#' @param clin If dat.icb is a data frame or matrix, then it contains clinical data (as data frame or matrix). By default, it is NULL.
#' @param geneSig A numeric vector of computed signature score.
#' @param n.cutoff Minimum number of samples included in the association analysis.
#' @param study Name of study.
#' @param sig.name Name of signature.
#' @param n0.cutoff Minimum number of non-responders in the analysis.
#' @param n1.cutoff Minimum number of responders in the analysis.
#' @param cancer.type Name of the cancer type for the given study.
#' @param treatment Name of the treatment for the given study. 
#'
#' @return A subset of results using an object of class logistic regression representing the fit. 
#' Outcome: Immunotherapy response outcome i.e., R (responder) and NR (non-responder).
#' Gene: Name of selected signature.
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
#' Assess the association of M1 Hwang signature and response. 
#' sig <- geneSigssGSEA(dat.icb = ICB_Liu, 
#'                      sig = M1_Hwang,
#'                      sig.name = 'M1_Hwang',
#'                      missing.perc = 0.5,
#'                      const.int = 0.001,
#'                      n.cutoff = 15,
#'                      sig.perc = 0.8, 
#'                      study = 'ICB_Liu')
#'             
#' geneSigLogReg(dat.icb = ICB_Liu,
#'               geneSig = sig,
#'               n.cutoff = 15,
#'               study =  'ICB_Liu',
#'               sig.name = 'M1_Hwang',
#'               n0.cutoff = 10,
#'               n1.cutoff = 10,
#'               cancer.type = 'Melanoma',
#'               treatment = 'PD-1/PD-L1')
#'                 
geneSigLogReg <- function(dat.icb, clin = NULL, geneSig, n.cutoff, study, sig.name, n0.cutoff, n1.cutoff, cancer.type, treatment){
  
  if( !class(dat.icb) %in% c("SummarizedExperiment", "MultiAssayExperiment", "data.frame", "matrix") ){
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
  
  if( length(dat_clin$response) >= n.cutoff &
      sum(dat_clin$response == "NR", na.rm = TRUE) >= n1.cutoff &
      sum(dat_clin$response == "R", na.rm = TRUE) >= n0.cutoff ){
    
    x <- ifelse( dat_clin$response %in% "R" , 0 ,
                 ifelse( dat_clin$response %in% "NR" , 1 , NA ) )
    
    fit <- glm( x ~ geneSig , family=binomial( link="logit" ) )
    
    res <- data.frame( Outcome = "R vs NR",
                       Gene = sig.name,
                       Study = study,
                       Coef = round( summary(fit)$coefficients[ "geneSig" , "Estimate"  ] , 3 ),
                       SE = round( summary(fit)$coefficients[ "geneSig" , "Std. Error" ] , 3 ),
                       N = length(x[!is.na(x)]),
                       Pval = summary(fit)$coefficients[ "geneSig" , "Pr(>|z|)" ],
                       Cancer_type = cancer.type,
                       Treatment = treatment)
    
    
  }else{  res <- data.frame( Outcome = "R vs NR",
                             Gene = sig.name,
                             Study = study,
                             Coef = NA,
                             SE = NA,
                             N = NA,
                             Pval = NA,
                             Cancer_type = NA,
                             Treatment = NA)
  
  message("lack of number of samples with known immunotherapy response outcome")
  
  }
  
  return(res)
  
}

