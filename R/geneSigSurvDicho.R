#' Fit Proportional Hazards Regression Model for Gene Signature: Dichotomous Variable
#' @description
#' Fits a Cox proportional hazards regression model with dichotomous signature data using the counting process formulation of Andersen and Gill.
#' 
#' @param dat.icb A MultiAssayExperiment (MAE) object, SummarizedExperiment (SE) object, or a data frame or matrix of gene expression data.
#' @param clin If dat.icb is a data frame or matrix, then it contains clinical data (as data frame or matrix). By default, it is NULL.
#' @param geneSig A numeric vector of computed signature score.
#' @param time.censor Possible censoring in months.
#' @param n.cutoff Minimum number of samples included in the association analysis.
#' @param n0.cutoff Minimum number of samples with status 0.
#' @param n1.cutoff Minimum number of samples with status 1.
#' @param study Name of study.
#' @param surv.outcome Overall survival (OS) or progression-free survival (PFS). 
#' @param sig.name Name of signature.
#' @param method The default method to convert a continuous variable into a dichotomous variable is the 'median' method. The first quartile (Q1) and third quartile (Q3) can also be applied. 
#' @param var.type If the variable is dichotomous (by default), then var.type is TRUE.
#' @param cancer.type Name of the cancer type for the given study.
#' @param treatment Name of the treatment for the given study. 
#'
#' @return A subset of results using an object of class 'coxph' representing the fit. 
#' Outcome: Immunotherapy time-to-event outcomes including overall survival (OS) and progression-free survival (PFS).
#' Gene: Name of selected signature.
#' Study: Name of study.
#' Coef: Estimate of treatment effect i.e., log hazard ratio.
#' SE: Standard error of treatment estimate.
#' N: Number of samples.
#' Pval: Estimated p-value.
#' Cancer_type: A character shows the cancer type.
#' Treatment: A character shows the treatment type.
#' 
#' @export
#'
#' @examples
#' Assess the association of EMT Thompson signature and PFS. 
#' sig <- geneSigMean(dat.icb = ICB_small_Liu, 
#'                    sig = EMT_Thompson,
#'                    sig.name = 'EMT_Thompson',
#'                    missing.perc = 0.5,
#'                    const.int = 0.001,
#'                    n.cutoff = 15,
#'                    sig.perc = 0.8, 
#'                    study = 'ICB_Liu')
#'             
#' geneSigSurvDicho(dat.icb = ICB_small_Liu,
#'                  geneSig = sig,
#'                  time.censor = 24,
#'                  n.cutoff = 15,
#'                  n0.cutoff = 5,
#'                  n1.cutoff = 5,
#'                  study =  'ICB_Liu',
#'                  surv.outcome = 'PFS',
#'                  sig.name = 'EMT_Thompson',
#'                  method = 'median',
#'                  var.type = FALSE,
#'                  cancer.type = 'Melanoma',
#'                  treatment = 'PD-1/PD-L1')
#'                 
geneSigSurvDicho <- function(dat.icb, clin = NULL, geneSig, time.censor, n.cutoff, n0.cutoff, n1.cutoff, study, surv.outcome, sig.name,
                             method = "median", var.type, cancer.type, treatment){
  
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
  
  ## association with OS
  if( surv.outcome == "OS"){
    
    if( length( dat_clin$event_occurred_os[ !is.na( dat_clin$event_occurred_os ) ] ) >= n.cutoff ){
      
      cox <- survDicho( status = dat_clin$event_occurred_os ,
                        time = dat_clin$survival_time_os ,
                        time.censor = time.censor ,
                        var = geneSig,
                        n0.cutoff = n0.cutoff,
                        n1.cutoff = n1.cutoff,
                        method = method,
                        var.type = var.type)
      
      res <- data.frame( Outcome = "OS",
                         Gene = sig.name,
                         Study = study,
                         Coef = round(cox["HR"], 3),
                         SE = round(cox["SE"], 3),
                         N = cox["N"],
                         Pval = cox["Pval"],
                         Cancer_type = cancer.type,
                         Treatment = treatment)
      
    }else{  res <- data.frame( Outcome = "OS",
                               Gene = sig.name,
                               Study = study,
                               Coef = NA,
                               SE = NA,
                               N = NA,
                               Pval = NA,
                               Cancer_type = NA,
                               Treatment = NA)
    
    message("lack of number of samples with known immunotherapy survival outcome")
    
    }
    
  }
  
  ## association with PFS
  if( surv.outcome == "PFS"){
    
    if( length( dat_clin$event_occurred_pfs[ !is.na( dat_clin$event_occurred_pfs ) ] ) >= n.cutoff ){
      
      cox <- survDicho( status = dat_clin$event_occurred_pfs ,
                        time = dat_clin$survival_time_pfs ,
                        time.censor = time.censor ,
                        var = geneSig,
                        n0.cutoff = n0.cutoff,
                        n1.cutoff = n1.cutoff,
                        method = method,
                        var.type = var.type)
      
      res <- data.frame( Outcome = "PFS",
                         Gene = sig.name,
                         Study = study,
                         Coef = round(cox["HR"], 3),
                         SE = round(cox["SE"], 3),
                         N = cox["N"],
                         Pval = cox["Pval"],
                         Cancer_type = cancer.type,
                         Treatment = treatment)
      
    }else{  res <- data.frame( Outcome = "PFS",
                               Gene = sig.name,
                               Study = study,
                               Coef = NA,
                               SE = NA,
                               N = NA,
                               Pval = NA,
                               Cancer_type = NA,
                               Treatment = NA)
    
    message("lack of number of samples with known immunotherapy survival outcome")
    
    }
    
  }
  
  return(res)
}
