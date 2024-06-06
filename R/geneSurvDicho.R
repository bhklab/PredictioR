#' Fit Proportional Hazards Regression Model for Genes: Dichotomous Expression Variable
#' @description
#' Fits a Cox proportional hazards regression model with dichotomous expression data for each gene using the counting process formulation of Andersen and Gill.
#'
#' @param dat.icb A MultiAssayExperiment (MAE) object, SummarizedExperiment (SE) object, or a data frame or matrix of gene expression data.
#' @param clin If dat.icb is a data frame or matrix, then it contains clinical data (as data frame or matrix). By default, it is NULL.
#' @param time.censor Possible censoring in months.
#' @param missing.perc A cutoff to remove genes with zero expression across samples.
#' @param const.int A pseudocount is added to the TPM (Transcripts Per Million) values before performing a log transformation.
#' @param n.cutoff Minimum number of samples included in the association analysis.
#' @param feature A vector of character strings for selected features.
#' @param study Name of study.
#' @param surv.outcome Overall survival (OS) or progression-free survival (PFS). 
#' @param n0.cutoff Minimum number of samples with status 0.
#' @param n1.cutoff Minimum number of samples with status 1.
#' @param method The default method to convert a continuous variable into a dichotomous variable is the 'median' method. The first quartile (Q1) and third quartile (Q3) can also be applied. 
#' @param var.type If the variable is dichotomous (by default), then var.type is TRUE.
#' @param cancer.type Name of the cancer type for the given study.
#' @param treatment Name of the treatment for the given study. 
#'
#' @return A subset of results using an object of class 'coxph' representing the fit. 
#' Outcome: Immunotherapy time-to-event outcomes including overall survival (OS) and progression-free survival (PFS).
#' Gene: Name of selected genes.
#' Study: Name of study.
#' Coef: Estimate of treatment effect i.e., log hazard ratio.
#' SE: Standard error of treatment estimate.
#' N: Number of samples.
#' Pval: Estimated p-value.
#' Cancer_type: A character shows the cancer type.
#' Treatment: A character shows the treatment type.
#' @export
#'
#' @examples
#' Assess the association between selected genes (i.e., dichotomous variable with median cutoff) and OS in immunotherapy.
#' geneSurvDicho(dat.icb = ICB_small_Liu,
#'              time.censor = 36,
#'              missing.perc = 0.5,
#'              const.int = 0.001,
#'              n.cutoff = 15,
#'              feature = c('CXCL9', 'CXCL10', 'TIGIT', 
#'                          'CD83', 'STAT1', 'CXCL11',
#'                          'CXCL13', 'CD8A', 'CTLA4'),
#'              study = 'ICB_Liu',
#'              surv.outcome = 'OS',
#'              n0.cutoff = 5, 
#'              n1.cutoff = 5, 
#'              method = "median",
#'              var.type = FALSE,
#'              cancer.type = 'Melanoma',
#'              treatment = 'PD-1/PD-L1')
#'              
geneSurvDicho <- function(dat.icb, clin = NULL, time.censor, missing.perc, const.int=0.001, n.cutoff, feature,
                          study, surv.outcome, n0.cutoff, n1.cutoff, method = "median", var.type, cancer.type, treatment){
  
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
  
  ## association with OS
  if( surv.outcome == "OS"){
    
    if( nrow(data) &  length( dat_clin$event_occurred_os[ !is.na( dat_clin$event_occurred_os ) ] ) >= n.cutoff ){
      
      res <- lapply(1:nrow(data), function(k){
        
        g <- as.numeric( scale( data[k , ] ) )
        names( g ) <- colnames( data )
        
        cox <- survDicho( status = dat_clin$event_occurred_os ,
                          time = dat_clin$survival_time_os ,
                          time.censor = time.censor ,
                          var = g,
                          n0.cutoff = n0.cutoff,
                          n1.cutoff = n1.cutoff,
                          method = method,
                          var.type = var.type)
        
        data.frame( Outcome = "OS",
                    Gene = rownames(data)[k],
                    Study = study,
                    Coef = cox["HR"],
                    SE = cox["SE"],
                    N = cox["N"],
                    Pval = cox["Pval"],
                    Cancer_type = cancer.type,
                    Treatment = treatment)
        
      })
      
      res <- do.call(rbind, res)
      rownames(res) <- NULL
      #  res$FDR <- p.adjust(res$Pval, method = "BH")
      
    }else{  res <- data.frame( Outcome = "OS",
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
    
  }
  
  ## association with PFS
  if( surv.outcome == "PFS"){
    
    
    if( nrow(data) & length( dat_clin$event_occurred_pfs[ !is.na( dat_clin$event_occurred_pfs  ) ] ) >= n.cutoff ){
      
      res <- lapply(1:nrow(data), function(k){
        
        g <- as.numeric( scale( data[k , ] ) )
        names( g ) <- colnames( data )
        
        cox <- survDicho( status = dat_clin$event_occurred_pfs ,
                          time = dat_clin$survival_time_pfs ,
                          time.censor= time.censor ,
                          var = g,
                          n0.cutoff = n0.cutoff,
                          n1.cutoff = n1.cutoff,
                          method = method,
                          var.type = var.type)
        
        data.frame( Outcome = "PFS",
                    Gene = rownames(data)[k],
                    Study = study,
                    Coef = cox["HR"],
                    SE = cox["SE"],
                    N = cox["N"],
                    Pval = cox["Pval"],
                    Cancer_type = cancer.type,
                    Treatment = treatment)
        
      })
      
      res <- do.call(rbind, res)
      rownames(res) <- NULL
      # res$FDR <- p.adjust(res$Pval, method = "BH")
      
    }else{  res <- data.frame( Outcome = "PFS",
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
    
  }
  
  return(res)
}
