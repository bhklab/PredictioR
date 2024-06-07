#' Meta-analysis via Generic Inverse Variance Method: Cancer-specific Study
#' @description
#' Random effect meta-analysis based on estimates (e.g. log hazard ratios) and their standard errors. The inverse variance method is used for pooling.
#'
#' @param coef A numeric vector of estimates of treatment effect, e.g., log hazard ratio or log odds ratio.
#' @param se A numeric vector of standard errors of treatment estimate.
#' @param study A vector of characters showing the studies' names. 
#' @param pval A numeric vector of estimated p-values.
#' @param n A numeric vector of number of samples.
#' @param cancer.type A vector of characters showing the cancer types.
#' @param treatment A vector of characters showing the treatment types.
#' @param cancer.spec By default, it is TRUE. 
#' @param treatment.spec By default, it is FALSE. 
#' @param feature Name of feature.
#'
#' @return A list of input data, meta-analysis and summary. 
#' Cancer_type: A character shows the cancer type.
#' Gene: Name of selected feature.
#' Coef: Estimate of treatment effect i.e., log odds ratio or log hazard ratio.
#' SE: Standard error of treatment estimate.
#' CI_lower: Lower bound of the 95% confidence interval.
#' CI_upper: Upper bound of the 95% confidence interval.
#' Pval: Estimated p-value.
#' I2: A percentage of total variation across studies that is due to heterogeneity.
#' Q_Pval: Estimated p-value associated with the Q statistic in meta-analysis.
#' @export
#'
#' @examples
#' expr <- list('ICB_Liu' = ICB_small_Liu, 'ICB_Padron' = ICB_small_Padron, 'ICB_Hugo' = ICB_small_Hugo, 
#'              'ICB_Mariathasan' = ICB_small_Mariathasan, 'ICB_Nathanson' = ICB_small_Nathanson, 
#'              'ICB_Riaz' = ICB_small_Riaz, 'ICB_Miao' = ICB_small_Miao, 'ICB_Van_Allen' = ICB_small_Van_Allen)
#' 
#' cancer_type <- c('Melanoma', 'Pancreas', 'Melanoma', 'Bladder', 'Melanoma', 'Melanoma', 'Kidney', 'Melanoma')
#' treatment_type <- c('PD-1/PD-L1', 'PD-1/PD-L1', 'PD-1/PD-L1', 'PD-1/PD-L1', 'CTLA4', 'IO+combo', 'PD-1/PD-L1', 'CTLA4')
#' 
#' assoc.res <- lapply(1:length(expr), function(k){
#' 
#'  geneSurvCont(dat.icb = expr[[k]],
#'  time.censor = 36,
#'  missing.perc = 0.5,
#'  const.int = 0.001,
#'  n.cutoff = 15,
#'  feature = 'CXCL9',
#'  study = names(expr)[k],
#'  surv.outcome = 'OS',
#'  cancer.type = cancer_type[k],
#'  treatment = treatment_type[k])
#' })
#' assoc.res <- do.call(rbind, assoc.res)
#' 
#' fit <- metaPerCanfun(coef = assoc.res$Coef, 
#'                      se = assoc.res$SE,
#'                      study = assoc.res$Study,
#'                      pval = assoc.res$Pval,
#'                      n = assoc.res$N,
#'                      cancer.type = assoc.res$Cancer_type,
#'                      treatment = assoc.res$Treatment,
#'                      feature = "CXCL9")
#'  
#'  fit$Melanoma$input_data
#'  fit$Melanoma$meta_output
#'  fit$Melanoma$meta_summery  
#'    
metaPerCanfun <- function(coef, se, study, pval, n, cancer.type, treatment, cancer.spec = TRUE, feature){
  
  data <- data.frame( Gene = feature,
                      Study = as.character( study ),
                      N = n,
                      Coef = as.numeric(as.character( coef )),
                      SE = as.numeric(as.character( se )),
                      Pval = as.numeric(as.character( pval )),
                      Cancer_type = as.character( cancer.type ),
                      Treatment = as.character(treatment))
  
  data <- data[ order( data$Coef ) , ]
  
  ## at least 3 studies needed to do the random effect meta-analyses
  if(nrow(data) >= 3){
    
    cancer <- cancerfun( data )
    data <- cancer[ order( cancer$Coef ) , ]
    
    group <- names(table(data$Cancer_type)[table(data$Cancer_type) >= 3])
    
    if( length(group) > 0 ){
      
      res <- lapply(1:length(group), function(i){
        
        sub_data <- data[data$Cancer_type == group[i], ]
        
        meta <- metagen( TE = Coef,
                         seTE = SE ,
                         data = sub_data ,
                         studlab = sub_data$Study ,
                         fixed = FALSE ,
                         random = TRUE ,
                         control = list( maxiter = 10000 , stepadj=0.5 ) )
        
        meta_res <- data.frame(Cancer_type = group[i],
                               Gene = feature,
                               Coef = meta$TE.random ,
                               SE = meta$seTE.random ,
                               CI_lower = meta$lower.random ,
                               CI_upper = meta$upper.random ,
                               Pval = meta$pval.random ,
                               I2 = meta$I2 ,
                               Q_Pval = meta$pval.Q )
        
        list(input_data = sub_data,
             meta_output = meta,
             meta_summery = meta_res)
        
      })
      
      names(res) <- group
      
    }else{
      
      message("not enough studies to do cancer-specific meta-analysis")
      res <- NA
      
    }
    
  }else{
    
    message("not enough studies to do cancer-specific meta-analysis")
    res <- NA
    
  }
  
  return(res)
  
}