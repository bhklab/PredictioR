#' Meta-analysis via Generic Inverse Variance Method
#' @description
#' Random effect meta-analysis based on estimates (e.g. log hazard ratios) and their standard errors. The inverse variance method is used for pooling.
#' 
#'
#' @param coef A numeric vector of estimates of treatment effect, e.g., log hazard ratio or log odds ratio.
#' @param se A numeric vector of standard errors of treatment estimate.
#' @param study A vector of characters showing the studies' names. 
#' @param pval A numeric vector of estimated p-values.
#' @param n A numeric vector of number of samples.
#' @param cancer.type A vector of characters showing the cancer types.
#' @param treatment A vector of characters showing the treatment types.
#' @param cancer.spec By default, it is FALSE. 
#' @param treatment.spec By default, it is FALSE. 
#' @param feature Name of feature.
#'
#' @return 
#' @export
#'
#' @examples
metafun <- function(coef, se, study, pval, n, cancer.type, treatment, cancer.spec = FALSE, treatment.spec = FALSE, feature){
  
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
    
    if( cancer.spec == TRUE ){
      
      cancer <- cancerfun( data )
      data <- cancer[ order( cancer$Coef ) , ]
      
    }
    
    if( treatment.spec == TRUE ){
      
      treatment <- treatmentfun( data )
      data <- treatment[ order( treatment$Coef ) , ]
      
    }
    
    meta <- metagen( TE = Coef,
                     seTE = SE ,
                     data = data ,
                     studlab = Study ,
                     fixed = FALSE ,
                     random = TRUE ,
                     control = list( maxiter = 10000 , stepadj=0.5 ) )
    
    meta_res <- data.frame(Gene = feature,
                           Coef = meta$TE.random ,
                           SE = meta$seTE.random ,
                           CI_lower = meta$lower.random ,
                           CI_upper = meta$upper.random ,
                           Pval = meta$pval.random ,
                           I2 = meta$I2 ,
                           Q_Pval = meta$pval.Q )
  }else{
    
    print("not enough studies to do meta-analysis")
    meta <- NA
    meta_res <- data.frame(Gene = feature,
                           Coef = NA,
                           SE =  NA,
                           CI_lower = NA,
                           CI_upper = NA,
                           Pval = NA,
                           I2= NA,
                           Q_Pval = NA)
  }
  
  return(list(input_data = data,
              meta_output = meta,
              meta_summery = meta_res))
}
