##############################################################
##############################################################
## create meta-analysis function
##############################################################
##############################################################

#' Title
#'
#' @param coef 
#' @param se 
#' @param study 
#' @param pval 
#' @param n 
#' @param cancer.type 
#' @param treatment 
#' @param cancer.spec 
#' @param treatment.spec 
#' @param feature 
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
