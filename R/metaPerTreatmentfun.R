##############################################################
##############################################################
## create meta-analysis per-treatment function
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
#' @param treatment.spec 
#' @param feature 
#'
#' @return
#' @export
#'
#' @examples
metaPerTreatmentfun <- function(coef, se, study, pval, n, cancer.type, treatment, treatment.spec = TRUE, feature){
  
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
    
    treatment <- treatmentfun( data )
    data <- treatment[ order( treatment$Coef ) , ]
    
    group <- names(table(data$Treatment)[table(data$Treatment) >= 3])
    
    if( length(group) > 0){
      
      res <- lapply(1:length(group), function(i){
        
        sub_data <- data[data$Treatment == group[i], ]
        
        meta <- metagen( TE = Coef,
                         seTE = SE ,
                         data = sub_data ,
                         studlab = sub_data$Study ,
                         fixed = FALSE ,
                         random = TRUE ,
                         control = list( maxiter = 10000 , stepadj=0.5 ) )
        
        meta_res <- data.frame(Treatment = group[i],
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

