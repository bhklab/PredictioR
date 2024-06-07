#' Treatment Specific Analysis
#' @description
#' Create a function to specify treatment specific analysis.
#' 
#' @param cancer A data frame containing a subset of association results obtained using logistic regression or 'coxph' models.
#'
#' @return Stratify studies within each treatment type that have at least 3 studies.
#' @export
#'
#' @examples
treatmentfun <- function( treatment ){
  
  treatment$Coef <- as.numeric( as.character( treatment$Coef ) )
  treatment$Study <- as.character( treatment$Study )
  treatment$Cancer_type <- as.character( treatment$Cancer_type )
  treatment$Pval <- as.numeric(as.character( treatment$Pval ))
  treatment$SE <- as.numeric(as.character( treatment$SE ))
  treatment$Study <- as.character(paste(treatment$Study, ", n = ", treatment$N, sep=""))
  
  tab <- table( treatment$Treatment)[ table(treatment$Treatment) %in% c(1,2) ]
  treatment$Treatment[ treatment$Treatment %in% names(tab) ] <- "Other"
  
  treatment
  
}
