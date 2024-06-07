#' Cancer Specific Analysis
#' @description
#' Create a function to specify cancer specific analysis.
#' 
#' @param cancer A data frame containing a subset of association results obtained using logistic regression or 'coxph' models.
#'
#' @return Stratify studies within each cancer type that have at least 3 studies.
#' @export
#'
#' @examples
cancerfun <- function( cancer ){
  
  cancer$Coef <- as.numeric( as.character( cancer$Coef ) )
  cancer$Study <- as.character( cancer$Study )
  cancer$Cancer_type <- as.character( cancer$Cancer_type )
  cancer$Pval <- as.numeric(as.character( cancer$Pval ))
  cancer$SE <- as.numeric(as.character( cancer$SE ))
  
  tab <- table( cancer$Cancer_type)[ table(cancer$Cancer_type) %in% c(1,2) ]
  cancer$Cancer_type[ cancer$Cancer_type %in% names(tab) ] <- "Other"
  cancer$Study <- as.character(paste(cancer$Study, ", n = ", cancer$N, sep=""))
  
  cancer
  
}