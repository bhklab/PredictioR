##############################################################
##############################################################
## create a function to specify treatment analyses
##############################################################
##############################################################

#' Title
#'
#' @param treatment 
#'
#' @return
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
