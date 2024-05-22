##################################################################################################
## Cox model: OS/PFS analyses and continuous expression, signature score, clinical data (var)
##################################################################################################
#' Title
#'
#' @param surv aaaa
#' @param time bbbb
#' @param time.censor cccc
#' @param var dddd
#'
#' @return hhhh
#' @export
#'
#' @examples
survCont <- function( surv , time , time.censor , var){
  
  data <- data.frame( surv=surv , time=time , variable=var )
  data <- data[!is.na(data$variable), ]
  data$time <- as.numeric(as.character(data$time))
  data$variable <- as.numeric( as.character(data$variable) )
  
  for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){
      
      data[ i , "time" ] <- time.censor
      data[ i , "surv" ] <- 0
      
    }
  }
  
  cox <- coxph( Surv( time , surv ) ~ variable, data=data )
  
  hr <- summary(cox)$coefficients[, "coef"]
  se <- summary(cox)$coefficients[, "se(coef)"]
  n <- round(summary(cox)$n)
  low <- summary(cox)$conf.int[, "lower .95"]
  up <- summary(cox)$conf.int[, "upper .95"]
  pval <- summary(cox)$coefficients[, "Pr(>|z|)"]
  
  res <- c( hr , se , n, low , up , pval )
  names(res) <- c( "HR" , "SE" , "N", "Low" , "Up" , "Pval")
  
  return(res)
  
}