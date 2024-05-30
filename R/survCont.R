#' Fit Proportional Hazards Regression Model: Continuous Expression Variable
#' @description
#' Fits a Cox proportional hazards regression model with continuous expression data using the counting process formulation of Andersen and Gill.
#'     
#' @param status A vector of 0 and 1, where 0 indicates 'sample was censored at time t' and 1 indicates 'sample had an event at time t'.
#' @param time A vector of time, in months, until endpoint or last follow-up. It is the follow up time (used with right censored data).
#' @param time.censor Possible censoring in months.
#' @param var A vector of continuous expression data. 
#'
#' @return A subset of results using an object of class 'coxph' representing the fit.
#' @export 
#'
#' @examples
survCont <- function( status , time , time.censor , var){
  
  data <- data.frame( status=status , time=time , variable=var )
  data <- data[!is.na(data$variable), ]
  data$time <- as.numeric(as.character(data$time))
  data$variable <- as.numeric( as.character(data$variable) )
  
  for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){
      
      data[ i , "time" ] <- time.censor
      data[ i , "status" ] <- 0
      
    }
  }
  
  cox <- coxph( Surv( time , status ) ~ variable, data=data )
  
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