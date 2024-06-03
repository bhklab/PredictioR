#' Fit Proportional Hazards Regression Model: Dichotomous Expression Variable
#' @description
#' Fits a Cox proportional hazards regression model with dichotomous expression data using the counting process formulation of Andersen and Gill.
#'
#'
#' @param status A vector of 0 and 1, where 0 indicates 'sample was censored at time t' and 1 indicates 'sample had an event at time t'.
#' @param time A vector of time, in months, until endpoint or last follow-up. It is the follow up time (used with right censored data).
#' @param time.censor Possible censoring in months.
#' @param var A vector of dichotomous or continuous expression data.
#' @param n0.cutoff Minimum number of samples with status 0.
#' @param n1.cutoff Minimum number of samples with status 1.
#' @param method The default method to convert a continuous variable into a dichotomous variable is the 'median' method. The first quartile (Q1) and third quartile (Q3) can also be applied. 
#' @param var.type If the variable is dichotomous (by default), then var.type is TRUE.
#'
#' @return A subset of results using an object of class 'coxph' representing the fit.
#' HR: Estimate of treatment effect i.e., log hazard ratio.
#' SE: Standard error of treatment estimate.
#' N: Number of samples.
#' Low: Lower bound of the 95% confidence interval.
#' Up: Upper bound of the 95% confidence interval.
#' Pval: Estimated p-value.
#' 
#' @export
#'
#' @examples
#' Assess the association between CXCL9 (i.e., dichotomous varibale with median cutoff) and OS in immunotherapy.
#' expr <- assay(ICB_Liu);
#' clin <- colData(ICB_Liu) %>% as.data.frame();
#' survDicho( status = clin$event_occurred_os ,
#'            time = clin$survival_time_os,
#'            time.censor= 36,
#'            var = as.numeric(expr["CXCL9", ]),
#'            n0.cutoff = 5,
#'            n1.cutoff = 5,
#'            method = "median",
#'            var.type = FALSE )
#' 
survDicho <- function(status , time , time.censor , var , n0.cutoff, n1.cutoff, method ="median", var.type = TRUE){
  
  data <- data.frame( status=status , time=time , variable=var )
  data <- data[!is.na(data$variable), ]
  data$time <- as.numeric(as.character(data$time))
  
  if(var.type != TRUE){
    
    if( method == "median"){
      data$variable <- ifelse( as.numeric(as.character(data$variable)) >= median(as.numeric(as.character(data$variable))) , 1 , 0 )
    }
    
    if( method == "Q1" ){
      data$variable <- ifelse( as.numeric(as.character(data$variable)) >= quantile(as.numeric(as.character(data$variable)))["25%"] , 1 , 0 )
    }
    
    if( method == "Q3" ){
      data$variable <- ifelse( as.numeric(as.character(data$variable)) >= quantile(as.numeric(as.character(data$variable)))["75%"] , 1 , 0 )
    }
    
  }
  
  for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){
      data[ i , "time" ] = time.censor
      data[ i , "status" ] = 0
      
    }
  }
  
  if( length( data$variable[ data$variable == 1 ] )>= n1.cutoff &
      length( data$variable[ data$variable == 0 ] ) >= n0.cutoff ){
    
    cox <- coxph( formula= Surv( time , status ) ~ variable , data=data )
    hr <- summary(cox)$coefficients[, "coef"]
    se <- summary(cox)$coefficients[, "se(coef)"]
    n <- round(summary(cox)$n)
    low <- summary(cox)$conf.int[, "lower .95"]
    up <- summary(cox)$conf.int[, "upper .95"]
    pval <- summary(cox)$coefficients[, "Pr(>|z|)"]
    
  } else{
    
    hr <- NA
    se = NA
    low <- NA
    up <- NA
    pval <- NA
    
  }
  
  res <- c( hr , se , n, low , up , pval )
  names(res) <- c( "HR" , "SE" , "N", "Low" , "Up" , "Pval")
  return(res)
  
}
