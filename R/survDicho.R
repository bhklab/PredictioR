########################################################################################
## Cox model: OS/PFS analyses and binary expression or signature score (low vs high)
########################################################################################
# n0.cutoff: minimum number of samples less than cutoff
# n1.cutoff: minimum number of samples greater than cutoff
# var.type: if variable (var) is dicho (by default), then var.type is TRUE

#' Title
#'
#' @param surv aaaaa
#' @param time bbbbb
#' @param time.censor ccccc
#' @param var dddd
#' @param n0.cutoff eeee
#' @param n1.cutoff ffff
#' @param method gggg
#' @param var.type hhhhh
#'
#' @return iiiii
#' @export
#'
#' @examples
survDicho <- function(surv , time , time.censor , var , n0.cutoff, n1.cutoff, method ="median", var.type = TRUE){
  
  data <- data.frame( surv=surv , time=time , variable=var )
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
      data[ i , "surv" ] = 0
      
    }
  }
  
  if( length( data$variable[ data$variable == 1 ] )>= n1.cutoff &
      length( data$variable[ data$variable == 0 ] ) >= n0.cutoff ){
    
    cox <- coxph( formula= Surv( time , surv ) ~ variable , data=data )
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
