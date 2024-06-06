#' Function to Plot Kaplan-Meier Survival Curves
#' @description
#' Function to plot several Kaplan-Meier survival curves with number of individuals at risk at some time points.
#' 
#'
#' @param status A vector of 0 and 1, where 0 indicates 'sample was censored at time t' and 1 indicates 'sample had an event at time t'.
#' @param time A vector of time, in months, until endpoint or last follow-up. It is the follow up time (used with right censored data).
#' @param time.censor Possible censoring in months.
#' @param var A vector of dichotomous or continuous expression data.
#' @param title Label for title of the Kaplan Meier's plot.
#' @param xlab Label for x-axis of the Kaplan Meier's plot.
#' @param ylab Label for y-axis of the Kaplan Meier's plot.
#' @param method The default method to convert a continuous variable into a dichotomous variable is the 'median' method. The first quartile (Q1) and third quartile (Q3) can also be applied. 
#' @param n0.cutoff Minimum number of samples with status 0.
#' @param n1.cutoff Minimum number of samples with status 1.
#' @param var.type If the variable is dichotomous (by default), then var.type is TRUE.
#'
#' @details
#' The original version of this function was kindly provided by Dr Christos Hatzis (January, 17th 2006).
#' 
#' @export
#'
#' @examples
#' Assess the association between CXCL9 (i.e., dichotomous varibale with median cutoff) and PFS in immunotherapy.
#' expr <- assay(ICB_small_Liu);
#' clin <- colData(ICB_small_Liu) %>% as.data.frame();
#' KMPlot(status = clin$event_occurred_pfs,
#'        time = clin$survival_time_pfs,
#'        time.censor = 36,
#'        var =  as.numeric(expr['CXCL9', ]), 
#'        title = "CXCL9 and PFS Association",
#'        xlab = "Time (Months)",
#'        ylab = "Progression-free Survival", 
#'        method = "median",
#'        n0.cutoff = 5, 
#'        n1.cutoff = 5,
#'        var.type = FALSE)

KMPlot <- function( status , time , time.censor , var , title , xlab, ylab,
                    method = "median", n0.cutoff, n1.cutoff, var.type = TRUE){
  
  data <- data.frame( status=status , time=time , variable=var )
  data <- data[!is.na(data$variable), ]
  data$time <- as.numeric(as.character(data$time))
  
  if( var.type != TRUE){
    
    if( method == "median"){
      bin.cutoff <- median(as.numeric(as.character(data$variable)))
      data$variable <- ifelse( as.numeric(as.character(data$variable)) >= bin.cutoff , 1 , 0 )
    }
    
    if( method == "Q1" ){
      bin.cutoff <- quantile(as.numeric(as.character(data$variable)))["25%"]
      data$variable <- ifelse( as.numeric(as.character(data$variable)) >= bin.cutoff , 1 , 0 )
    }
    
    if( method == "Q3" ){
      bin.cutoff <- quantile(as.numeric(as.character(data$variable)))["75%"]
      data$variable <- ifelse( as.numeric(as.character(data$variable)) >= bin.cutoff , 1 , 0 )
    }
    
  }
  
  for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){
      data[ i , "time" ] = time.censor
      data[ i , "status" ] = 0
      
    }
  }
  
  if( length( data$variable[ data$variable == 1 ] )>= n1.cutoff & length( data$variable[ data$variable == 0 ] ) >= n0.cutoff ){
    
    km.coxph.plot( Surv( time, status ) ~ variable , data.s = data, x.label = xlab, y.label = ylab,
                   main.title = paste( title , "\n(cutoff=" , round( bin.cutoff , 2 ) , ")" , sep="" ) ,
                   sub.title = "",
                   leg.text = c( "Low" , "High"),
                   leg.pos = "topright",
                   .col = c( "#b2182b","#2166ac"),
                   show.n.risk = TRUE,
                   n.risk.step = 7,
                   n.risk.cex = 0.70,
                   ylim = c(0,1),
                   leg.inset = 0,
                   .lwd = 2 ,
                   verbose=FALSE )
  }
  
}

