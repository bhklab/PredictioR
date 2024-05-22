##################################################################################################
## Kaplan-Meier (KM) plot: OS/PFS analyses and binary expression or signature score (low vs high)
##################################################################################################
# n0.cutoff: minimum number of samples less than cutoff
# n1.cutoff: minimum number of samples greater than cutoff
# var.type: if variable (var) is dicho (by default), then var.type is TRUE

#' Title
#'
#' @param surv aaaa
#' @param time bbbbb
#' @param time.censor ccccc
#' @param var dddd
#' @param title eeee
#' @param xlab fffff
#' @param ylab gggg
#' @param method hhhh
#' @param n0.cutoff iiii
#' @param n1.cutoff jjjjj
#' @param var.type kkkkk
#'
#' @return llll
#' @export
#'
#' @examples
KMPlot <- function( surv , time , time.censor , var , title , xlab, ylab,
                    method = "median", n0.cutoff, n1.cutoff, var.type = TRUE){
  
  data <- data.frame( surv=surv , time=time , variable=var )
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
      data[ i , "surv" ] = 0
      
    }
  }
  
  if( length( data$variable[ data$variable == 1 ] )>= n1.cutoff & length( data$variable[ data$variable == 0 ] ) >= n0.cutoff ){
    
    km.coxph.plot( Surv( time, surv ) ~ variable , data.s = data, x.label = xlab, y.label = ylab,
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

