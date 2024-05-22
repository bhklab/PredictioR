##############################################################
##############################################################
## create forestplot function: Pan cancer analysis
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
#' @param xlab 
#' @param label 
#' @param feature 
#'
#' @return
#' @export
#'
#' @examples
forestPlot <- function(coef, se, study, pval, n , cancer.type, treatment, xlab , label, feature){
  
  study <- paste( study , ", n = " , n , sep= "" )
  res <- metafun(coef, se, study, pval, n, cancer.type, treatment, cancer.spec = FALSE, treatment.spec = FALSE, feature)
  data <- res$input_data
  meta <- res$meta_output
  
  m <- c( min( c( 0 , data$Coef ) , na.rm=TRUE) - .5 , ( max( c( 0 , abs(data$Coef) ) , na.rm=TRUE) ) + .5 )
  
  forest( meta ,
          leftcols = c("studlab", "effect.ci" , "Pval" ),
          leftlabs= c( "Study" , paste(label, "[95%CI]", sep = " ") , "P-value" ) ,
          xlab = xlab,
          digits.se = 2 ,
          colgap.forest=unit(10, "mm") ,
          plotwidth = unit( 30 , "mm") ,
          pooled.totals = TRUE,
          smlab = " ",
          comb.random =TRUE,
          comb.fixed = FALSE,
          text.fixed.w = FALSE,
          layout = "JAMA",
          print.I2.ci = TRUE,
          print.Q = FALSE,
          print.pval.Q = TRUE,
          print.I2 = TRUE,
          print.tau2 = FALSE,
          resid.hetstat = FALSE,
          test.overall.random = TRUE,
          test.overall.fixed = FALSE,
          xlim = m ,
          col.square= "#4d4d4d" ,
          col.study= "#4d4d4d" ,
          col.square.lines = "#4d4d4d" ,
          col.diamond.random  = "#2166ac"  ,
          col.diamond.lines.random  ="#2166ac" ,
          col.by = "#2166ac",
          addrow.subgroups=TRUE )
  
}
