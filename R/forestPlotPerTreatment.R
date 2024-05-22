#########################################################
#########################################################
## create forestplot function: Per treatment analysis
#########################################################
#########################################################

#' Title
#'
#' @param coef 
#' @param se 
#' @param study 
#' @param pval 
#' @param n 
#' @param cancer.type 
#' @param treatment 
#' @param feature 
#' @param xlab 
#' @param label 
#'
#' @return
#' @export
#'
#' @examples
forestPlotPerTreatment <- function( coef, se, study, pval, n , cancer.type, treatment, feature, xlab , label){
  
  res <- metafun(coef, se, study, pval, n, cancer.type, treatment, cancer.spec = FALSE, treatment.spec = TRUE, feature)
  treatment <- treatmentfun(res$input_data)
  
  remove <- names( table( treatment$Treatment )[ table(treatment$Treatment) %in% c(1,2) ] )
  
  if(length(table( treatment$Treatment )[ table(treatment$Treatment) >= 3 ]) > 1){
    
    if( length( unique( treatment$Treatment[ !treatment$Treatment %in% remove ] ) ) > 1 ){
      
      m <- c( min( c( 0 , treatment$Coef ) , na.rm=TRUE) - .5 , ( max( c( 0 , abs(treatment$Coef) ) , na.rm=TRUE) ) + .5 )
      meta <- res$meta_output
      
      if( length(remove) > 0 ){
        
        meta.subgroup <- update(meta ,
                                byvar = Treatment ,
                                exclude = treatment$Treatment %in% remove ,
                                fixed = FALSE ,
                                random = TRUE ,
                                control = list( maxiter = 10000 , stepadj=0.5 ) )
      } else{
        meta.subgroup <- update(meta ,
                                byvar = Treatment ,
                                comb.random = TRUE ,
                                fixed = FALSE ,
                                random = TRUE ,
                                control = list( maxiter = 10000 , stepadj=0.5 ) )
      }
      
    }
    
  }else{
    
    stop("not enough studies to do treatment-specific sub-group meta-analysis")
    
  }
  
  forest( meta.subgroup,
          #leftcols = c("studlab", "effect.ci", "Pval"),
          #leftlabs= c( "Study" , paste(label, "[95%CI]", sep = " ") , "P-value" ) ,
          digits.se = 2,
          colgap.forest=unit(10, "mm") ,
          plotwidth = unit( 30 , "mm") ,
          xlab = xlab,
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



