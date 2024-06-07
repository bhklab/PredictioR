#' Forest Plot
#' @description
#' Draws a forest plot in the active graphics window (using grid graphics system).
#' 
#' @param coef A numeric vector of estimates of treatment effect, e.g., log hazard ratio or log odds ratio.
#' @param se A numeric vector of standard errors of treatment estimate.
#' @param study A vector of characters showing the studies' names. 
#' @param pval A numeric vector of estimated p-values.
#' @param n A numeric vector of number of samples.
#' @param cancer.type A vector of characters showing the cancer types.
#' @param treatment A vector of characters showing the treatment types.
#' @param xlab A label for x-axis.
#' @param label A label for estimates of treatment effect, e.g., log hazard ratio (logHR).
#' @param feature Name of feature.
#'
#' @return
#' @export
#'
#' @examples
#' expr <- list('ICB_Liu' = ICB_small_Liu, 'ICB_Padron' = ICB_small_Padron, 'ICB_Hugo' = ICB_small_Hugo, 
#'              'ICB_Mariathasan' = ICB_small_Mariathasan, 'ICB_Nathanson' = ICB_small_Nathanson, 
#'              'ICB_Riaz' = ICB_small_Riaz, 'ICB_Miao' = ICB_small_Miao, 'ICB_Van_Allen' = ICB_small_Van_Allen)
#' 
#' cancer_type <- c('Melanoma', 'Pancreas', 'Melanoma', 'Bladder', 'Melanoma', 'Melanoma', 'Kidney', 'Melanoma')
#' treatment_type <- c('PD-1/PD-L1', 'PD-1/PD-L1', 'PD-1/PD-L1', 'PD-1/PD-L1', 'CTLA4', 'IO+combo', 'PD-1/PD-L1', 'CTLA4')
#' 
#' assoc.res <- lapply(1:length(expr), function(k){
#' 
#'  geneSurvCont(dat.icb = expr[[k]],
#'  time.censor = 36,
#'  missing.perc = 0.5,
#'  const.int = 0.001,
#'  n.cutoff = 15,
#'  feature = 'CXCL9',
#'  study = names(expr)[k],
#'  surv.outcome = 'OS',
#'  cancer.type = cancer_type[k],
#'  treatment = treatment_type[k])
#'  
#' })
#' assoc.res <- do.call(rbind, assoc.res)
#' 
#' forestPlot(coef = assoc.res$Coef, 
#'            se = assoc.res$SE,
#'            study = assoc.res$Study,
#'            pval = assoc.res$Pval,
#'            n = assoc.res$N,
#'            cancer.type = assoc.res$Cancer_type,
#'            treatment = assoc.res$Treatment,
#'            xlab = 'logOR estimate',
#'            label = 'logOR',
#'            feature = "CXCL9")
#' 
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
