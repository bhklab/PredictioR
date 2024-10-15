#' Forest Plot: Cancer-specific Analysis
#' @description
#' Draws a forest plot in the active graphics window (using grid graphics system) across cancer types.
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
#' forestPlotPerCan(coef = assoc.res$Coef, 
#'                  se = assoc.res$SE,
#'                  study = assoc.res$Study,
#'                  pval = assoc.res$Pval,
#'                  n = assoc.res$N,
#'                  cancer.type = assoc.res$Cancer_type,
#'                  treatment = assoc.res$Treatment,
#'                  xlab = 'logOR estimate',
#'                  label = 'logOR',
#'                  feature = "CXCL9")
#'            
forestPlotPerCan <- function( coef, se, study, pval, n, cancer.type, treatment, xlab , label, feature){
  
  res <- metafun(coef, se, study, pval, n, cancer.type, treatment, cancer.spec = TRUE, treatment.spec = FALSE, feature)
  cancer <- cancerfun(res$input_data)
  
  remove <- names( table( cancer$Cancer_type)[ table(cancer$Cancer_type) %in% c(1,2) ] )
  
  if(length(table( cancer$Cancer_type )[ table(cancer$Cancer_type) >= 3 ]) > 1){
    
    if( length( unique( cancer$Cancer_type[ !cancer$Cancer_type %in% remove ] ) ) > 1 ){
      
      m <- c( min( c( 0 , cancer$Coef ) , na.rm=TRUE) - .5 , ( max( c( 0 , abs(cancer$Coef) ) , na.rm=TRUE) ) + .5 )
      meta <- res$meta_output
      
      if( length(remove) > 0 ){
        
        meta.subgroup <- update(meta ,
                                byvar = Cancer_type ,
                                exclude = cancer$cancer_type %in% remove ,
                                fixed = FALSE ,
                                random = TRUE ,
                                control = list( maxiter = 10000 , stepadj=0.5 ) )
      } else{
        meta.subgroup <- update(meta ,
                                byvar = Cancer_type ,
                                comb.random = TRUE ,
                                fixed = FALSE ,
                                random = TRUE ,
                                control = list( maxiter = 10000 , stepadj=0.5 ) )
      }
      
    }
    
  }else{
    
    stop("not enough studies to do cancer-specific sub-group meta-analysis")
    
  }
  
  forest( meta.subgroup,
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