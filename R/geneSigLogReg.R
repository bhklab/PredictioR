#####################################################################
#####################################################################
## Get gene association (as continuous) with response (R vs NR)
#####################################################################
#####################################################################
# n1.cutoff: cutoff for NR (== 1) samples
# n0.cutoff: cutoff for R (== 0) samples

#' Title
#'
#' @param dat.icb aaaaa
#' @param clin bbbbb
#' @param geneSig ccccc
#' @param n.cutoff ddddd
#' @param study eeeee
#' @param sig.name fffff
#' @param n0.cutoff ggggg
#' @param n1.cutoff hhhhh
#' @param cancer.type iiiii
#' @param treatment jjjjjj
#'
#' @return kkkkk
#' @export
#'
#' @examples
geneSigLogReg <- function(dat.icb, clin = NULL, geneSig, n.cutoff, study, sig.name, n0.cutoff, n1.cutoff, cancer.type, treatment){
  
  if( !class(dat.icb) %in% c("SummarizedExperiment", "MultiAssayExperiment", "data.frame", "matrix") ){
    stop(message("function requires SummarizedExperiment, MultiAssayExperiment, data.frame, or matrix class of data"))
  }
  
  if( class(dat.icb) == "MultiAssayExperiment"){
    dat <- createSE(dat.icb)
    dat_expr <- assay(dat)
    dat_clin <- colData(dat)
  }
  
  if( class(dat.icb) == "SummarizedExperiment"){
    dat_expr <- assay(dat.icb)
    dat_clin <- colData(dat.icb)
  }
  
  if( sum(nrow(clin)) > 0  ){
    
    dat_expr <- dat.icb
    dat_clin <- clin
    
  }
  
  #cancer_type <- names( table( dat_clin$cancer_type )[ table( dat_clin$cancer_type ) >= n.cutoff ] )
  #message(paste(study, cancer_type, sep="/"))
  
  if( length(dat_clin$response) >= n.cutoff &
      sum(dat_clin$response == "NR", na.rm = TRUE) >= n1.cutoff &
      sum(dat_clin$response == "R", na.rm = TRUE) >= n0.cutoff ){
    
    x <- ifelse( dat_clin$response %in% "R" , 0 ,
                 ifelse( dat_clin$response %in% "NR" , 1 , NA ) )
    
    fit <- glm( x ~ geneSig , family=binomial( link="logit" ) )
    
    res <- data.frame( Outcome = "R vs NR",
                       Gene = sig.name,
                       Study = study,
                       Coef = round( summary(fit)$coefficients[ "geneSig" , "Estimate"  ] , 3 ),
                       SE = round( summary(fit)$coefficients[ "geneSig" , "Std. Error" ] , 3 ),
                       N = length(x[!is.na(x)]),
                       Pval = summary(fit)$coefficients[ "geneSig" , "Pr(>|z|)" ],
                       Cancer_type = cancer.type,
                       Treatment = treatment)
    
    
  }else{  res <- data.frame( Outcome = "R vs NR",
                             Gene = sig.name,
                             Study = study,
                             Coef = NA,
                             SE = NA,
                             N = NA,
                             Pval = NA,
                             Cancer_type = NA,
                             Treatment = NA)
  
  message("lack of number of samples with known immunotherapy response outcome")
  
  }
  
  return(res)
  
}

