#########################################################################
#########################################################################
## Get gene association (as binary) with survival outcome (OS/PFS)
#########################################################################
#########################################################################

#' Title
#'
#' @param dat.icb aaaaa
#' @param clin bbbbb
#' @param geneSig ccccc
#' @param time.censor ddddd
#' @param n.cutoff eeeee
#' @param n0.cutoff fffff
#' @param n1.cutoff ggggg
#' @param study hhhhh
#' @param surv.outcome iiii
#' @param sig.name jjjjj
#' @param method kkkkk
#' @param var.type lllll
#' @param cancer.type mmmm
#' @param treatment nnnnn
#'
#' @return ooooo
#' @export
#'
#' @examples
geneSigSurvDicho <- function(dat.icb, clin = NULL, geneSig, time.censor, n.cutoff, n0.cutoff, n1.cutoff, study, surv.outcome, sig.name,
                             method = "median", var.type, cancer.type, treatment){
  
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
  
  ## association with OS
  if( surv.outcome == "OS"){
    
    if( length( dat_clin$event_occurred_os[ !is.na( dat_clin$event_occurred_os ) ] ) >= n.cutoff ){
      
      cox <- survDicho( surv = dat_clin$event_occurred_os ,
                        time = dat_clin$survival_time_os ,
                        time.censor = time.censor ,
                        var = geneSig,
                        n0.cutoff = n0.cutoff,
                        n1.cutoff = n1.cutoff,
                        method = method,
                        var.type = var.type)
      
      res <- data.frame( Outcome = "OS",
                         Gene = sig.name,
                         Study = study,
                         Coef = round(cox["HR"], 3),
                         SE = round(cox["SE"], 3),
                         N = cox["N"],
                         Pval = cox["Pval"],
                         Cancer_type = cancer.type,
                         Treatment = treatment)
      
    }else{  res <- data.frame( Outcome = "OS",
                               Gene = sig.name,
                               Study = study,
                               Coef = NA,
                               SE = NA,
                               N = NA,
                               Pval = NA,
                               Cancer_type = NA,
                               Treatment = NA)
    
    message("lack of number of samples with known immunotherapy survival outcome")
    
    }
    
  }
  
  ## association with PFS
  if( surv.outcome == "PFS"){
    
    if( length( dat_clin$event_occurred_pfs[ !is.na( dat_clin$event_occurred_pfs ) ] ) >= n.cutoff ){
      
      cox <- survDicho( surv = dat_clin$event_occurred_pfs ,
                        time = dat_clin$survival_time_pfs ,
                        time.censor = time.censor ,
                        var = geneSig,
                        n0.cutoff = n0.cutoff,
                        n1.cutoff = n1.cutoff,
                        method = method,
                        var.type = var.type)
      
      res <- data.frame( Outcome = "PFS",
                         Gene = sig.name,
                         Study = study,
                         Coef = round(cox["HR"], 3),
                         SE = round(cox["SE"], 3),
                         N = cox["N"],
                         Pval = cox["Pval"],
                         Cancer_type = cancer.type,
                         Treatment = treatment)
      
    }else{  res <- data.frame( Outcome = "PFS",
                               Gene = sig.name,
                               Study = study,
                               Coef = NA,
                               SE = NA,
                               N = NA,
                               Pval = NA,
                               Cancer_type = NA,
                               Treatment = NA)
    
    message("lack of number of samples with known immunotherapy survival outcome")
    
    }
    
  }
  
  return(res)
}
