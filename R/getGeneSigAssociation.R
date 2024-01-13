source("C:/PredictioR/R/getHR.R")
source("C:/PredictioR/R/getSummarizedExperiment.R")

#########################################################################
#########################################################################
## Get gene association (as continuous) with survival outcome (OS/PFS)
#########################################################################
#########################################################################

geneSigSurvCont <- function(dat.icb, geneSig, time.censor, n.cutoff, study, surv.outcome, sig.name){

         if( !class(dat.icb) %in% c("SummarizedExperiment", "MultiAssayExperiment") ){
           stop(message("function requires SummarizedExperiment or MultiAssayExperiment class of data"))
            }

        if( class(dat.icb) == "MultiAssayExperiment"){
          dat <- SummarizedExperiment(dat.icb)
          dat_expr <- assay(dat)
          dat_clin <- colData(dat)
          }

       if( class(dat.icb) == "SummarizedExperiment"){
         dat_expr <- assay(dat.icb)
         dat_clin <- colData(dat.icb)
        }

        cancer_type <- names( table( dat_clin$cancer_type )[ table( dat_clin$cancer_type ) >= n.cutoff ] )

        message(paste(study, cancer_type, sep="/"))

        ## association with OS
        if( surv.outcome == "OS"){

          if( !( cancer_type %in% "Lymph_node" ) &
              length( dat_clin$event_occurred_os[ !is.na( dat_clin$event_occurred_os ) & dat_clin$cancer_type %in% cancer_type ] ) >= n.cutoff ){

              cox <- survCont( surv = dat_clin$event_occurred_os[ dat_clin$cancer_type %in% cancer_type ] ,
                               time = dat_clin$survival_time_os[ dat_clin$cancer_type %in% cancer_type ] ,
                               time.censor = time.censor ,
                               var = geneSig )

              res <- data.frame( Outcome = "OS",
                          Gene = sig.name,
                          Study = study,
                          Coef = cox["HR"],
                          SE = cox["SE"],
                          N = cox["N"],
                          Pval = cox["Pval"],
                          Treatment = unique(dat_clin$treatment))

          }else{  res <- data.frame( Outcome = "OS",
                                     Gene = sig.name,
                                     Study = study,
                                     Coef = NA,
                                     SE = NA,
                                     N = NA,
                                     Pval = NA,
                                     Treatment = NA)

          message("lack of number of samples with known immunotherapy survival outcome")

          }

        }

        ## association with PFS
        if( surv.outcome == "PFS"){

          if( !( cancer_type %in% "Lymph_node" ) &
              length( dat_clin$event_occurred_pfs[ !is.na( dat_clin$event_occurred_pfs ) & dat_clin$cancer_type %in% cancer_type ] ) >= n.cutoff ){

            cox <- survCont( surv = dat_clin$event_occurred_pfs[ dat_clin$cancer_type %in% cancer_type ] ,
                             time = dat_clin$survival_time_pfs[ dat_clin$cancer_type %in% cancer_type ] ,
                             time.censor = time.censor ,
                             var = geneSig )

            res <- data.frame( Outcome = "PFS",
                               Gene = sig.name,
                               Study = study,
                               Coef = cox["HR"],
                               SE = cox["SE"],
                               N = cox["N"],
                               Pval = cox["Pval"],
                               Treatment = unique(dat_clin$treatment))

          }else{  res <- data.frame( Outcome = "PFS",
                                     Gene = sig.name,
                                     Study = study,
                                     Coef = NA,
                                     SE = NA,
                                     N = NA,
                                     Pval = NA,
                                     Treatment = NA)

          message("lack of number of samples with known immunotherapy survival outcome")

          }

        }




  return(res)
}



#########################################################################
#########################################################################
## Get gene association (as binary) with survival outcome (OS/PFS)
#########################################################################
#########################################################################

geneSigSurvDicho <- function(dat.icb, geneSig, time.censor, n.cutoff, n0.cutoff, n1.cutoff, study, surv.outcome, sig.name, method = "median"){

  if( !class(dat.icb) %in% c("SummarizedExperiment", "MultiAssayExperiment") ){
    stop(message("function requires SummarizedExperiment or MultiAssayExperiment class of data"))
  }

  if( class(dat.icb) == "MultiAssayExperiment"){
    dat <- SummarizedExperiment(dat.icb)
    dat_expr <- assay(dat)
    dat_clin <- colData(dat)
  }

  if( class(dat.icb) == "SummarizedExperiment"){
    dat_expr <- assay(dat.icb)
    dat_clin <- colData(dat.icb)
  }

  cancer_type <- names( table( dat_clin$cancer_type )[ table( dat_clin$cancer_type ) >= n.cutoff ] )

  message(paste(study, cancer_type, sep="/"))

  ## association with OS
  if( surv.outcome == "OS"){

    if( !( cancer_type %in% "Lymph_node" ) &
        length( dat_clin$event_occurred_os[ !is.na( dat_clin$event_occurred_os ) & dat_clin$cancer_type %in% cancer_type ] ) >= n.cutoff ){

      cox <- survDicho( surv = dat_clin$event_occurred_os[ dat_clin$cancer_type %in% cancer_type ] ,
                        time = dat_clin$survival_time_os[ dat_clin$cancer_type %in% cancer_type ] ,
                        time.censor = time.censor ,
                        var = geneSig,
                        n0.cutoff = n0.cutoff,
                        n1.cutoff = n1.cutoff,
                        method = method)

      res <- data.frame( Outcome = "OS",
                         Gene = sig.name,
                         Study = study,
                         Coef = cox["HR"],
                         SE = cox["SE"],
                         N = cox["N"],
                         Pval = cox["Pval"],
                         Treatment = unique(dat_clin$treatment))

    }else{  res <- data.frame( Outcome = "OS",
                               Gene = sig.name,
                               Study = study,
                               Coef = NA,
                               SE = NA,
                               N = NA,
                               Pval = NA,
                               Treatment = NA)

    message("lack of number of samples with known immunotherapy survival outcome")

    }

  }

  ## association with PFS
  if( surv.outcome == "PFS"){

    if( !( cancer_type %in% "Lymph_node" ) &
        length( dat_clin$event_occurred_pfs[ !is.na( dat_clin$event_occurred_pfs ) & dat_clin$cancer_type %in% cancer_type ] ) >= n.cutoff ){

      cox <- survDicho( surv = dat_clin$event_occurred_pfs[ dat_clin$cancer_type %in% cancer_type ] ,
                        time = dat_clin$survival_time_pfs[ dat_clin$cancer_type %in% cancer_type ] ,
                        time.censor = time.censor ,
                        var = geneSig,
                        n0.cutoff = n0.cutoff,
                        n1.cutoff = n1.cutoff,
                        method = method)

      res <- data.frame( Outcome = "PFS",
                         Gene = sig.name,
                         Study = study,
                         Coef = cox["HR"],
                         SE = cox["SE"],
                         N = cox["N"],
                         Pval = cox["Pval"],
                         Treatment = unique(dat_clin$treatment))

    }else{  res <- data.frame( Outcome = "PFS",
                               Gene = sig.name,
                               Study = study,
                               Coef = NA,
                               SE = NA,
                               N = NA,
                               Pval = NA,
                               Treatment = NA)

    message("lack of number of samples with known immunotherapy survival outcome")

    }

  }

  return(res)
}

#####################################################################
#####################################################################
## Get gene association (as continuous) with response (R vs NR)
#####################################################################
#####################################################################

geneSigLogReg <- function(dat.icb, geneSig, n.cutoff, study, sig.name){

       if( !class(dat.icb) %in% c("SummarizedExperiment", "MultiAssayExperiment") ){
          stop(message("function requires SummarizedExperiment or MultiAssayExperiment class of data"))
         }

       if( class(dat.icb) == "MultiAssayExperiment"){
         dat <- SummarizedExperiment(dat.icb)
         dat_expr <- assay(dat)
         dat_clin <- colData(dat)
        }

      if( class(dat.icb) == "SummarizedExperiment"){
        dat_expr <- assay(dat.icb)
        dat_clin <- colData(dat.icb)
       }

    cancer_type <- names( table( dat_clin$cancer_type )[ table( dat_clin$cancer_type ) >= n.cutoff ] )

    message(paste(study, cancer_type, sep="/"))

    dat_clin <- dat_clin[!is.na(dat_clin$response), ]
    geneSig <- geneSig[names(geneSig) %in% dat_clin$patientid]

    if( !( cancer_type %in% "Lymph_node" ) & length(dat_clin$response) >= n.cutoff ){

        x <- ifelse( dat_clin$response[ dat_clin$cancer_type %in% cancer_type ] %in% "R" , 0 ,
                     ifelse( dat_clin$response[ dat_clin$cancer_type %in% cancer_type ] %in% "NR" , 1 , NA ) )

        fit <- glm( x ~ geneSig , family=binomial( link="logit" ) )

        res <- data.frame( Outcome = "R vs NR",
                           Gene = sig.name,
                           Study = study,
                           Coef = round( summary(fit)$coefficients[ "geneSig" , "Estimate"  ] , 3 ),
                           SE = round( summary(fit)$coefficients[ "geneSig" , "Std. Error" ] , 3 ),
                           N = length(x[!is.na(x)]),
                           Pval = summary(fit)$coefficients[ "geneSig" , "Pr(>|z|)" ],
                           Treatment = unique(dat_clin$treatment))

    }else{  res <- data.frame( Outcome = "R vs NR",
                               Gene = sig.name,
                               Study = study,
                               Coef = NA,
                               SE = NA,
                               N = NA,
                               Pval = NA,
                               Treatment = NA)

    message("lack of number of samples with known immunotherapy response outcome")

    }

  return(res)

}

