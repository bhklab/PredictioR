source("C:/PredictioR/R/getHR.R")
source("C:/PredictioR/R/getSummarizedExperiment.R")

#########################################################################
#########################################################################
## Get gene association (as continuous) with survival outcome (OS/PFS)
#########################################################################
#########################################################################

getGeneSigAssociationSurvival <- function(dat_icb, geneSig, time_censor, cutoff_n, study, survival_outcome, signature_name){

         if( !class(dat_icb) %in% c("SummarizedExperiment", "MultiAssayExperiment") ){
           stop(message("function requires SummarizedExperiment or MultiAssayExperiment class of data"))
            }

        if( class(dat_icb) == "MultiAssayExperiment"){
          dat <- getSummarizedExperiment(dat_icb)
          dat_expr <- assay(dat)
          dat_clin <- colData(dat)
          }

       if( class(dat_icb) == "SummarizedExperiment"){
         dat_expr <- assay(dat_icb)
         dat_clin <- colData(dat_icb)
        }

        cancer_type <- names( table( dat_clin$cancer_type )[ table( dat_clin$cancer_type ) >= cutoff_n ] )

        message(paste(study, cancer_type, sep="/"))

        ## association with OS
        if( survival_outcome == "OS"){

          if( !( cancer_type %in% "Lymph_node" ) &
              length( dat_clin$event_occurred_os[ !is.na( dat_clin$event_occurred_os ) & dat_clin$cancer_type %in% cancer_type ] ) >= cutoff_n ){

              cox <- getHRcontinous( surv = dat_clin$event_occurred_os[ dat_clin$cancer_type %in% cancer_type ] ,
                                     time = dat_clin$survival_time_os[ dat_clin$cancer_type %in% cancer_type ] ,
                                     time_censor= time_censor ,
                                     variable= geneSig )

              res <- data.frame( Outcome = "OS",
                          Gene = signature_name,
                          Study = study,
                          Coef = cox["HR"],
                          SE = cox["SE"],
                          N = cox["N"],
                          Pval = cox["Pval"],
                          Seq = unique(dat_clin$rna) )

          }else{  res <- data.frame( Outcome = "OS",
                                     Gene = signature_name,
                                     Study = study,
                                     Coef = NA,
                                     SE = NA,
                                     N = NA,
                                     Pval = NA,
                                     Seq = NA)

          message("lack of number of samples with known immunotherapy survival outcome")

          }

        }

        ## association with PFS
        if( survival_outcome == "PFS"){

          if( !( cancer_type %in% "Lymph_node" ) &
              length( dat_clin$event_occurred_pfs[ !is.na( dat_clin$event_occurred_pfs ) & dat_clin$cancer_type %in% cancer_type ] ) >= cutoff_n ){

            cox <- getHRcontinous( surv = dat_clin$event_occurred_pfs[ dat_clin$cancer_type %in% cancer_type ] ,
                                   time = dat_clin$survival_time_pfs[ dat_clin$cancer_type %in% cancer_type ] ,
                                   time_censor= time_censor ,
                                   variable= geneSig )

            res <- data.frame( Outcome = "PFS",
                               Gene = signature_name,
                               Study = study,
                               Coef = cox["HR"],
                               SE = cox["SE"],
                               N = cox["N"],
                               Pval = cox["Pval"],
                               Seq = unique(dat_clin$rna) )

          }else{  res <- data.frame( Outcome = "PFS",
                                     Gene = signature_name,
                                     Study = study,
                                     Coef = NA,
                                     SE = NA,
                                     N = NA,
                                     Pval = NA,
                                     Seq = NA)

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

getGeneSigAssociationSurvivalDicho <- function(dat_icb, geneSig, time_censor, cutoff_n, cutoff_n0, cutoff_n1, study, survival_outcome, signature_name){

  if( !class(dat_icb) %in% c("SummarizedExperiment", "MultiAssayExperiment") ){
    stop(message("function requires SummarizedExperiment or MultiAssayExperiment class of data"))
  }

  if( class(dat_icb) == "MultiAssayExperiment"){
    dat <- getSummarizedExperiment(dat_icb)
    dat_expr <- assay(dat)
    dat_clin <- colData(dat)
  }

  if( class(dat_icb) == "SummarizedExperiment"){
    dat_expr <- assay(dat_icb)
    dat_clin <- colData(dat_icb)
  }

  cancer_type <- names( table( dat_clin$cancer_type )[ table( dat_clin$cancer_type ) >= cutoff_n ] )

  message(paste(study, cancer_type, sep="/"))

  ## association with OS
  if( survival_outcome == "OS"){

    if( !( cancer_type %in% "Lymph_node" ) &
        length( dat_clin$event_occurred_os[ !is.na( dat_clin$event_occurred_os ) & dat_clin$cancer_type %in% cancer_type ] ) >= cutoff_n ){

      cox <- getHRdicho( surv = dat_clin$event_occurred_os[ dat_clin$cancer_type %in% cancer_type ] ,
                         time = dat_clin$survival_time_os[ dat_clin$cancer_type %in% cancer_type ] ,
                         time_censor= time_censor ,
                         variable= geneSig,
                         cutoff_n0 = cutoff_n0,
                         cutoff_n1 = cutoff_n1)

      res <- data.frame( Outcome = "OS",
                         Gene = signature_name,
                         Study = study,
                         Coef = cox["HR"],
                         SE = cox["SE"],
                         N = cox["N"],
                         Pval = cox["Pval"],
                         Seq = unique(dat_clin$rna) )

    }else{  res <- data.frame( Outcome = "OS",
                               Gene = signature_name,
                               Study = study,
                               Coef = NA,
                               SE = NA,
                               N = NA,
                               Pval = NA,
                               Seq = NA)

    message("lack of number of samples with known immunotherapy survival outcome")

    }

  }

  ## association with PFS
  if( survival_outcome == "PFS"){

    if( !( cancer_type %in% "Lymph_node" ) &
        length( dat_clin$event_occurred_pfs[ !is.na( dat_clin$event_occurred_pfs ) & dat_clin$cancer_type %in% cancer_type ] ) >= cutoff_n ){

      cox <- getHRdicho( surv = dat_clin$event_occurred_pfs[ dat_clin$cancer_type %in% cancer_type ] ,
                         time = dat_clin$survival_time_pfs[ dat_clin$cancer_type %in% cancer_type ] ,
                         time_censor= time_censor ,
                         variable= geneSig,
                         cutoff_n0 = cutoff_n0,
                         cutoff_n1 = cutoff_n1)

      res <- data.frame( Outcome = "PFS",
                         Gene = signature_name,
                         Study = study,
                         Coef = cox["HR"],
                         SE = cox["SE"],
                         N = cox["N"],
                         Pval = cox["Pval"],
                         Seq = unique(dat_clin$rna) )

    }else{  res <- data.frame( Outcome = "PFS",
                               Gene = signature_name,
                               Study = study,
                               Coef = NA,
                               SE = NA,
                               N = NA,
                               Pval = NA,
                               Seq = NA)

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

getGeneSigAssociationLogReg <- function(dat, geneSig, cutoff_n, study, signature_name){

       if( !class(dat_icb) %in% c("SummarizedExperiment", "MultiAssayExperiment") ){
          stop(message("function requires SummarizedExperiment or MultiAssayExperiment class of data"))
         }

       if( class(dat_icb) == "MultiAssayExperiment"){
         dat <- getSummarizedExperiment(dat_icb)
         dat_expr <- assay(dat)
         dat_clin <- colData(dat)
        }

      if( class(dat_icb) == "SummarizedExperiment"){
        dat_expr <- assay(dat_icb)
        dat_clin <- colData(dat_icb)
       }

    cancer_type <- names( table( dat_clin$cancer_type )[ table( dat_clin$cancer_type ) >= cutoff_n ] )

    message(paste(study, cancer_type, sep="/"))

    if( !( cancer_type %in% "Lymph_node" ) ){

        x <- ifelse( dat_clin$response[ dat_clin$cancer_type %in% cancer_type ] %in% "R" , 0 ,
                     ifelse( dat_clin$response[ dat_clin$cancer_type %in% cancer_type ] %in% "NR" , 1 , NA ) )

        fit <- glm( x ~ geneSig , family=binomial( link="logit" ) )

        res <- data.frame( Outcome = "R vs NR",
                    Gene = signature_name,
                    Study = paste( study, cancer_type, sep="__" ),
                    Coef = round( summary(fit)$coefficients[ "geneSig" , "Estimate"  ] , 3 ),
                    SE = round( summary(fit)$coefficients[ "geneSig" , "Std. Error" ] , 3 ),
                    N = length(x[!is.na(x)]),
                    Pval = summary(fit)$coefficients[ "geneSig" , "Pr(>|z|)" ],
                    Seq = unique(dat_clin$rna) )

    }else{  res <- data.frame( Outcome = "R vs NR",
                               Gene = signature_name,
                               Study = study,
                               Coef = NA,
                               SE = NA,
                               N = NA,
                               Pval = NA,
                               Seq = NA)

    message("lack of number of samples with known immunotherapy response outcome")

    }

  return(res)

}

