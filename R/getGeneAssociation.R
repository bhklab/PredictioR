## load functions and libraries

source("C:/PredictioR/R/getHR.R")
source("C:/PredictioR/R/getSummarizedExperiment.R")

############################################################
############################################################
## Get gene association with survival outcome (OS/PFS)
############################################################
############################################################

getGeneAssociationSurvival <- function(dat_icb, time_censor, cutoff_n, feature, study, survival_outcome){


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

        data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type & dat_clin$rna %in% c( "fpkm" , "tpm" )]
        remove <- rem(data)

        if( length(remove) ){
          data <- data[-remove,]
        }

        data <- as.matrix( data[ rownames(data) %in% feature , ] )

        ## association with OS
        if( survival_outcome == "OS"){

          if( nrow(data) & !( cancer_type %in% "Lymph_node" ) &
              length( dat_clin$event_occurred_os[ !is.na( dat_clin$event_occurred_os ) & dat_clin$cancer_type %in% cancer_type ] ) >= cutoff_n ){

            res <- lapply(1:nrow(data), function(k){

              g <- as.numeric( scale( data[k , ] ) )
              names( g ) <- colnames( data )

              cox <- getHRcontinous( surv = dat_clin$event_occurred_os[ dat_clin$cancer_type %in% cancer_type ] ,
                                     time = dat_clin$survival_time_os[ dat_clin$cancer_type %in% cancer_type ] ,
                                     time_censor= time_censor ,
                                     variable= g )

              data.frame( Outcome = "OS",
                          Gene = rownames(data)[k],
                          Study = study,
                          Coef = cox["HR"],
                          SE = cox["SE"],
                          N = cox["N"],
                          Pval = cox["Pval"],
                          Seq = unique(dat_clin$rna) )

            })

            res <- do.call(rbind, res)
            rownames(res) <- NULL
            res$FDR <- p.adjust(res$Pval, method = "BH")

          }else{  res <- data.frame( Outcome = "OS",
                                     Gene = NA,
                                     Study = study,
                                     Coef = NA,
                                     SE = NA,
                                     N = NA,
                                     Pval = NA,
                                     Seq = NA,
                                     FDR =NA)

          message("none of the genes were found in the study and/or lack of number of samples with known immunotherapy survival outcome")

          }

        }

        ## association with PFS
        if( survival_outcome == "PFS"){


          if( nrow(data) & !( cancer_type %in% "Lymph_node" ) &
              length( dat_clin$event_occurred_pfs[ !is.na( dat_clin$event_occurred_pfs ) & dat_clin$cancer_type %in% cancer_type ] ) >= cutoff_n ){

            res <- lapply(1:nrow(data), function(k){

              g <- as.numeric( scale( data[k , ] ) )
              names( g ) <- colnames( data )

              cox <- getHRcontinous( surv = dat_clin$event_occurred_pfs[ dat_clin$cancer_type %in% cancer_type ] ,
                                     time = dat_clin$survival_time_pfs[ dat_clin$cancer_type %in% cancer_type ] ,
                                     time_censor= time_censor , variable= g )

              data.frame( Outcome = "PFS",
                          Gene = rownames(data)[k],
                          Study = study,
                          Coef = cox["HR"],
                          SE = cox["SE"],
                          N = cox["N"],
                          Pval = cox["Pval"],
                          Seq = unique(dat_clin$rna) )

            })

            res <- do.call(rbind, res)
            rownames(res) <- NULL
            res$FDR <- p.adjust(res$Pval, method = "BH")

          }else{  res <- data.frame( Outcome = "PFS",
                                     Gene = NA,
                                     Study = study,
                                     Coef = NA,
                                     SE = NA,
                                     N = NA,
                                     Pval = NA,
                                     Seq = NA,
                                     FDR = NA)

          message("none of the genes were found in the study and/or lack of number of samples with known immunotherapy survival outcome")

          }

        }

  return(res)
}

############################################################
############################################################
## Get gene association with response (R vs NR)
############################################################
############################################################

getGeneAssociationLogReg <- function(dat_icb, cutoff_n, feature, study){

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

    data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type & dat_clin$rna %in% c( "fpkm" , "tpm" )]
    remove <- rem(data)

    if( length(remove) ){
      data <- data[-remove,]
    }

    data <- as.matrix( data[ rownames(data) %in% feature , ] )

    if( nrow(data) & !( cancer_type %in% "Lymph_node" ) ){

      res <- lapply(1:nrow(data), function(k){

        g <- as.numeric( scale( data[k , ] ) )
        names( g ) <- colnames( data )

        ## QUESTION: any filtration based on length of x?
        x <- ifelse( dat_clin$response[ dat_clin$cancer_type %in% cancer_type ] %in% "R" , 0 ,
                     ifelse( dat_clin$response[ dat_clin$cancer_type %in% cancer_type ] %in% "NR" , 1 , NA ) )

        fit <- glm( x ~ g , family=binomial( link="logit" ) )

        data.frame( Outcome = "R vs NR",
                    Gene = rownames(data)[k],
                    Study = study,
                    Coef = round( summary(fit)$coefficients[ "g" , "Estimate"  ] , 3 ),
                    SE = round( summary(fit)$coefficients[ "g" , "Std. Error" ] , 3 ),
                    N = length(x[!is.na(x)]),
                    Pval = summary(fit)$coefficients[ "g" , "Pr(>|z|)" ],
                    Seq = unique(dat_clin$rna) )

      })

      res <- do.call(rbind, res)
      rownames(res) <- NULL
      res$FDR <- p.adjust(res$Pval, method = "BH")

    }else{  res <- data.frame( Outcome = "R vs NR",
                               Gene = NA,
                               Study = study,
                               Coef = NA,
                               SE = NA,
                               N = NA,
                               Pval = NA,
                               Seq = NA,
                               FDR = NA)

    message("none of the genes were found in the study and/or lack of number of samples with known immunotherapy survival outcome")

    }

  return(res)

}

