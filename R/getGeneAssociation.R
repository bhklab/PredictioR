## load functions and libraries

source("C:/PredictioR/R/getHR.R")
source("C:/PredictioR/R/getSummarizedExperiment.R")

##########################################################################
##########################################################################
## Get gene association (as continuous) with survival outcome (OS/PFS)
##########################################################################
##########################################################################

geneSurvCont <- function(dat.icb, time.censor, missing.perc, const.int=0.001, n.cutoff, feature, study, surv.outcome){

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

        message(paste(study))

        data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type ]
        remove <- rem(data, missing.perc, const.int)

        if( length(remove) ){
          data <- data[-remove,]
        }

        data <- as.matrix( data[ rownames(data) %in% feature , ] )

        ## association with OS
        if( surv.outcome == "OS"){

          if( nrow(data) & !( cancer_type %in% "Lymph_node" ) &
              length( dat_clin$event_occurred_os[ !is.na( dat_clin$event_occurred_os ) & dat_clin$cancer_type %in% cancer_type ] ) >= n.cutoff ){

            res <- lapply(1:nrow(data), function(k){

              g <- as.numeric( scale( data[k , ] ) )
              names( g ) <- colnames( data )

              cox <- survCont( surv = dat_clin$event_occurred_os[ dat_clin$cancer_type %in% cancer_type ] ,
                               time = dat_clin$survival_time_os[ dat_clin$cancer_type %in% cancer_type ] ,
                               time.censor = time.censor ,
                               var = g )

              data.frame( Outcome = "OS",
                          Gene = rownames(data)[k],
                          Study = study,
                          Coef = cox["HR"],
                          SE = cox["SE"],
                          N = cox["N"],
                          Pval = cox["Pval"],
                          Treatment = unique(dat_clin$treatment))

            })

            res <- do.call(rbind, res)
            rownames(res) <- NULL
            # res$FDR <- p.adjust(res$Pval, method = "BH")

          }else{  res <- data.frame( Outcome = "OS",
                                     Gene = NA,
                                     Study = study,
                                     Coef = NA,
                                     SE = NA,
                                     N = NA,
                                     Pval = NA,
                                     Treatment = NA)

          message("none of the genes were found in the study and/or lack of number of samples with known immunotherapy survival outcome")

          }

        }

        ## association with PFS
        if( surv.outcome == "PFS"){


          if( nrow(data) & !( cancer_type %in% "Lymph_node" ) &
              length( dat_clin$event_occurred_pfs[ !is.na( dat_clin$event_occurred_pfs ) & dat_clin$cancer_type %in% cancer_type ] ) >= n.cutoff ){

            res <- lapply(1:nrow(data), function(k){

              g <- as.numeric( scale( data[k , ] ) )
              names( g ) <- colnames( data )

              cox <- survCont( surv = dat_clin$event_occurred_pfs[ dat_clin$cancer_type %in% cancer_type ] ,
                               time = dat_clin$survival_time_pfs[ dat_clin$cancer_type %in% cancer_type ] ,
                               time.censor = time.censor ,
                               var = g )

              data.frame( Outcome = "PFS",
                          Gene = rownames(data)[k],
                          Study = study,
                          Coef = cox["HR"],
                          SE = cox["SE"],
                          N = cox["N"],
                          Pval = cox["Pval"],
                          Treatment = unique(dat_clin$treatment))

            })

            res <- do.call(rbind, res)
            rownames(res) <- NULL
           # res$FDR <- p.adjust(res$Pval, method = "BH")

          }else{  res <- data.frame( Outcome = "PFS",
                                     Gene = NA,
                                     Study = study,
                                     Coef = NA,
                                     SE = NA,
                                     N = NA,
                                     Pval = NA,
                                     Treatment = NA)

          message("none of the genes were found in the study and/or lack of number of samples with known immunotherapy survival outcome")

          }

        }

  return(res)
}


##########################################################################
##########################################################################
## Get gene association (as continuous) with survival outcome (OS/PFS)
##########################################################################
##########################################################################

geneSurvDicho <- function(dat.icb, time.censor, missing.perc, const.int=0.001, n.cutoff, feature, study, surv.outcome, n0.cutoff, n1.cutoff, method = "median"){

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

  data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type]
  remove <- rem(data, missing.perc, const.int)

  if( length(remove) ){
    data <- data[-remove,]
  }

  data <- as.matrix( data[ rownames(data) %in% feature , ] )

  ## association with OS
  if( surv.outcome == "OS"){

    if( nrow(data) & !( cancer_type %in% "Lymph_node" ) &
        length( dat_clin$event_occurred_os[ !is.na( dat_clin$event_occurred_os ) & dat_clin$cancer_type %in% cancer_type ] ) >= n.cutoff ){

      res <- lapply(1:nrow(data), function(k){

        g <- as.numeric( scale( data[k , ] ) )
        names( g ) <- colnames( data )

        cox <- survDicho( surv = dat_clin$event_occurred_os[ dat_clin$cancer_type %in% cancer_type ] ,
                          time = dat_clin$survival_time_os[ dat_clin$cancer_type %in% cancer_type ] ,
                          time.censor = time.censor ,
                          var = g,
                          n0.cutoff = n0.cutoff,
                          n1.cutoff = n1.cutoff,
                          method = method)

        data.frame( Outcome = "OS",
                    Gene = rownames(data)[k],
                    Study = study,
                    Coef = cox["HR"],
                    SE = cox["SE"],
                    N = cox["N"],
                    Pval = cox["Pval"],
                    Treatment = unique(dat_clin$treatment))

      })

       res <- do.call(rbind, res)
       rownames(res) <- NULL
     #  res$FDR <- p.adjust(res$Pval, method = "BH")

     }else{  res <- data.frame( Outcome = "OS",
                               Gene = NA,
                               Study = study,
                               Coef = NA,
                               SE = NA,
                               N = NA,
                               Pval = NA,
                               Treatment = NA)

    message("none of the genes were found in the study and/or lack of number of samples with known immunotherapy survival outcome")

    }

  }

  ## association with PFS
  if( surv.outcome == "PFS"){


    if( nrow(data) & !( cancer_type %in% "Lymph_node" ) &
        length( dat_clin$event_occurred_pfs[ !is.na( dat_clin$event_occurred_pfs ) & dat_clin$cancer_type %in% cancer_type ] ) >= n.cutoff ){

      res <- lapply(1:nrow(data), function(k){

        g <- as.numeric( scale( data[k , ] ) )
        names( g ) <- colnames( data )

        cox <- survDicho( surv = dat_clin$event_occurred_pfs[ dat_clin$cancer_type %in% cancer_type ] ,
                          time = dat_clin$survival_time_pfs[ dat_clin$cancer_type %in% cancer_type ] ,
                          time.censor= time.censor ,
                          var = g,
                          n0.cutoff = n0.cutoff,
                          n1.cutoff = n1.cutoff,
                          method = method)

        data.frame( Outcome = "PFS",
                    Gene = rownames(data)[k],
                    Study = study,
                    Coef = cox["HR"],
                    SE = cox["SE"],
                    N = cox["N"],
                    Pval = cox["Pval"],
                    Treatment = unique(dat_clin$treatment))

      })

      res <- do.call(rbind, res)
      rownames(res) <- NULL
     # res$FDR <- p.adjust(res$Pval, method = "BH")

    }else{  res <- data.frame( Outcome = "PFS",
                               Gene = NA,
                               Study = study,
                               Coef = NA,
                               SE = NA,
                               N = NA,
                               Pval = NA,
                               Treatment = NA)

    message("none of the genes were found in the study and/or lack of number of samples with known immunotherapy survival outcome")

    }

  }

  return(res)
}

#################################################################
#################################################################
## Get gene association (as continuous) with response (R vs NR)
#################################################################
#################################################################

geneLogReg <- function(dat.icb, missing.perc, const.int=0.001, n.cutoff, feature, study){

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

    data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type]
    remove <- rem(data, missing.perc, const.int)

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
                    Treatment = unique(dat_clin$treatment))

      })

      res <- do.call(rbind, res)
      rownames(res) <- NULL
    #  res$FDR <- p.adjust(res$Pval, method = "BH")

    }else{  res <- data.frame( Outcome = "R vs NR",
                               Gene = NA,
                               Study = study,
                               Coef = NA,
                               SE = NA,
                               N = NA,
                               Pval = NA,
                               Treatment = NA)

    message("none of the genes were found in the study and/or lack of number of samples with known immunotherapy survival outcome")

    }

  return(res)

}

