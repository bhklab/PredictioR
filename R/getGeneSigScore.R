## load libraries

library(GSVA)

########################################################
## create NULL function
########################################################

if_NULL <- function( x , colnames , study ){

  if( is.null( x ) ){

    x = as.data.frame( cbind( study , matrix( NA , nrow= length(study) , ncol=length( colnames )-1 ) ) )

  }
  x

}

##########################################################
##########################################################
# scale data
##########################################################
##########################################################

getScale <- function( x ){

  rid = rownames(x)
  cid = colnames(x)
  out = t( apply( x , 1 , scale ) )
  rownames(out) = rid
  colnames(out) = cid
  out

}

#####################################################################
## Get signature score: GSVA
#####################################################################

getGeneSigGSVA <- function(dat_icb, sig, signature_name, cutoff_n, cutoff_sig, study){


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

    data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type & dat_clin$rna %in% c( "fpkm" , "tpm" )]
    remove <- rem(data)

    if( length(remove) ){
      data <- data[-remove,]
    }

    if( ifelse( is.null( nrow( data[ rownames(data) %in% sig$gene_name , ]) ) , 1 , nrow( data[ rownames(data) %in% sig$gene_name , ] ) ) / length( sig$gene_name ) > cutoff_sig & ncol(data) > cutoff_n ){

      #print( paste( signature_name , "|" , "GSVA" , sep=" " ) )

      geneSig <- NULL
      geneSig <- as.numeric(gsva( getScale( x=data ) , list(sig$gene_name) , verbose=FALSE ) )
      names( geneSig ) <- colnames(data)

    }else{

      geneSig <- NA
      message("not enough samples and/or genes in a data")

    }

  return(geneSig)
}

#####################################################################
## Get signature score: weighted mean
#####################################################################

getGeneSigMean <- function(dat_icb, sig, signature_name, cutoff_n, cutoff_sig, study){

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

    #message(paste(study))

    data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type & dat_clin$rna %in% c( "fpkm" , "tpm" )]
    remove <- rem(data)

    if( length(remove) ){
      data <- data[-remove,]
    }

    if( ifelse( is.null( nrow( data[ rownames(data) %in% sig$gene_name , ]) ) , 1 , nrow( data[ rownames(data) %in% sig$gene_name , ] ) ) / length( sig$gene_name ) > cutoff_sig & ncol(data) > cutoff_n ){

      #print( paste( signature_name , "|" , "Weighted Mean" , sep=" " ) )

      gene <- intersect( rownames(data) , sig$gene_name)
      s <- sig[ sig$gene_name %in% gene, ]

      scaled_dat <- getScale( x= data[ gene , ] )

      remove <- which(is.na(scaled_dat))
      if(length(remove)){
        scaled_dat <- scaled_dat[-which(is.na(scaled_dat)), ]
      }else{ scaled_dat }

      geneSig <- NULL
      geneSig <- apply( scaled_dat , 2 , function(x) ( sum( ( x * s$weight ) ) /  nrow( s ) ) )
      names( geneSig ) = colnames(data)

    }else{

      geneSig <- NA
      message("not enough samples and/or genes in a data")

    }

  return(geneSig)
}

#####################################################################
## Get signature score: PredictIO
#####################################################################

getGeneSigPredictIO <- function(dat_icb, sig, signature_name, cutoff_n, cutoff_sig, study){

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

    #message(paste(study))

    data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type & dat_clin$rna %in% c( "fpkm" , "tpm" )]
    remove <- rem(data)

    if( length(remove) ){
      data <- data[-remove,]
    }

    geneSig <- NULL
    if( ncol(data)){

      #print( paste( signature_name , "|" , "Specific" , sep=" " ) )

      sensitive <- sig[ sig$weight == "sensitive" , ]$gene_name
      resistance <- sig[ sig$weight == "resistance", ]$gene_name

      IO_resistance <- NULL
      if( ifelse( is.null( nrow( data[ rownames(data) %in% resistance , ]) ) , 1 , nrow( data[ rownames(data) %in% resistance , ] ) ) / length( resistance ) > cutoff_sig ){
        IO_resistance <- as.numeric( gsva( getScale( x= data ) , list(resistance) , verbose=FALSE ) )
      }

      IO_sensitive <- NULL
      if( ifelse( is.null( nrow( data[ rownames(data) %in% sensitive , ]) ) , 1 , nrow( data[ rownames(data) %in% sensitive , ] ) ) / length( sensitive ) > cutoff_sig ){
        IO_sensitive <- as.numeric( gsva( getScale( x= data ) , list(sensitive) , verbose=FALSE ) )
      }

      if( !is.null( IO_resistance ) & !is.null( IO_sensitive ) ){

        geneSig <- IO_sensitive - IO_resistance
        names(geneSig) <- colnames(data)

          }else{

            geneSig <- NA
            message("not enough samples and/or genes in a data")

        }

    }else{

      geneSig <- NA
      message("not enough samples in a data")

    }

  return(geneSig)
}


#####################################################################
## Get signature score: COX_IS
#####################################################################

getGeneSigCOX_IS <- function(dat_icb, sig, signature_name, cutoff_n, cutoff_sig, study){


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

    #message(paste(study))

    data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type & dat_clin$rna %in% c( "fpkm" , "tpm" )]
    remove <- rem(data)

    if( length(remove) ){
      data <- data[-remove,]
    }

    if( ifelse( is.null( nrow( data[ rownames(data) %in% sig$gene_name , ]) ) , 1 , nrow( data[ rownames(data) %in% sig$gene_name , ] ) ) / length( sig$gene_name ) > cutoff_sig & ncol(data) > cutoff_n ){

      #print( paste( signature_name , "|" , "Specific" , sep=" " ) )

      geneSig <- NULL
      pos <- apply( data[ rownames(data) %in% sig[ sig$weight %in% 1 , ]$gene_name , ] , 2 , function(x){ ( sum( x ) /  length( x ) ) } )
      neg <- apply( data[ rownames(data) %in% sig[ sig$weight %in% -1 , ]$gene_name , ] , 2 , function(x){ ( sum( x ) /  length( x ) ) } )

      geneSig <- as.numeric( scale( pos / neg ) )
      names(geneSig) <- colnames(data)

      if( length( geneSig[ !is.nan( geneSig ) ] ) < cutoff_n ){

        geneSig <- NA
        message("not enough samples and/or genes in a data")

       }

    }else{

      geneSig <- NA
      message("not enough samples and/or genes in a data")

    }

  return(geneSig)
}

#####################################################################
## Get signature score: IPS
#####################################################################
## calculate Immunophenoscore
ipsmap <- function (x) {
  if (x<=0) {
    ips<-0
  } else {
    if (x>=3) {
      ips<-10
    } else {
      ips<-round(x*10/3, digits=0)
    }
  }
  return(ips)
}

getIPS <- function( expr, sig ){

  expr <- as.data.frame(expr)
  sample_names <- names(expr)

  unique_ips_genes <- as.vector(unique(sig$name))

  IPS <- NULL
  MHC <- NULL
  CP <- NULL
  EC <- NULL
  SC <- NULL
  AZ <- NULL

  # Gene names in expression file
  GVEC <- row.names( expr )
  # Genes names in IPS genes file
  VEC <- as.vector( sig$gene_name )
  # Match IPS genes with genes in expression file
  ind <- which( is.na( match( VEC , GVEC ) ) )

  if( length( ind ) ){
    sig <- sig[-ind,]
  }

  for (i in 1:length(sample_names)) {
    GE <- expr[[i]]
    mGE <- mean(GE, na.rm=T)
    sGE <- sd(GE, na.rm=T)
    Z1 <- (expr[as.vector(sig$gene_name),i] - mGE)/sGE
    W1 <- sig$weight
    coef <- NULL
    MIG <- NULL
    k <- 1

    for (gen in unique_ips_genes) {

      MIG[k] <- mean(Z1[which (as.vector(sig$name)==gen)], na.rm=TRUE)
      coef[k] <- mean(W1[which (as.vector(sig$name)==gen)], na.rm=T)
      k <- k+1

    }

    WG <- MIG * coef
    MHC[i] <- mean(WG[1:10], na.rm=T)
    CP[i] <- mean(WG[11:20], na.rm=T)
    EC[i] <- mean(WG[21:24], na.rm=T)
    SC[i] <- mean(WG[25:26], na.rm=T)
    AZ[i] <- sum(MHC[i], CP[i], EC[i], SC[i])
    IPS[i] <- ipsmap(AZ[i])
  }

  AZ
}

getGeneSigIPS <- function(dat_icb, sig, signature_name, cutoff_n, cutoff_sig, study){


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

    #message(paste(study))

    data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type & dat_clin$rna %in% c( "fpkm" , "tpm" )]
    remove <- rem(data)

    if( length(remove) ){
      data <- data[-remove,]
    }


    geneSig = NULL
    if( ncol(data) & nrow(data) > 10000 ){ # Question: no need to consider the 80% of genes in data?

        #print( paste( signature_name , "|" , "Specific" , sep=" " ) )

        geneSig <- as.numeric( scale( getIPS( expr = data , sig = sig) ) )
        names( geneSig ) <- colnames(data)

      }else{

      geneSig <- NA
      message("not enough samples and/or genes in a data")

      }

  return(geneSig)
}

#####################################################################
## Get signature score: IPRES
#####################################################################

getGeneSigIPRES <- function(dat_icb, sig, signature_name, cutoff_n, cutoff_sig, study){

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


    IPRES <- sig

    sig <- list()
    for( j in 1:length(IPRES)){
      sig[[j]] <-  IPRES[[j]]$gene_name
    }
    names(sig) <- names(IPRES)

    cancer_type <- names( table( dat_clin$cancer_type )[ table( dat_clin$cancer_type ) >= cutoff_n ] )

    #message(paste(study))

    data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type & dat_clin$rna %in% c( "fpkm" , "tpm" )]
    remove <- rem(data)

    if( length(remove) ){
      data <- data[-remove,]
    }

      #print( paste( signature_name , "|" , "Specific" , sep=" " ) )


      scale_data <- getScale( x=data )
      geneSig <- NULL

      for(k in 1:length(sig)){

        if( ifelse( is.null( nrow( scale_data[ rownames(scale_data) %in% sig[[k]] , ]) ) , 1 , nrow( scale_data[ rownames(scale_data) %in% sig[[k]] , ] ) ) / length( sig[[k]] ) >= cutoff_sig & ncol(scale_data) > cutoff_n ){

          geneSig[[k]] <- gsva(scale_data , list(sig[[k]]) , verbose=FALSE)

         }else{

          geneSig[[k]] <- NA

         }

    }

    if( sum(is.na(geneSig)) > 0 & sum(is.na(geneSig)) != length(geneSig) ){ geneSig <- geneSig[!is.na(geneSig)]   }else{

        if(sum(is.na(geneSig)) > 0 & sum(is.na(geneSig)) == length(geneSig)){ geneSig <- NA }else{

          geneSig <- geneSig
        }
      }

      if( sum(!is.na(geneSig)) != 0){

        geneSig <- do.call(rbind, geneSig)

        ## should be mean or is there any other way to integrate signatures? Check the paper.
        geneSig <- apply( geneSig , 2 , mean , na.rm=TRUE )
        names(geneSig) <- colnames(data)

      }else{

        geneSig <- NA
        message("not enough samples and/or genes in a data")
      }

  return(geneSig)
}


