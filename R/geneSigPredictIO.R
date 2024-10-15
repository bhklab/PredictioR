#' Compute PredictIO Signature Score 
#' @description
#' Consider the specific method explained at ... to compute PredictIO signature score. 
#' 
#' @param dat.icb A MultiAssayExperiment (MAE) object, SummarizedExperiment (SE) object, or a data frame or matrix of gene expression data.
#' @param sig A data frame of list of genes' symbol named 'gene_name'. 
#' @param sig.name Name of signature.
#' @param missing.perc A cutoff to remove genes with zero expression across samples.
#' @param const.int A pseudocount is added to the TPM (Transcripts Per Million) values before performing a log transformation.
#' @param n.cutoff Minimum number of samples included in the association analysis.
#' @param sig.perc Minimum percentage of genes in a given expression data. 
#' @param study Name of study. 
#' @param gene_annot Specify gene annotation including gene symbol (i.e., gene_name), ENTREZ ID (i.e., entrez_id), and ENSEMBL gene ID (i.e., gene_id).
#'
#' @return A numeric vector of computed signature score.
#' @export
#'
#' @examples
#' Compute the signature score for PredictIO Bareche signature using the specific method. 
#' geneSigPredictIO(dat.icb = ICB_small_Mariathasan, 
#'                  sig = PredictIO_Bareche,
#'                  sig.name = 'PredictIO_Bareche',
#'                  missing.perc = 0.5,
#'                  const.int = 0.001,
#'                  n.cutoff = 15,
#'                  sig.perc = 0.8, 
#'                  study = 'ICB_Mariathasan')
#'              
geneSigPredictIO <- function(dat.icb, sig, sig.name, missing.perc, const.int =0.001, n.cutoff, sig.perc, study, gene_annot = "gene_name"){
  
  if( !class(dat.icb) %in% c("SummarizedExperiment", "MultiAssayExperiment", "data.frame", "matrix") ){
    stop(message("function requires SummarizedExperiment, MultiAssayExperiment, data.frame, or matrix class of data"))
  }
  
  if( class(dat.icb) == "MultiAssayExperiment"){
    
    dat <- createSE(dat.icb)
    dat_expr <- assay(dat)
    
  }
  
  if( class(dat.icb) == "SummarizedExperiment"){
    
    dat_expr <- assay(dat.icb)
    
  }
  
  if( class(dat.icb) %in% c("data.frame", "matrix")  ){
    
    dat_expr <- dat.icb
    
  }
  
  data <- dat_expr
  remove <- rem(data, missing.perc, const.int)
  
  if( length(remove) ){
    data <- data[-remove,]
  }
  
  geneSig <- NULL
  if( ncol(data)){
    
    #print( paste( signature_name , "|" , "Specific" , sep=" " ) )
    
    if( gene_annot == "gene_name" ){   
      sensitive <- sig[ sig$weight == "sensitive" , ]$gene_name
      resistance <- sig[ sig$weight == "resistance", ]$gene_name
    }
    
    
    if( gene_annot == "entrez_id" ){   
      sensitive <- sig[ sig$weight == "sensitive" , ]$entrez_id
      resistance <- sig[ sig$weight == "resistance", ]$entrez_id
    }
    
    
    if( gene_annot == "gene_id" ){   
      sensitive <- sig[ sig$weight == "sensitive" , ]$gene_id
      resistance <- sig[ sig$weight == "resistance", ]$gene_id
    }
    
    
    IO_resistance <- NULL
    if( ifelse( is.null( nrow( data[ rownames(data) %in% resistance , ]) ) , 1 , nrow( data[ rownames(data) %in% resistance , ] ) ) / length( resistance ) > sig.perc ){
      gsvaPar <- gsvaParam(scalefun( x=data ) , list(resistance))
      IO_resistance <- as.numeric(gsva(gsvaPar, verbose=FALSE))
    }
    
    IO_sensitive <- NULL
    if( ifelse( is.null( nrow( data[ rownames(data) %in% sensitive , ]) ) , 1 , nrow( data[ rownames(data) %in% sensitive , ] ) ) / length( sensitive ) > sig.perc ){
      gsvaPar <- gsvaParam(scalefun( x=data ) , list(sensitive))
      IO_sensitive <- as.numeric(gsva(gsvaPar, verbose=FALSE))
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

