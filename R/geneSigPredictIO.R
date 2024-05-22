#####################################################################
## Get signature score: PredictIO
#####################################################################

#' Title
#'
#' @param dat.icb 
#' @param clin 
#' @param sig 
#' @param sig.name 
#' @param missing.perc 
#' @param const.int 
#' @param n.cutoff 
#' @param sig.perc 
#' @param study 
#'
#' @return
#' @export
#'
#' @examples
geneSigPredictIO <- function(dat.icb, clin = NULL, sig, sig.name, missing.perc, const.int =0.001, n.cutoff, sig.perc, study){
  
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
  #message(paste(study))
  #data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type ]
  data <- dat_expr
  remove <- rem(data, missing.perc, const.int)
  
  if( length(remove) ){
    data <- data[-remove,]
  }
  
  geneSig <- NULL
  if( ncol(data)){
    
    #print( paste( signature_name , "|" , "Specific" , sep=" " ) )
    
    sensitive <- sig[ sig$weight == "sensitive" , ]$gene_name
    resistance <- sig[ sig$weight == "resistance", ]$gene_name
    
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

