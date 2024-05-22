#####################################################################
## Get signature score: COX_IS
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
geneSigCOX_IS <- function(dat.icb, clin = NULL, sig, sig.name, missing.perc, const.int =0.001, n.cutoff, sig.perc, study){
  
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
  
  if( ifelse( is.null( nrow( data[ rownames(data) %in% sig$gene_name , ]) ) , 1 , nrow( data[ rownames(data) %in% sig$gene_name , ] ) ) / length( sig$gene_name ) > sig.perc & ncol(data) >= n.cutoff ){
    
    #print( paste( signature_name , "|" , "Specific" , sep=" " ) )
    
    geneSig <- NULL
    pos <- apply( data[ rownames(data) %in% sig[ sig$weight %in% 1 , ]$gene_name , ] , 2 , function(x){ ( sum( x) /  length( x) ) } )
    neg <- apply( data[ rownames(data) %in% sig[ sig$weight %in% -1 , ]$gene_name , ] , 2 , function(x){ ( sum( x) /  length( x) ) } )
    
    geneSig <- as.numeric( scale( pos / neg ) )
    names(geneSig) <- colnames(data)
    
    if( length( geneSig[ !is.nan( geneSig ) ] ) < n.cutoff ){
      
      geneSig <- NA
      message("not enough samples and/or genes in a data")
      
    }
    
  }else{
    
    geneSig <- NA
    message("not enough samples and/or genes in a data")
    
  }
  
  return(geneSig)
}
