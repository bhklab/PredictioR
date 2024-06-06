#' Compute IPSOV Signature Score 
#' @description
#' Consider the specific method explained at ... to compute IPSOV signature score. 
#' 
#' @param dat.icb A MultiAssayExperiment (MAE) object, SummarizedExperiment (SE) object, or a data frame or matrix of gene expression data.
#' @param clin If dat.icb is a data frame or matrix, then it contains clinical data (as data frame or matrix). By default, it is NULL.
#' @param sig A data frame of list of genes' symbol named 'gene_name'. 
#' @param sig.name Name of signature.
#' @param missing.perc A cutoff to remove genes with zero expression across samples.
#' @param const.int A pseudocount is added to the TPM (Transcripts Per Million) values before performing a log transformation.
#' @param n.cutoff Minimum number of samples included in the association analysis.
#' @param sig.perc Minimum percentage of genes in a given expression data. 
#' @param study Name of study. 
#'
#' @return A numeric vector of computed signature score.
#' @export
#'
#' @examples
#' Compute the signature score for IPSOV Shen signature using the specific method. 
#' geneSigIPSOV(dat.icb = ICB_small_Liu, 
#'              sig = IPSOV_Shen,
#'              sig.name = 'IPSOV_Shen',
#'              missing.perc = 0.5,
#'              const.int = 0.001,
#'              n.cutoff = 15,
#'              sig.perc = 0.8, 
#'              study = 'ICB_Liu')
#'              
geneSigIPSOV <- function(dat.icb, clin = NULL, sig, sig.name, missing.perc, const.int =0.001, n.cutoff, sig.perc, study){
  
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
  
  IPSOV.dat <- sig
  
  sig <- list()
  for( j in 1:length(IPSOV.dat)){
    sig[[j]] <-  IPSOV.dat[[j]]$gene_name
  }
  names(sig) <- names(IPSOV.dat)
  
  #cancer_type <- names( table( dat_clin$cancer_type )[ table( dat_clin$cancer_type ) >= n.cutoff ] )
  #message(paste(study))
  #data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type ]
  data <- dat_expr
  remove <- rem(data, missing.perc, const.int)
  
  if( length(remove) ){
    data <- data[-remove,]
  }
  
  #print( paste( signature_name , "|" , "Specific" , sep=" " ) )
  scale.data <- scalefun( data )
  geneSig <- NULL
  
  for(k in 1:length(sig)){
    
    if( ifelse( is.null( nrow( scale.data[ rownames(scale.data) %in% sig[[k]] , ]) ) , 1 , nrow( scale.data[ rownames(scale.data) %in% sig[[k]] , ] ) ) / length( sig[[k]] ) >= sig.perc & ncol(scale.data) >= n.cutoff ){
      
      geneSig[[k]] <- gsva(scale.data , list(sig[[k]]) , method = "ssgsea", verbose=FALSE)
      
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
    geneSig <- apply( geneSig , 2 , mean , na.rm=TRUE )
    names(geneSig) <- colnames(data)
    
  }else{
    
    geneSig <- NA
    message("not enough samples and/or genes in a data")
  }
  
  return(geneSig)
}
