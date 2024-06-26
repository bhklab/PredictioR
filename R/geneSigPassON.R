#' Compute PassON Signature Score 
#' @description
#' Consider the specific method explained at ... to compute PassON signature score. 
#' 
#' @param dat.icb A MultiAssayExperiment (MAE) object, SummarizedExperiment (SE) object, or a data frame or matrix of gene expression data.
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
#' Compute the signature score for PassON Du signature using the specific method. 
#' geneSigPassON(dat.icb = ICB_small_Mariathasan, 
#'               sig = PassON_Du,
#'               sig.name = 'PassON_Du',
#'               missing.perc = 0.5,
#'               const.int = 0.001,
#'               n.cutoff = 15,
#'               sig.perc = 0.8, 
#'               study = 'ICB_Mariathasan')
#'              
geneSigPassON <- function(dat.icb, sig, sig.name, missing.perc, const.int =0.001, n.cutoff, sig.perc, study){
  
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
  
  PassON.dat <- sig
  
  sig <- list()
  for( j in 1:length(PassON.dat)){
    sig[[j]] <-  PassON.dat[[j]]$gene_name
  }
  names(sig) <- names(PassON.dat)
 
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
      
      gsvaPar <- ssgseaParam(scale.data , list(sig[[k]]))
      genesig <- gsva(gsvaPar, verbose=FALSE)
      geneSig[[k]] <- genesig[1, ]
      
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