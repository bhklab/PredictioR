#' Compute Signature Score using Single Sample Gene Set Enrichment Analysis (ssGSEA)
#' @description
#' Consider the single sample gene set enrichment analysis to compute signature score. 
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
#' Compute the signature score for M1 Hwang signature using the ssGSEA method. 
#' geneSigssGSEA(dat.icb = ICB_small_Mariathasan, 
#'             sig = M1_Hwang,
#'             sig.name = 'M1_Hwang',
#'             missing.perc = 0.5,
#'             const.int = 0.001,
#'             n.cutoff = 15,
#'             sig.perc = 0.8, 
#'             study = 'ICB_Mariathasan')
#' 
geneSigssGSEA <- function(dat.icb, sig, sig.name, missing.perc, const.int = 0.001, n.cutoff, sig.perc, study){
  
  if( !class(dat.icb) %in% c("SummarizedExperiment", "MultiAssayExperiment", "data.frame", "matrix") ){
    stop(message("function requires SummarizedExperiment, MultiAssayExperiment, data.frame, or matrix class of data"))
  }
  
  if( class(dat.icb) == "MultiAssayExperiment"){
    
    dat <- createSE(dat.icb)
    dat_expr <- assay(dat)
    
  }
  
  if( class(dat.icb) == "SummarizedExperiment"){
    
    dat_expr <- assay(dat.icb)
    dat_clin <- colData(dat.icb)
    
  }
  
  if( class(dat.icb) %in% c("data.frame", "matrix")  ){
    
    dat_expr <- dat.icb
    
  }
  
  data <- dat_expr
  remove <- rem(data, missing.perc, const.int)
  
  if( length(remove) ){
    data <- data[-remove,]
  }
  
  if( ifelse( is.null( nrow( data[ rownames(data) %in% sig$gene_name , ]) ) , 1 , nrow( data[ rownames(data) %in% sig$gene_name , ] ) ) / length( sig$gene_name ) > sig.perc & ncol(data) >= n.cutoff ){
    
    #print( paste( signature_name , "|" , "ssGSEA" , sep=" " ) )
    
    geneSig <- NULL
    gsvaPar <- ssgseaParam(scalefun( x=data ) , list(sig$gene_name))
    geneSig <- gsva(gsvaPar, verbose=FALSE)
    
  }else{
    
    geneSig <- NA
    message("not enough samples and/or genes in a data")
    
  }
  
  return(geneSig)
}
