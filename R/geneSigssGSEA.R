#' Compute Signature Score using Single Sample Gene Set Enrichment Analysis (ssGSEA)
#' @description
#' Consider the single sample gene set enrichment analysis to compute signature score. 
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
#' Compute the signature score for M1 Hwang signature using the ssGSEA method. 
#' geneSigssGSEA(dat.icb = ICB_small_Liu, 
#'             sig = M1_Hwang,
#'             sig.name = 'M1_Hwang',
#'             missing.perc = 0.5,
#'             const.int = 0.001,
#'             n.cutoff = 15,
#'             sig.perc = 0.8, 
#'             study = 'ICB_Liu')
#' 
geneSigssGSEA <- function(dat.icb, clin = NULL, sig, sig.name, missing.perc, const.int = 0.001, n.cutoff, sig.perc, study){
  
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
  #data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type]
  
  data <- dat_expr
  remove <- rem(data, missing.perc, const.int)
  
  if( length(remove) ){
    data <- data[-remove,]
  }
  
  if( ifelse( is.null( nrow( data[ rownames(data) %in% sig$gene_name , ]) ) , 1 , nrow( data[ rownames(data) %in% sig$gene_name , ] ) ) / length( sig$gene_name ) > sig.perc & ncol(data) >= n.cutoff ){
    
    #print( paste( signature_name , "|" , "ssGSEA" , sep=" " ) )
    
    geneSig <- NULL
    geneSig <- as.numeric(gsva( scalefun( x=data ) , list(sig$gene_name) ,
                                method = "ssgsea", kcdf = "Gaussian", verbose=FALSE ) )
    names( geneSig ) <- colnames(data)
    
  }else{
    
    geneSig <- NA
    message("not enough samples and/or genes in a data")
    
  }
  
  return(geneSig)
}
