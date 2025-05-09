#' Compute IPS Signature Score 
#' @description
#' Consider the specific method explained at ... to compute IPS signature score. 
#' 
#' @param dat.icb A MultiAssayExperiment (MAE) object, SummarizedExperiment (SE) object, or a data frame or matrix of gene expression data.
#' @param sig A data frame of list of genes' symbol named 'gene_name'. 
#' @param sig.name Name of signature.
#' @param missing.perc A cutoff to remove genes with zero expression across samples.
#' @param const.int A pseudocount is added to the TPM (Transcripts Per Million) values before performing a log transformation.
#' @param n.cutoff Minimum number of samples included in the association analysis.
#' @param study Name of study. 
#'
#' @return A numeric vector of computed signature score.
#' @export
#'
#' @examples
#' Compute the signature score for IPS Charoentong signature using the specific method. 
#' geneSigIPS(dat.icb = ICB_small_Mariathasan, 
#'            sig = IPS_Charoentong,
#'            sig.name = 'IPS_Charoentong',
#'            missing.perc = 0.5,
#'            const.int = 0.001,
#'            n.cutoff = 15,
#'            study = 'ICB_Mariathasan')
#'              
geneSigIPS <- function(dat.icb, sig, sig.name, missing.perc, const.int =0.001, n.cutoff, study, gene.annot = "gene_name"){
  
  
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
  
  
  geneSig = NULL
  if( ncol(data) & nrow(data) > 10000  ){ 
    
    #print( paste( signature_name , "|" , "Specific" , sep=" " ) )
    
    geneSig <- as.numeric( scale( IPSfun( expr = data , sig = sig, gene.annot = gene.annot) ) )
    names( geneSig ) <- colnames(data)
    
  }else{
    
    geneSig <- NA
    message("not enough samples and/or genes in a data")
    
  }
  
  return(geneSig)
}
