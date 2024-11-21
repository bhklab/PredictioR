#' Compute Signature Score using Weighted Mean Average
#' @description
#' Consider the weighted mean average to compute signature score. 
#' 
#' @param dat.icb A MultiAssayExperiment (MAE) object, SummarizedExperiment (SE) object, or a data frame or matrix of gene expression data.
#' @param sig A data frame of list of genes' symbol named 'gene_name'. 
#' @param sig.name Name of signature.
#' @param missing.perc A cutoff to remove genes with zero expression across samples.
#' @param const.int A pseudocount is added to the TPM (Transcripts Per Million) values before performing a log transformation.
#' @param n.cutoff Minimum number of samples included in the association analysis.
#' @param sig.perc Minimum percentage of genes in a given expression data. 
#' @param study Name of study. 
#' @param gene.annot Specify gene annotation including gene symbol (i.e., gene_name), ENTREZ ID (i.e., entrez_id), and ENSEMBL gene ID (i.e., gene_id).
#'
#' @return A numeric vector of computed signature score.
#' @export
#'
#' @examples
#' Compute the signature score for EMT Thompson signature using the weighted mean method. 
#' geneSigMean(dat.icb = ICB_small_Mariathasan, 
#'             sig = EMT_Thompson,
#'             sig.name = 'EMT_Thompson',
#'             missing.perc = 0.5,
#'             const.int = 0.001,
#'             n.cutoff = 15,
#'             sig.perc = 0.8, 
#'             study = 'ICB_Mariathasan')
#'              
geneSigMean <- function(dat.icb, sig, sig.name, missing.perc, const.int =0.001, n.cutoff, sig.perc, study, gene.annot = "gene_name"){
  
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
  
  if( gene.annot == "gene_name" ){ genes <- sig$gene_name  }
  if( gene.annot == "entrez_id" ){ genes <- sig$entrez_id  }
  if( gene.annot == "gene_id" ){ genes <- sig$gene_id  }
  
  if( ifelse( is.null( nrow( data[ rownames(data) %in% genes , ]) ) , 1 , nrow( data[ rownames(data) %in% genes , ] ) ) / length( genes ) > sig.perc & ncol(data) >= n.cutoff ){
    
    #print( paste( signature_name , "|" , "Weighted Mean" , sep=" " ) )
    
    gene <- intersect( rownames(data) , genes)
    s <- sig[ genes %in% gene, ]
    
    scaled_dat <- scalefun( x= data[ gene , ] )
    
    remove <- which(is.na(scaled_dat))
    if(length(remove)){
      scaled_dat <- scaled_dat[-which(is.na(scaled_dat)), ]
    }else{ scaled_dat }
    
    geneSig <- NULL
    geneSig <- apply( scaled_dat , 2 , function(x) ( sum( x * s$weight, na.rm=TRUE  )  /  nrow( s ) ) )
    names( geneSig ) <- colnames(data)
    
  }else{
    
    geneSig <- NA
    message("not enough samples and/or genes in a data")
    
  }
  
  return(geneSig)
}

