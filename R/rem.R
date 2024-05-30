#' Exclude genes with zero or low expression
#'
#' @param x The expression matrix or data frame with genes as rows and samples as columns. 
#' @param missing.perc A cutoff to remove genes with zero expression across samples.
#' @param const.int A pseudocount is added to the TPM (Transcripts Per Million) values before performing a log transformation.
#' @return Processed expression data after removing low-expressed genes.
#' @export
#'
#' @examples
rem <- function(x, missing.perc, const.int){
  
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  
  # data is log2(TPM+0.001)
  r <- as.numeric(apply(x, 1, function(i) sum(round(i, 6) == round(log2(const.int), 6)) ))
  remove <- which(r > dim(x)[2]* missing.perc)
  return(remove)
  
}