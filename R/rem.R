############################################################################
## Remove genes if with expression zero in 50% (or missing.perc) of sample
############################################################################

#' Title
#'
#' @param x An expression data matrix.
#' @param missing.perc A cutoff for percentage of missing variables across samples.
#' @param const.int A constant value to do log-transformed using TPM data.
#'
#' @return An expression data after removing low-expressed genes.
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