#' Scaling and Centering the Matrix of Expression Data
#' 
#' @description
#' A generic function whose default method centers and scales the rows (or features) of a numeric expression matrix.
#' 
#' @param x A numeric matrix of expression data.
#'
#' @return
#' @export
#'
#' @examples
scalefun <- function( x ){
  
  rid = rownames(x)
  cid = colnames(x)
  out = t( apply( x , 1 , scale ) )
  rownames(out) = rid
  colnames(out) = cid
  out
  
  }
