##########################################################
##########################################################
# scale data
##########################################################
##########################################################

#' Title
#'
#' @param x 
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
