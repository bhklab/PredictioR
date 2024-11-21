#' IPS Function
#'
#' @param expr A numeric matrix of expression data.
#' @param sig A data frame of list of genes. 
#' @param gene.annot Specify gene annotation including gene symbol (i.e., gene_name), ENTREZ ID (i.e., entrez_id), and ENSEMBL gene ID (i.e., gene_id).
#' @export
#'
IPSfun <- function( expr, sig, gene.annot ){
  
  expr <- as.data.frame(expr)
  sample_names <- names(expr)
  
  unique_ips_genes <- as.vector(unique(sig$name))
  
  IPS <- NULL
  MHC <- NULL
  CP <- NULL
  EC <- NULL
  SC <- NULL
  AZ <- NULL
  
  # Gene names in expression file
  GVEC <- row.names( expr )
  
  # Genes names in IPS genes file
  if(gene.annot == "gene_name"){ VEC <- as.vector( sig$gene_name ) }
  if( gene.annot == "entrez_id" ){ VEC <- as.vector( sig$entrez_id )  }
  if( gene.annot == "gene_id" ){ VEC <- as.vector( sig$gene_id )  }
  
  # Match IPS genes with genes in expression file
  ind <- which( is.na( match( VEC , GVEC ) ) )
  
  if( length( ind ) ){
    sig <- sig[-ind,]
  }
  
  for (i in 1:length(sample_names)) {
    GE <- expr[[i]]
    mGE <- mean(GE, na.rm=TRUE)
    sGE <- sd(GE, na.rm=TRUE)
    
    if(gene.annot == "gene_name"){ Z1 <- (expr[as.vector(sig$gene_name),i] - mGE)/sGE }
    if( gene.annot == "entrez_id" ){ Z1 <- (expr[as.vector(sig$entrez_id),i] - mGE)/sGE  }
    if( gene.annot == "gene_id" ){ Z1 <- (expr[as.vector(sig$gene_id),i] - mGE)/sGE  }
    
    W1 <- sig$weight
    coef <- NULL
    MIG <- NULL
    k <- 1
    
    for (gen in unique_ips_genes) {
      
      MIG[k] <- mean(Z1[which (as.vector(sig$name)==gen)], na.rm=TRUE)
      coef[k] <- mean(W1[which (as.vector(sig$name)==gen)], na.rm=TRUE)
      k <- k+1
      
    }
    
    WG <- MIG * coef
    MHC[i] <- mean(WG[1:10], na.rm=TRUE)
    CP[i] <- mean(WG[11:20], na.rm=TRUE)
    EC[i] <- mean(WG[21:24], na.rm=TRUE)
    SC[i] <- mean(WG[25:26], na.rm=TRUE)
    AZ[i] <- sum(MHC[i], CP[i], EC[i], SC[i])
    IPS[i] <- ipsmap(AZ[i])
  }
  
  AZ
}
