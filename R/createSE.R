#' Convert MultiAssayExperiment Object into SummarizedExperiment Object 
#' @description
#' Convert MultiAssayExperiment (MAE) object into SummarizedExperiment (SE) object which is limited to protein-coding genes. 
#'  
#'  
#'
#' @param dat.icb The MAE object including RNA, DNA, clinical and annotation data.
#'
#' @return A SummarizedExperiment object.
#' @export
#'
createSE <- function(dat.icb){
  
  ## slot names for expression data
  dat_type <- c("expr", "expr_gene_tpm")
  
  dat_type_ix <- which(dat_type %in% names(dat.icb))
  expr <- assay(dat.icb[[dat_type[dat_type_ix]]])
  clin <- as.data.frame(colData(dat.icb[[dat_type[dat_type_ix]]]))
  annot <- as.data.frame(rowData(dat.icb[[dat_type[dat_type_ix]]]))
  
  ## limit to protein-coding genes
  annot <- annot[annot$gene_type == "protein_coding", , drop=FALSE]
  
  ## remove ENSG*Par_Y
  remove_Par_Y <- grep("PAR_Y",rownames(annot))
  if(length(remove_Par_Y) > 0){
    annot <- annot[-remove_Par_Y, ]
  }
  
  ## keep the smallest ENSG for each gene symbol
  annot <- annot[order(rownames(annot)), , drop=FALSE]
  annot <- annot[!duplicated(annot$gene_name), , drop=FALSE]
  
  ## subset gene expression data
  expr <- expr[rownames(annot), , drop=FALSE]
  ## change feature names from ENSG ids to gene symbols
  rownames(expr) <- rownames(annot) <- annot$gene_name
  
  ## create SummarizedExperiments
  eset <- SummarizedExperiment(assay= list("gene_expression"=expr),
                               rowData=annot,
                               colData=clin)
  
  return(eset)
  
}
