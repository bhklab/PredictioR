#' Volcano Plot
#' @description
#' Creates a volcano plot of log-odds or log-hazard ratios of changes versus the log p-values from association analysis.
#' 
#' @param feature A vector of characters showing the features' names.
#' @param coef A numeric vector of log-odds or log-hazard ratio estimates.
#' @param pval A numeric vector of estimated p-values.
#' @param padj A numeric vector of adjusted p-values to correct for multiple tests using the false discovery rate (FDR).
#' @param pos.cutoff Number of selected top features with positive log-odds or log-hazard ratio estimates.
#' @param neg.cutoff Number of selected top features with negative log-odds or log-hazard ratio estimates.
#' @param x.lab Label of x-axis.
#' @param padj.label If the adjusted p-values (FDR) are used to select significantly associated features, then it is TRUE.
#' @param cutoff Cut-off for adjusted p-values (FDR) or p-values.
#' @param colors A vector of colors to identify positive, negative, and not associated features.
#' @param coef.cutoff Cut-off for estimated coef i.e., log-odds or log-hazard ratio.
#'
#' @export
#'
volcanoPlot <- function(feature, coef, pval, padj, pos.cutoff, neg.cutoff, x.lab, padj.label, cutoff, colors, coef.cutoff){
  
  data <- data.frame(feature = feature,
                     coef = coef,
                     pval = pval,
                     FDR = padj)
  
  if( padj.label == FALSE){
    
    data$diffexpressed <- "NO"
    data$diffexpressed[data$coef > coef.cutoff & data$pval < cutoff] <- paste(paste("Pval < ", cutoff, sep=""), paste("Coef > ", coef.cutoff, sep=""), sep=", ")
    data$diffexpressed[data$coef < coef.cutoff & data$pval < cutoff] <- paste(paste("Pval < ", cutoff, sep=""), paste("Coef < ", coef.cutoff, sep=""), sep=", ")
    
    mycolors <- colors
    names(mycolors) <- c(paste(paste("Pval < ", cutoff, sep=""), paste("Coef > ", coef.cutoff, sep=""), sep=", "),
                         paste(paste("Pval < ", cutoff, sep=""), paste("Coef < ", coef.cutoff, sep=""), sep=", "),
                         "NO")
    
    data$delabel <- NA
    data <- data[order(data$pval, decreasing = FALSE), ]
    id_pos <- data[data$coef > coef.cutoff , "feature"][1:pos.cutoff]
    id_neg <- data[data$coef < coef.cutoff , "feature"][1:neg.cutoff]
    id <- c(id_pos, id_neg)
    
    for(j in 1:length(id)){
      k <- which(data$feature == id[j])
      data$delabel[k] <- data[k, ]$feature
    }
    
  }else{
    
    data$diffexpressed <- "NO"
    data$diffexpressed[data$coef > coef.cutoff & data$FDR < cutoff] <- paste(paste("FDR < ", cutoff, sep=""), paste("Coef > ", coef.cutoff, sep=""), sep=", ")
    data$diffexpressed[data$coef < coef.cutoff & data$FDR < cutoff] <- paste(paste("FDR < ", cutoff, sep=""), paste("Coef < ", coef.cutoff, sep=""), sep=", ")
    
    mycolors <- colors
    names(mycolors) <- c(paste(paste("FDR < ", cutoff, sep=""), paste("Coef > ", coef.cutoff, sep=""), sep=", "),
                         paste(paste("FDR < ", cutoff, sep=""), paste("Coef < ", coef.cutoff, sep=""), sep=", "),
                         "NO")
    
    data$delabel <- NA
    data <- data[order(data$FDR, decreasing = FALSE), ]
    id_pos <- data[data$coef > coef.cutoff , "feature"][1:pos.cutoff]
    id_neg <- data[data$coef < coef.cutoff , "feature"][1:neg.cutoff]
    id <- c(id_pos, id_neg)
    
    for(j in 1:length(id)){
      k <- which(data$feature == id[j])
      data$delabel[k] <- data[k, ]$feature
    }
    
  }
  
  ggplot(data=data, aes(x=coef, y=-log10(pval), col= diffexpressed)) +
    geom_point(size = 2.5) + theme_minimal() +
    ylab("-log10 P value") +
    xlab(x.lab) +
    scale_colour_manual(values = mycolors) +
    theme(
      axis.text.x=element_text(size=10,  face="bold"),
      axis.title=element_text(size=10,face="bold"),
      axis.text.y=element_text(size=10, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position="bottom",
      legend.text = element_text(size = 7, face="bold"),
      legend.title = element_blank()) +
    geom_text_repel(aes(label= delabel),
                    size = 2.4,
                    color = "black",
                    min.segment.length = 0,
                    na.rm = TRUE, direction = "both",
                    seed = 2356,
                    fontface= "bold",
                    max.overlaps = 50)
  
 }
