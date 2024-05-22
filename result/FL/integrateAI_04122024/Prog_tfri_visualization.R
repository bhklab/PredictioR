## librries

library(stringr)
library(survcomp)
library(GSVA)
library(dplyr)
library(meta)
library(metafor)
library(forestplot)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(data.table)
library(kableExtra)
library(summarytools)
library(MultiAssayExperiment)
library(pheatmap)
library(RColorBrewer)
library(ggVennDiagram)
library(ggvenn)
library(VennDiagram)

app_dir <- str_split(rstudioapi::getActiveDocumentContext()$path,'PredictioR')[[1]][1]

source(file.path(app_dir, 'PredictioR', "R/getHR.R"))
source(file.path(app_dir, 'PredictioR', "R/getSummarizedExperiment.R"))
source(file.path(app_dir, 'PredictioR', "R/getMetaAnalysis.R"))

## load data
dir <- "C:/PredictioR/result/FL/integrateAI_04122024"

res.pan <- read.csv(file.path(dir, "pan_cancer.csv"))
res.percan <- read.csv(file.path(dir, "per_cancer.csv"))
res.pertreatment <- read.csv(file.path(dir, "per_treatment.csv"))

volcanoPlot <- function(feature, coef, pval, padj, pos.cutoff, neg.cutoff, x.lab, padj.label, cutoff){

  data <- data.frame(feature = feature,
                     coef = coef,
                     pval = pval,
                     FDR = padj)

  if( padj.label == FALSE){

    data$diffexpressed <- "NO"
    data$diffexpressed[data$coef > 0 & data$pval < cutoff] <- paste(paste("Pval < ", cutoff, sep=""), "Coef > 0", sep=", ")
    data$diffexpressed[data$coef < 0 & data$pval < cutoff] <- paste(paste("Pval < ", cutoff, sep=""), "Coef < 0", sep=", ")

    mycolors <- c( "#4d004b","#014636", "#999999")
    names(mycolors) <- c(paste(paste("Pval < ", cutoff, sep=""), "Coef > 0", sep=", "),
                         paste(paste("Pval < ", cutoff, sep=""), "Coef < 0", sep=", "),
                         "NO")

    data$delabel <- NA
    data <- data[order(data$pval, decreasing = FALSE), ]
    id_pos <- data[data$coef > 0 , "feature"][1:pos.cutoff]
    id_neg <- data[data$coef < 0 , "feature"][1:neg.cutoff]
    id <- c(id_pos, id_neg)

    for(j in 1:length(id)){
      k <- which(data$feature == id[j])
      data$delabel[k] <- data[k, ]$feature
    }

  }else{

    data$diffexpressed <- "NO"
    data$diffexpressed[data$coef > 0 & data$FDR < cutoff] <- paste(paste("FDR < ", cutoff, sep=""), "Coef > 0", sep=", ")
    data$diffexpressed[data$coef < 0 & data$FDR < cutoff] <- paste(paste("FDR < ", cutoff, sep=""), "Coef < 0", sep=", ")

    mycolors <- c( "#4d004b","#014636", "#999999")
    names(mycolors) <- c(paste(paste("FDR < ", cutoff, sep=""), "Coef > 0", sep=", "),
                         paste(paste("FDR < ", cutoff, sep=""), "Coef < 0", sep=", "),
                         "NO")

    data$delabel <- NA
    data <- data[order(data$FDR, decreasing = FALSE), ]
    id_pos <- data[data$coef > 0 , "feature"][1:pos.cutoff]
    id_neg <- data[data$coef < 0 , "feature"][1:neg.cutoff]
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
      legend.text = element_text(size = 6, face="bold"),
      legend.title = element_blank()) +
    geom_text_repel(aes(label= delabel),
                    size = 2.2,
                    color = "black",
                    min.segment.length = 0,
                    na.rm = TRUE, direction = "both",
                    seed = 2356,
                    fontface= "bold",
                    max.overlaps = 50)
}

##################################################################################
################################### per cancer ###################################
##################################################################################
# Kidney
df <- res.percan[res.percan$Cancer_type == "Kidney", ]
df$FDR.FL <- p.adjust(df$Pval.FL, method = "BH")

jpeg(file = file.path("C:/PredictioR/result/FL/integrateAI_04122024/Fig",
                        paste(paste("volcano_response", "Kidney", sep="_"), ".jpeg", sep="")),
       width = 500, height = 450, res = 150)

pdf(file = file.path("C:/PredictioR/result/FL/integrateAI_04122024/Fig",
                     paste(paste("volcano_response", "Kidney", sep="_"), ".pdf", sep="")),
    width = 3.5, height = 3)

volcanoPlot(feature = df$Gene,
                   coef = df$coef.FL,
                   pval = df$Pval.FL,
                   padj = df$FDR.FL,
                   neg.cutoff = 3,
                   pos.cutoff = 3,
                   x.lab = "logOR estimate",
                   padj.label = FALSE,
                   cutoff = 0.05)

dev.off()

# Melanoma
df <- res.percan[res.percan$Cancer_type =="Melanoma", ]
df$FDR.FL <- p.adjust(df$Pval.FL, method = "BH")

jpeg(file = file.path("C:/PredictioR/result/FL/integrateAI_04122024/Fig",
                        paste(paste("volcano_response", "Melanoma", sep="_"), ".jpeg", sep="")),
       width = 500, height = 450, res = 150)

pdf(file = file.path("C:/PredictioR/result/FL/integrateAI_04122024/Fig",
                     paste(paste("volcano_response", "Melanoma", sep="_"), ".pdf", sep="")),
    3.5, height = 3)

volcanoPlot(feature = df$Gene,
              coef = df$coef.FL,
              pval = df$Pval.FL,
              padj = df$FDR.FL,
              neg.cutoff = 10,
              pos.cutoff = 1,
              x.lab = "logOR estimate",
              padj.label = FALSE,
              cutoff = 0.05)

dev.off()


# Lung
df <- res.percan[res.percan$Cancer_type =="Lung", ]
df$FDR.FL <- p.adjust(df$Pval.FL, method = "BH")

jpeg(file = file.path("C:/PredictioR/result/FL/integrateAI_04122024/Fig",
                      paste(paste("volcano_response", "Lung", sep="_"), ".jpeg", sep="")),
     width = 500, height = 450, res = 150)

pdf(file = file.path("C:/PredictioR/result/FL/integrateAI_04122024/Fig",
                     paste(paste("volcano_response", "Lung", sep="_"), ".pdf", sep="")),
    width = 3.5, height = 3)

volcanoPlot(feature = df$Gene,
            coef = df$coef.FL,
            pval = df$Pval.FL,
            padj = df$FDR.FL,
            neg.cutoff = 1,
            pos.cutoff = 1,
            x.lab = "logOR estimate",
            padj.label = FALSE,
            cutoff = 0.05)

dev.off()

##################################################################################
################################### pan cancer ###################################
##################################################################################
df <- res.pan
df$FDR.FL <- p.adjust(df$Pval.FL, method = "BH")

jpeg(file = file.path("C:/PredictioR/result/FL/integrateAI_04122024/Fig",
                      paste(paste("volcano_response", "pan", sep="_"), ".jpeg", sep="")),
     width = 500, height = 450, res = 150)

pdf(file = file.path("C:/PredictioR/result/FL/integrateAI_04122024/Fig",
                     paste(paste("volcano_response", "pan", sep="_"), ".pdf", sep="")),
    width = 3.5, height = 3)

volcanoPlot(feature = df$Gene,
            coef = df$coef.FL,
            pval = df$Pval.FL,
            padj = df$FDR.FL,
            neg.cutoff = 10,
            pos.cutoff = 1,
            x.lab = "logOR estimate",
            padj.label = TRUE,
            cutoff = 0.05)

dev.off()

##################################################################################
################################### per treatment ################################
##################################################################################
df <- res.pertreatment
df$FDR.FL <- p.adjust(df$Pval.FL, method = "BH")

jpeg(file = file.path("C:/PredictioR/result/FL/integrateAI_04122024/Fig",
                      paste(paste("volcano_response", "PD-(L)1", sep="_"), ".jpeg", sep="")),
     width = 500, height = 450, res = 150)

pdf(file = file.path("C:/PredictioR/result/FL/integrateAI_04122024/Fig",
                      paste(paste("volcano_response", "PD-(L)1", sep="_"), ".pdf", sep="")),
     width = 3.5, height = 3)

volcanoPlot(feature = df$Gene,
            coef = df$coef.FL,
            pval = df$Pval.FL,
            padj = df$FDR.FL,
            neg.cutoff = 10,
            pos.cutoff = 1,
            x.lab = "logOR estimate",
            padj.label = TRUE,
            cutoff = 0.05)

dev.off()

####################################################################################
## Forestplot
####################################################################################
## merge association results

dir <- "C:/PredictioR/result/FL/DHDP_03272024/result/association"
files <- list.files(dir)

res.assoc <- lapply(1:length(files), function(k){

  read.csv(file.path(dir, files[k]))

})

res.assoc <- do.call(rbind, res.assoc)

# pan-cancer
res.meta <- res.pan[order(res.pan$Pval.FL), ]
res.meta <- res.meta[res.meta$I2 < 0.50, ]
res.meta <- res.meta[res.meta$Pval.FL < 0.05 & res.meta$Pval < 0.05,]
selected.sig <- res.meta$Gene
dat <- res.assoc[res.assoc$Gene %in% selected.sig, ]

for(i in 1:length(selected.sig)){

  df <- dat[dat$Gene == selected.sig[i], ]

  jpeg(file=file.path("C:/PredictioR/result/FL/integrateAI_04122024/Fig",
                      paste(selected.sig[i], "forestplot_pan_logreg.jpeg", sep="_")),
       width = 1050, height = 800, res = 150)

  #pdf(file=file.path("C:/PredictioR/result/FL/integrateAI_04122024/Fig",
  #                    paste(selected.sig[i], "forestplot_pan_logreg.pdf", sep="_")),
  #     width = 8, height = 8)

  forestPlot(coef = df$Coef,
                   se = df$SE,
                   study  = df$Study,
                   pval = df$Pval,
                   n = df$N,
                   cancer.type = do.call(rbind, strsplit(df$Study, split='__', fixed=TRUE))[, 2],
                   treatment = df$Treatment,
                   feature = unique(df$Gene),
                   xlab = "logOR estimate",
                   label = "logOR")

  dev.off()

}

# per-cancer
res.meta <- res.percan[order(res.percan$Pval.FL), ]
res.meta <- res.meta[res.meta$I2 < 0.50, ]
res.meta <- res.meta[res.meta$Pval.FL < 0.05 & res.meta$Pval < 0.05,]
selected.sig <- res.meta[res.meta$Cancer_type == "Melanoma", "Gene"]
dat <- res.assoc[res.assoc$Gene %in% selected.sig, ]
dat <- dat[dat$Cancer_type == "Melanoma", ]
dat <- dat[!is.na(dat$Gene), ]

for(i in 1:length(selected.sig)){

  df <- dat[dat$Gene == selected.sig[i], ]

    jpeg(file=file.path("C:/PredictioR/result/FL/integrateAI_04122024/Fig",
                        paste(selected.sig[i], "forestplot_percan_logreg.jpeg", sep="_")),
         width = 1050, height = 600, res = 150)

    #pdf(file=file.path("C:/PredictioR/result/FL/integrateAI_04122024/Fig",
    #                   paste(selected.sig[i], "forestplot_percan_logreg.pdf", sep="_")),
    #    width = 8, height = 7)

    forestPlot(coef = df$Coef,
               se = df$SE,
               study  = df$Study,
               pval = df$Pval,
               n = df$N,
               cancer.type = do.call(rbind, strsplit(df$Study, split='__', fixed=TRUE))[, 2],
               treatment = df$Treatment,
               feature = unique(df$Gene),
               xlab = "logOR estimate",
               label = "logOR")

    dev.off()

}

# per-treatment
res.meta <- res.pertreatment[order(res.pertreatment$Pval.FL), ]
res.meta <- res.meta[res.meta$I2 < 0.50, ]
res.meta <- res.meta[res.meta$Pval.FL < 0.05 & res.meta$Pval < 0.05,]
selected.sig <- res.meta$Gene[1:10]
dat <- res.assoc[res.assoc$Gene %in% selected.sig, ]
dat <- dat[dat$Treatment == "PD-(L)1", ]

for(i in 1:length(selected.sig)){

  df <- dat[dat$Gene == selected.sig[i], ]

  jpeg(file=file.path("C:/PredictioR/result/FL/integrateAI_04122024/Fig",
                      paste(selected.sig[i], "forestplot_percan_logreg.jpeg", sep="_")),
       width = 1050, height = 700, res = 150)

  #pdf(file=file.path("C:/PredictioR/result/FL/integrateAI_04122024/Fig",
  #                  paste(selected.sig[i], "forestplot_pertreatment_logreg.pdf", sep="_")),
  #    width = 8, height = 7)

  forestPlot(coef = df$Coef,
                         se = df$SE,
                         study  = df$Study,
                         pval = df$Pval,
                         n = df$N,
                         cancer.type = do.call(rbind, strsplit(df$Study, split='__', fixed=TRUE))[, 2],
                         treatment = df$Treatment,
                         feature = unique(df$Gene),
                         xlab = "logOR estimate",
                         label = "logOR")

  dev.off()

}

####################################################################################
## Heatmap
####################################################################################

os <- res.pan
os$FDR.FL <- p.adjust(os$Pval.FL, method = "BH")
os <- os[os$FDR.FL < 0.05 | os$Pval.FL < 0.05, ]
os <- os[os$FDR < 0.05 | os$Pval < 0.05, ]
os <- os[!is.na(os$Coef), ]

pfs <- res.percan[res.percan$Cancer_type == "Melanoma", ]
pfs$FDR.FL <- p.adjust(pfs$Pval.FL, method = "BH")
pfs <- pfs[pfs$FDR.FL < 0.05 | pfs$Pval.FL < 0.05, ]
pfs <- pfs[pfs$FDR < 0.05 | pfs$Pval < 0.05, ]
pfs <- pfs[!is.na(pfs$Coef), ]

logreg <- res.pertreatment
logreg$FDR.FL <- p.adjust(logreg$Pval.FL, method = "BH")
logreg <- logreg[logreg$FDR.FL < 0.05 | logreg$Pval.FL < 0.05, ]
logreg <- logreg[logreg$FDR < 0.05 | logreg$Pval < 0.05, ]
logreg <- logreg[!is.na(logreg$Coef), ]

int <- union(union(pfs$Gene, logreg$Gene), os$Gene)

os_mod <- lapply(1:length(int), function(k){

  if( sum(os$Gene %in% int[k])>0 ){  res <- os[os$Gene %in% int[k], c("Gene", "coef.FL", "se.FL", "Pval.FL", "FDR.FL")] }else{

    res <-data.frame(Gene = int[k],
                     coef.FL = NA,
                     se.FL =NA,
                     Pval.FL = NA,
                     FDR.FL = NA)
  }

  res
})

os <- do.call(rbind, os_mod)
os <- os[order(os$Gene), ]

pfs_mod <- lapply(1:length(int), function(k){

  if( sum(pfs$Gene %in% int[k])>0 ){  res <- pfs[pfs$Gene %in% int[k], c("Gene", "coef.FL", "se.FL", "Pval.FL", "FDR.FL")] }else{

    res <-data.frame(Gene = int[k],
                     coef.FL = NA,
                     se.FL =NA,
                     Pval.FL = NA,
                     FDR.FL = NA)
  }

  res
})

pfs <- do.call(rbind, pfs_mod)
pfs <- pfs[order(pfs$Gene), ]

logreg_mod <- lapply(1:length(int), function(k){

  if( sum(logreg$Gene %in% int[k])>0 ){  res <- logreg[logreg$Gene %in% int[k], c("Gene", "coef.FL", "se.FL", "Pval.FL", "FDR.FL")] }else{

    res <-data.frame(Gene = int[k],
                     coef.FL = NA,
                     se.FL =NA,
                     Pval.FL = NA,
                     FDR.FL = NA)
  }

  res
})

logreg <- do.call(rbind, logreg_mod)
logreg <- logreg[order(logreg$Gene), ]

data = cbind( logreg$coef.FL , pfs$coef.FL, os$coef.FL )
rownames(data) = logreg$Gene
colnames(data) = c( "Treatment" , "Melanoma", "PanCancer" )

pval = cbind( logreg$Pval.FL , pfs$Pval.FL, os$Pval.FL )
rownames(pval) = logreg$Gene
colnames(pval) = c( "Treatment" , "Melanoma", "PanCancer" )
pval[ is.na(pval) ] = 1

padj = cbind( logreg$FDR.FL , pfs$FDR.FL, os$FDR.FL )
rownames(padj) = logreg$Gene
colnames(padj) = c( "Treatment" , "Melanoma", "PanCancer" )
padj[ is.na(padj) ] = 1

annot_col = data.frame(
  "PanCancer" = factor( ifelse( round( padj[ , 'PanCancer' ] , 2 ) <= .05 , 'FDR' , ifelse( round( pval[ , "PanCancer" ] , 2 ) <= 0.05 , "Pvalue" ,  'NS' ) ) ) ,
  Melanoma = factor( ifelse( round( padj[ , 'Melanoma' ] , 2 ) <= .05 , 'FDR' , ifelse( round( pval[ , "Melanoma" ] , 2 ) <= 0.05 , "Pvalue" , 'NS' ) ) ) ,
  "Treatment" = factor( ifelse( round( padj[ , 'Treatment' ] , 2 ) <= .05 , 'FDR' , ifelse( round( pval[ , "Treatment" ] , 2 ) <= 0.05 , "Pvalue" , 'NS' ) ) )
)

rownames(annot_col) = rownames(data)

ann_colors = list(
  "Treatment" = c( FDR = "#636363", Pvalue = "#bdbdbd", NS = "#f0f0f0" ) ,
  Melanoma = c( FDR = "#636363", Pvalue = "#bdbdbd", NS = "#f0f0f0" ),
  "PanCancer" = c( FDR = "#636363", Pvalue = "#bdbdbd", NS = "#f0f0f0" )
)

neg = seq( round( min( data , na.rm=TRUE ) , 1 ) , 0 , by=.05 )
neg = neg[ -length(neg)]
pos = seq( 0 , round( max( data , na.rm=TRUE ) , 1 ) , by=.05 )

col = c( colorRampPalette( rev( brewer.pal(4, "Blues") ) )( length(neg) ) ,
         colorRampPalette( brewer.pal(4, "OrRd") )( length(pos) )
)

#"#4d004b","#014636"

# filename="C://PredictioR/result/Roche/Fig/Heatmap_pan.jpeg", height=4, width=4
jpeg(file = "C:/PredictioR/result/FL/integrateAI_04122024/Fig/heatmap_pan_melanoma_treatment.jpeg",
     width = 800, height = 350, res = 150)

#pdf(file = "C:/PredictioR/result/FL/integrateAI_04122024/Fig/heatmap_pan_melanoma_treatment.pdf",
#     width = 6, height = 2.2)


df <- t( data[ order( rowSums(data)) , ] )
#annot_row <- annot_row[rownames(df), ]
annot_col <- annot_col[colnames(df), ]



pheatmap( df , cluster_rows=FALSE , cluster_cols=FALSE , scale="none" ,
          #annotation_row = annot_row,
          annotation_col = annot_col,
          annotation_colors = ann_colors,
          name= "Coef",
          col = col , breaks = c( neg , pos ) , na_col="#f0f0f0" , border_color="#424242",
          number_color="black" , show_colnames = T, show_rownames = T,
          fontsize = 7, fontsize_number = 6)

dev.off()
















