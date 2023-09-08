## load libraries

library(survcomp)
library(GSVA)
library(MultiAssayExperiment)
library(qs)
library(data.table)

############################################################
## Remove genes if with expression zero in 50% of sample
############################################################

rem <- function(x){
  
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  
  # data is log2(TPM+0.001)
  r <- as.numeric(apply(x, 1, function(i) sum(round(i, 6) == round(log2(0.001), 6)) )) 
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
  
 }

########################################################################################
## Cox model: OS/PFS analyses and continuous expression or signature score
########################################################################################

getHRcontinous <- function( surv , time , time_censor , variable){
	
  data <- data.frame( surv=surv , time=time , variable=variable )
	data$time <- as.numeric(as.character(data$time))
	data$variable <- as.numeric( as.character(data$variable) )

  	for(i in 1:nrow(data)){
  	  
	    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time_censor ){
	     
	       data[ i , "time" ] <- time_censor
	       data[ i , "surv" ] <- 0
    	
	       }
  	}

	cox <- coxph( Surv( time , surv ) ~ variable, data=data )
	
	hr <- summary(cox)$coefficients[, "coef"] 
	se <- summary(cox)$coefficients[, "se(coef)"]
	n <- round(summary(cox)$n)
	low <- summary(cox)$conf.int[, "lower .95"] 
	up <- summary(cox)$conf.int[, "upper .95"]
	pval <- summary(cox)$coefficients[, "Pr(>|z|)"]

	res <- c( hr , se , n, low , up , pval )
	names(res) <- c( "HR" , "SE" , "N", "Low" , "Up" , "Pval")
	
	return(res)
	
 }

########################################################################################
## Cox model: OS/PFS analyses and binary expression or signature score (low vs high)
########################################################################################

getHRdicho <- function( surv , time , time_censor , variable , cutoff ){
 
  data <- data.frame( surv=surv , time=time , variable=variable )
  data$time <- as.numeric(as.character(data$time))
  data$variable <- ifelse( as.numeric(as.character(data$variable)) >= cutoff , 1 , 0 )
  
  for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time_censor ){
      data[ i , "time" ] = time_censor
      data[ i , "surv" ] = 0
      
    }
  }
  
  if( length( data$variable[ data$variable == 1 ] )>=5 & length( data$variable[ data$variable == 0 ] ) >= 5 ){
    
    cox <- coxph( formula= Surv( time , surv ) ~ variable , data=data )
    hr <- summary(cox)$coefficients[, "coef"] 
    se <- summary(cox)$coefficients[, "se(coef)"]
    n <- round(summary(cox)$n)
    low <- summary(cox)$conf.int[, "lower .95"] 
    up <- summary(cox)$conf.int[, "upper .95"]
    pval <- summary(cox)$coefficients[, "Pr(>|z|)"]
    
  } else{
    
    hr <- NA
    se = NA
    low <- NA
    up <- NA
    pval <- NA	
    
  }	
  
  res <- c( hr , se , n, low , up , pval )
  names(res) <- c( "HR" , "SE" , "N", "Low" , "Up" , "Pval")
  return(res)
  
 }

##################################################################################################
## Kaplan-Meier (KM) plot: OS/PFS analyses and binary expression or signature score (low vs high)
##################################################################################################

getKMplot <- function( surv , time , time_censor , variable , cutoff , title , xlab, ylab ){
  
  data <- data.frame( surv=surv , time=time , variable=variable )
  data$time <- as.numeric(as.character(data$time))
  data$variable <- ifelse( as.numeric(as.character(data$variable)) >= cutoff , 1 , 0 )
  
  for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time_censor ){
      data[ i , "time" ] = time_censor
      data[ i , "surv" ] = 0
      
    }
  }
  
  if( length( data$variable[ data$variable == 1 ] )>= 5 & length( data$variable[ data$variable == 0 ] ) >= 5 ){
    
    km.coxph.plot( Surv( time, surv ) ~ variable , data.s = data, x.label = xlab, y.label = ylab, 
                   main.title = paste( title , "\n(cutoff=" , round( cutoff , 2 ) , ")" , sep="" ) ,
                   sub.title = "", 
                   leg.text = c( "Low" , "High"), 
                   leg.pos = "topright", 
                   .col = c( "#d73027","#4575b4"),  
                   show.n.risk = TRUE, 
                   n.risk.step = 7, 
                   n.risk.cex = 0.70, 
                   ylim = c(0,1), 
                   leg.inset = 0,
                   .lwd = 2 , 
                   verbose=FALSE )
   }
 
}


##########################################################################################
## volcano plot for signatures or genes association with immunotherapy responses results 
##########################################################################################

getVolcanoPlot <- function(feature, coef, pval, padj, cutoff_pos, cutoff_neg, x_lab, padj_label){
  
  data <- data.frame(feature = feature,
                     coef = coef,
                     pval = pval,
                     FDR = padj)
  
  if( padj_label == FALSE){
    
    data$diffexpressed <- "NO"
    data$diffexpressed[data$coef > 0 & data$pval < 0.05] <- "Pval < 0.05, Coef > 0"
    data$diffexpressed[data$coef < (0) & data$pval < 0.05] <- "Pval < 0.05, Coef < 0"
    
    mycolors <- c( "#fee08b","#74add1", "#bababa")
    names(mycolors) <- c("Pval < 0.05, Coef > 0", 
                         "Pval < 0.05, Coef < 0",
                         "NO")
    
    data$delabel <- NA
    data <- data[order(data$pval, decreasing = FALSE), ]
    id_pos <- data[data$coef > 0 , "feature"][1:cutoff_pos]
    id_neg <- data[data$coef < (0) , "feature"][1:cutoff_neg]
    id <- c(id_pos, id_neg)
    
    for(j in 1:length(id)){
      k <- which(data$feature == id[j])
      data$delabel[k] <- data[k, ]$feature
    }
    
  }else{
    
    data$diffexpressed <- "NO"
    data$diffexpressed[data$coef > 0 & data$FDR < 0.05] <- "FDR < 0.05, Coef > 0"
    data$diffexpressed[data$coef < (0) & data$FDR < 0.05] <- "FDR < 0.05, Coef < 0"
    
    mycolors <- c( "#fee08b","#74add1", "#bababa")
    names(mycolors) <- c("FDR < 0.05, Coef > 0", 
                         "FDR < 0.05, Coef < 0",
                         "NO")
    
    data$delabel <- NA
    data <- data[order(data$FDR, decreasing = FALSE), ]
    id_pos <- data[data$coef > 0 , "feature"][1:cutoff_pos]
    id_neg <- data[data$coef < (0) , "feature"][1:cutoff_neg]
    id <- c(id_pos, id_neg)
    
    for(j in 1:length(id)){
      k <- which(data$feature == id[j])
      data$delabel[k] <- data[k, ]$feature
    }
    
  }
 
  ggplot(data=data, aes(x=coef, y=-log10(pval), col= diffexpressed)) + 
    geom_point(size = 2.7) + theme_minimal() +
    ylab("-log10 P value") +  
    xlab(x_lab) +
    scale_colour_manual(values = mycolors) +
    theme(
      axis.text.x=element_text(size=12,  face="bold"),
      axis.title=element_text(size=12,face="bold"),
      axis.text.y=element_text(size=12, face = "bold"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      legend.position="bottom",
      legend.text = element_text(size = 9, face="bold"),
      legend.title = element_blank()) +
      geom_text_repel(aes(label= delabel),
                    size = 2.7,
                    color = "black",
                    min.segment.length = 0,
                    na.rm = TRUE, direction = "both", 
                    seed = 2356,
                    fontface= "bold",
                    max.overlaps = 50)
  
  }

