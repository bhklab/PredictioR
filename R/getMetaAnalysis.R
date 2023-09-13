## load libraries

library(qs)
library(meta)
library(metafor)
library(forestplot)
library(ggplot2)
library(ggrepel)

##############################################################
##############################################################
## create a function to specify cancer specific analyses
##############################################################
##############################################################

getCancer <- function( cancer ){
  
  if( length( grep( "Coef" , colnames( cancer ) ) ) ){
    
    cancer$Coef <- as.numeric( as.character( cancer$Coef ) )
    
  }
  
  cancer$Study <- as.character( cancer$Study )
  cancer$Cancer_type <- as.character( cancer$Cancer_type )
  cancer$Pval <- as.numeric(as.character( cancer$Pval )) 
  cancer$SE <- as.numeric(as.character( cancer$SE )) 
  
  tab <- table( cancer$Cancer_type)[ table(cancer$Cancer_type) %in% c(1,2) ]
  
  cancer$Study <- paste( cancer$Study , ", n = " , cancer$N , sep= "" ) 
  cancer$Cancer_type[ cancer$Cancer_type %in% names(tab) ] <- "Other"
  
  cancer
  
}

##############################################################
##############################################################
## create meta-analysis function
##############################################################
##############################################################

getMeta <- function(Coef, SE, Study, Pval, N, cancer_dat, gene){
  
  data <- data.frame( Gene = gene,
                      Study = as.character( Study ),
                      N = N, 
                      Coef = as.numeric(as.character( Coef )),
                      SE = as.numeric(as.character( SE )),
                      Pval = as.numeric(as.character( Pval )),
                      Cancer_type = as.character( do.call(rbind, strsplit(Study, split='__', fixed=TRUE))[, 2] ))
  
  data <- data[ order( data$Coef ) , ]
  
  ## at least 3 studies needed to do the random effect meta-analyses
  if(nrow(data) < 3){ stop("not enough studies to do meta-analysis") }else{
    
    if( cancer_dat == TRUE ){
      
      cancer <- getCancer( data )
      data <- cancer[ order( cancer$Coef ) , ]
      
    }
    
    meta <- metagen( TE = Coef,
                     seTE = SE ,
                     data = data ,
                     studlab = Study ,
                     fixed = FALSE ,
                     random = TRUE ,
                     control = list( maxiter = 10000 , stepadj=0.5 ) )
    
    meta_res <- data.frame(Gene = gene,
                           Coef = meta$TE.random ,  
                           SE = meta$seTE.random ,   
                           CI_lower = meta$lower.random ,   
                           CI_upper = meta$upper.random ,
                           Pval = meta$pval.random , 
                           I2 = meta$I2 ,
                           Q_Pval = meta$pval.Q )
    
    return(list(input_data = data,
                meta_output = meta, 
                meta_summery = meta_res))
    
      }
  }

##############################################################
##############################################################
## create forestplot function
##############################################################
##############################################################

getForestplotPanCancer <- function(Coef, SE, Study, Pval, N , xlab , label, cancer_dat, gene){
 
  res <- getMeta(Coef, SE, Study, Pval, N, cancer_dat = FALSE, gene)
  data <- res$input_data
  meta <- res$meta_output
  
  m <- c( min( c( 0 , data$Coef ) , na.rm=TRUE) - .5 , ( max( c( 0 , abs(data$Coef) ) , na.rm=TRUE) ) + .5 )
  data$Study <- paste( data$Study , ", n = " , data$N , sep= "" ) 
  
  forest( meta , 
          leftcols = c("studlab", "effect.ci" , "Pval" ),
          leftlabs= c( "Study" , paste(label, "[95%CI]", sep = " ") , "P-value" ) , 
          xlab = xlab, 
          digits.se = 2 ,
          colgap.forest=unit(10, "mm") ,
          plotwidth = unit( 30 , "mm") , 
          pooled.totals = TRUE,
          smlab = " ",
          comb.random =TRUE,
          comb.fixed = FALSE,
          text.fixed.w = FALSE,
          layout = "JAMA",
          print.I2.ci = TRUE,
          print.Q = FALSE,
          print.pval.Q = TRUE,
          print.I2 = TRUE,
          print.tau2 = FALSE,
          resid.hetstat = FALSE,
          test.overall.random = TRUE,
          test.overall.fixed = FALSE,
          xlim = m ,  
          col.square= "black" ,  
          col.study= "black" ,  
          col.square.lines = "black" ,
          col.diamond.random  = "#1565c0"  ,
          col.diamond.lines.random  ="#1565c0" ,
          col.by = "#1565c0",
          addrow.subgroups=TRUE )
  
 }


#########################################################
## Per cancer
#########################################################

getForestplotPerCancer <- function( Coef, SE, Study, Pval, N , gene, cancer_dat, xlab , label){
  
  res <- getMeta(Coef, SE, Study, Pval, N, cancer_dat = TRUE, gene)
  cancer <- getCancer(res$input_data)
  
  remove <- names( table( cancer$Cancer_type)[ table(cancer$Cancer_type) %in% c(1,2) ] )
  
  if( length( unique( cancer$Cancer_type[ !cancer$Cancer_type %in% remove ] ) ) > 1 ){
    
    m <- c( min( c( 0 , cancer$Coef ) , na.rm=TRUE) - .5 , ( max( c( 0 , abs(cancer$Coef) ) , na.rm=TRUE) ) + .5 )
    meta <- res$meta_output
    
    if( length(remove) > 0 ){
      
      meta.subgroup <- update.meta(meta , 
                                   byvar= Cancer_type , 
                                   exclude = cancer$cancer_type %in% remove ,
                                   fixed = FALSE ,
                                   random = TRUE ,
                                   control = list( maxiter = 10000 , stepadj=0.5 ) )
    } else{
      meta.subgroup <- update.meta(meta , 
                                   byvar=Cancer_type , 
                                   comb.random = TRUE ,
                                   fixed = FALSE ,
                                   random = TRUE ,
                                   control = list( maxiter = 10000 , stepadj=0.5 ) )
    }
   
  }
  
  ## needs to be improved !!!!!!
  forest( meta.subgroup, 
          #leftcols = c("studlab", "effect.ci", "Pval"),
          #leftlabs= c( "Study" , paste(label, "[95%CI]", sep = " ") , "P-value" ) , 
          digits.se = 2,
          colgap.forest=unit(10, "mm") ,
          plotwidth = unit( 30 , "mm") , 
          xlab = xlab,
          pooled.totals = TRUE,
          smlab = " ",
          comb.random =TRUE,
          comb.fixed = FALSE,
          text.fixed.w = FALSE,
          layout = "JAMA",
          print.I2.ci = TRUE,
          print.Q = FALSE,
          print.pval.Q = TRUE,
          print.I2 = TRUE,
          print.tau2 = FALSE,
          resid.hetstat = FALSE,
          test.overall.random = TRUE,
          test.overall.fixed = FALSE,
          xlim = m , 
          col.square= "black" ,  
          col.study= "black" ,  
          col.square.lines = "black" ,
          col.diamond.random  = "#1565c0"  ,
          col.diamond.lines.random  ="#1565c0" ,
          col.by = "#1565c0",
          addrow.subgroups=TRUE )
  
}


