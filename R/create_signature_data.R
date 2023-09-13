## load libraries

library(qs)
library(data.table)

#############################################
## check gene annotations + add ENS IDs
#############################################
signature <- read.csv( file= "C:/PredictIO_CodeOcean_gene/PredictIO_package/data/signature_INFO.csv")
signature$Signature <- as.character( signature$Signature )
signature$method <- as.character( signature$method )

load("C:/PredictIO_CodeOcean_gene/PredictIO_package/data/Gencode.v40.annotation.RData")
gene_annot <- features_gene
  
dir <- "C:/PredictIO_CodeOcean_gene/PredictIO_package/data/signatures/"
sig_label <- list.files(dir)[-c(12, 19, 34)]
sig_name <- substr(sig_label, 1, nchar(sig_label)-4)

for(i in 1:length(sig_label)){
  
  sig <- read.csv(paste(dir, sig_label[i], sep=""))
 
  ### GSVA 
  if(signature[signature$Signature == sig_name[i], "method"] == "GSVA"){
    
    sig <- lapply(1:nrow(sig), function(j){
      
      if( sig$gene[j] %in%  gene_annot$gene_name){
        
        data.frame(signature_name = sig$signature_name[1],
                   gene_id = gene_annot[gene_annot$gene_name == sig$gene[j], "gene_id"],
                   gene_type = gene_annot[gene_annot$gene_name == sig$gene[j], "gene_type"],
                   gene_name = gene_annot[gene_annot$gene_name == sig$gene[j], "gene_name"],
                   weight = NA)
        
      }else{
        
        data.frame(signature_name = sig$signature_name[1],
                   gene_id = NA,
                   gene_type = NA,
                   gene_name = sig$gene[j],
                   weight = NA)
      }  
      
    })
    
  }
  
  ### weighted mean
  
  if(signature[signature$Signature == sig_name[i], "method"] == "weighted_mean"){
    
    sig <- lapply(1:nrow(sig), function(j){
      
      if( sig$gene[j] %in%  gene_annot$gene_name){
        
        data.frame(signature_name = sig$signature_name[1],
                   gene_id = gene_annot[gene_annot$gene_name == sig$gene[j], "gene_id"],
                   gene_type = gene_annot[gene_annot$gene_name == sig$gene[j], "gene_type"],
                   gene_name = gene_annot[gene_annot$gene_name == sig$gene[j], "gene_name"],
                   weight = sig$coef[j])
        
      }else{
        
        data.frame(signature_name = sig$signature_name[1],
                   gene_id = NA,
                   gene_type = NA,
                   gene_name = sig$gene[j],
                   weight = sig$coef[j])
      }  
      
    })
    
  }
  
  ### PredictIO
  
  if(sig_name[i] == "PredictIO"){
    
    sig <- lapply(1:nrow(sig), function(j){
      
      if( sig$gene[j] %in%  gene_annot$gene_name){
        
        data.frame(signature_name = sig$signature_name[1],
                   gene_id = gene_annot[gene_annot$gene_name == sig$gene[j], "gene_id"],
                   gene_type = gene_annot[gene_annot$gene_name == sig$gene[j], "gene_type"],
                   gene_name = gene_annot[gene_annot$gene_name == sig$gene[j], "gene_name"],
                   weight = sig$coef[j])
        
      }else{
        
        data.frame(signature_name = sig$signature_name[1],
                   gene_id = NA,
                   gene_type = NA,
                   gene_name = sig$gene[j],
                   weight = sig$coef[j])
      }  
      
    })
    
  }
  
  
  ### COX_IS
  
  if(sig_name[i] == "COX_IS"){
    
    sig <- lapply(1:nrow(sig), function(j){
      
      if( sig$gene[j] %in%  gene_annot$gene_name){
        
        data.frame(signature_name = sig$signature_name[1],
                   gene_id = gene_annot[gene_annot$gene_name == sig$gene[j], "gene_id"],
                   gene_type = gene_annot[gene_annot$gene_name == sig$gene[j], "gene_type"],
                   gene_name = gene_annot[gene_annot$gene_name == sig$gene[j], "gene_name"],
                   weight = sig$coef[j])
        
      }else{
        
        data.frame(signature_name = sig$signature_name[1],
                   gene_id = NA,
                   gene_type = NA,
                   gene_name = sig$gene[j],
                   weight = sig$coef[j])
      }  
      
    })
    
  }
  
  ### IPS
  
  if(sig_name[i] == "IPS"){
    
    sig <- lapply(1:nrow(sig), function(j){
      
      if( sig$gene[j] %in%  gene_annot$gene_name){
        
        data.frame(signature_name = sig$signature_name[1],
                   gene_id = gene_annot[gene_annot$gene_name == sig$gene[j], "gene_id"],
                   gene_type = gene_annot[gene_annot$gene_name == sig$gene[j], "gene_type"],
                   gene_name = gene_annot[gene_annot$gene_name == sig$gene[j], "gene_name"],
                   name = sig$name[j],
                   class = sig$class[j],
                   weight = sig$coef[j])
        
      }else{
        
        data.frame(signature_name = sig$signature_name[1],
                   gene_id = NA,
                   gene_type = NA,
                   gene_name = sig$gene[j],
                   name = sig$name[j],
                   class = sig$class[j],
                   weight = sig$coef[j])
      }  
      
    })
    
  }
  
  sig <- do.call(rbind, sig)
  write.csv(sig, file = paste(paste(dir, "csv/", sep=""), 
                              paste(substr(sig_label[i], 1, nchar(sig_label[i])-4), ".csv", sep=""),
                              sep=""))
  qsave(sig, paste(paste(dir, "qs/", sep=""), 
                   paste(substr(sig_label[i], 1, nchar(sig_label[i])-4), ".qs", sep=""),
                   sep=""))
}

## IPRES (19 signatures)
dir <- "C:/PredictIO_CodeOcean_gene/PredictIO_package/data/signatures/IPRES/"
sig_label <- list.files(dir)
sig_name <- substr(sig_label, 1, nchar(sig_label)-4)

for(i in 1:length(sig_label)){
  
  sig <- read.csv(paste(dir, sig_label[i], sep=""))
  
  sig <- lapply(1:nrow(sig), function(j){
    
    if( sig$gene[j] %in%  gene_annot$gene_name){
      
      data.frame(signature_name = sig$signature_name[1],
                 sub_signature_name = sig$sub_signature_name[j],
                 gene_id = gene_annot[gene_annot$gene_name == sig$gene[j], "gene_id"],
                 gene_type = gene_annot[gene_annot$gene_name == sig$gene[j], "gene_type"],
                 gene_name = gene_annot[gene_annot$gene_name == sig$gene[j], "gene_name"],
                 weight = NA)
      
    }else{
      
      data.frame(signature_name = sig$signature_name[1],
                 sub_signature_name = sig$sub_signature_name[j],
                 gene_id = NA,
                 gene_type = NA,
                 gene_name = sig$gene[j],
                 weight = NA)
    }  
    
  })
  
  sig <- do.call(rbind, sig)
  write.csv(sig, file = paste(dir, 
                              paste(substr(sig_label[i], 1, nchar(sig_label[i])-4), ".csv", sep=""),
                              sep=""))
}


## qs file IPRES (19 signatures) 

dir <- "C:/PredictIO_CodeOcean_gene/PredictIO_package/data/signatures/IPRES/"
sig_label <- list.files(dir)
sig_name <- substr(sig_label, 1, nchar(sig_label)-4)

IPRES <- lapply(1:length(sig_name), function(i){
  
  sig <- read.csv(paste(dir, paste(substr(sig_label[i], 1, nchar(sig_label[i])-4), ".csv", sep=""),
                                sep=""))
  
  sig
  
})

names(IPRES) <- sig_name
qsave(IPRES, file= "C:/PredictIO_CodeOcean_gene/PredictIO_package/data/signatures/qs/IPRES.qs")





