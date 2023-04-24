library(edgeR)
library(dplyr)
library(tidyverse)

cpm_normalization = function(exp_data, clinical_data,out_dir,project){
  
  raw<-read.csv(exp_data, header=TRUE, sep = ",")
  
  head(raw)
  
  if(length(which(duplicated(raw$gene_id)))>=1){
    print("Duplicated rows found and removing them")
    raw<-distinct(raw,gene_id, .keep_all = TRUE)
    print("Duplicated rows removed")
  } else{
    print("No duplicated rows found")
  }
  
  raw_mat<-raw %>% remove_rownames %>% column_to_rownames(var="gene_id")
  
  ##Read the clinical data
  clinical_data<-read.csv(clinical_data, header=TRUE, sep = ",")
  
  ##Pre-processing the metadata for edgeR##
  ## Selecting the metadata information of gene-expression columns 
  col<-data.frame(colnames(raw_mat)) 
  colnames(col)<-"icgc_specimen_id"
  metadata<-left_join(x=col,y=clinical_data, by ="icgc_specimen_id")
  
  ##Creating a simple list-based data object, called DGEList for edgeR.
  ## It will store integer read counts, with rows corresponding to genes and columns to sample/library.
  dgList<-DGEList(counts = raw_mat, genes = row.names(raw_mat), group = factor(metadata$specimen_type))
  dgList$samples
  
  ##Filter out the genes with low counts and normalize counts
  keep <- filterByExpr(y = dgList)
  dge <- dgList[keep, , keep.lib.sizes=FALSE]
  
  ##Normalize the raw read counts using TMM methods 
  dge <- calcNormFactors(object = dge,method = "TMM")
  dge$samples
  tmm <- cpm(dge)
  
  ## Save the normalized counts
  write.csv(tmm,file=paste0(out_dir,project,"_normalized_matrix.csv"))
  
}

