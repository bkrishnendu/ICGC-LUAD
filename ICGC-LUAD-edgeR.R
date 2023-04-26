library(edgeR)
library(dplyr)
library(tidyverse)

##Defining the data directory
data_dir<-"F:/Bose/ICGA-LUAD/mRNA/"
##Defining the metadata directory
meta_dir<-"F:/Bose/ICGA-LUAD/Sample/"
out_dir = "F:/Bose/ICGA-LUAD/Output/"

##Pre-processing the data frame for edgeR
raw<-read.csv(paste0(data_dir,"mRNA_LUAD_raw.csv"), header=TRUE, sep = ",")
head(raw)

if(length(which(duplicated(raw$gene_id)))>=1){
  print("Duplicated rows found and removing them")
  raw<-distinct(raw,gene_id, .keep_all = TRUE)
  print("Duplicated rows removed")
} else{
  print("No duplicated rows found")
}


raw_mat<-raw %>% remove_rownames %>% column_to_rownames(var="gene_id")


specimen<-read.table(paste0(meta_dir,"specimen.tsv"), header = TRUE, sep = ",")
specimen<-specimen[,c(1,4,5,7)]
donor = read.table(paste0(meta_dir, "donor.tsv"), header = TRUE, sep = ",")
##Prepare the complete clinical_data
clinical_data<-left_join(x=specimen,y=donor, by ="icgc_donor_id")
write.csv(clinical_data, file=paste0(meta_dir,"LUAD_Clinical_Data.csv"))

##Pre-processing the metadata for edgeR##
## Selecting the metadata information of gene-expression columns 
col<-data.frame(colnames(raw_mat)) 
colnames(col)<-"icgc_specimen_id"
metadata<-left_join(x=col,y=clinical_data, by ="icgc_specimen_id")

## Remove "Recurrent tumor samples"
which(metadata$specimen_type=="Recurrent tumour - solid tissue")
raw_mat=raw_mat[,c(-231,-289)]
metadata=metadata[c(-231,-289),]



##Creating a simple list-based data object, called DGEList for edgeR.
## It will store integer read counts, with rows corresponding to genes and columns to sample/library.
dgList<-DGEList(counts = raw_mat, genes = row.names(raw_mat), group = factor(metadata$specimen_type))
dgList$samples

##Filter out the genes with low counts and normalize counts
keep <- filterByExpr(y = dgList)
dge <- dgList[keep, , keep.lib.sizes=FALSE]

dge <- calcNormFactors(object = dge)
dge$samples
dge <- estimateDisp(y = dge)
plotBCV(dge)



## Design the model matrix
design <- model.matrix(~ dge$samples$group)
design
rownames(design) <- colnames(dge)

## DEG Analysis
fit <- glmFit(dge, design)
fit
colnames(fit)
lrt<-glmLRT(fit)
deg= as.data.frame(topTags(lrt, n="Inf"))
deg=deg[deg$FDR <0.05,]

write.csv(deg, file= paste0(out_dir,"ICGC-LUAD_DEGs.csv"))






