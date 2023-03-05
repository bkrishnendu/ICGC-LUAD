#EdgeR normalize gene expression data and perform differential expression	

source("ICGC-LUAD_DEGs.R")

project="ICGC-LUAD"	

exp_data="F:/Bose/ICGA-LUAD/mRNA/mRNA_LUAD_raw.csv"

clinical_data="F:/Bose/ICGA-LUAD/Sample/specimen.tsv"	

out_dir="F:/Bose/ICGA-LUAD/Output/"	

calculate_deg(exp_data,clinical_data,out_dir = out_dir,project = project)
