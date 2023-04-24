
# ICGC-LUAD

The project aims to normalize gene expression data and perform differential expression analysis from ICGC-Lung adenocarcinoma (LUAD) project.

1. Normalize gene expression data and perform differential expression using edgeR: 

        
        source("ICGC-LUAD_DEGs.R") # Call the function
        
        project="ICGC-LUAD" # project name
        
        exp_data="F:/Bose/ICGA-LUAD/mRNA/mRNA_LUAD_raw.csv" # Raw read count directory
        
        clinical_data="F:/Bose/ICGA-LUAD/Sample/specimen.tsv" # Metadata directory
       
        out_dir="F:/Bose/ICGA-LUAD/Output/"  # Output directory
    
        cpm_normalization(exp_data,clinical_data,out_dir = out_dir,project = project) # perform analysis
