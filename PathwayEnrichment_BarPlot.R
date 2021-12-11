library(dplyr)
library(ggplot2)

### We have multile txt files and we want to  create a column, named "FDR". 
### This calulates and returns -log10 values of HypergeoPval columns.
### Next step, I have sorted the files based on high to low values of FDR.
### Then, I have filterd out rows which have FDR <-1.35
### Next step, I have removed the rows which have FDR = Inf.

filenames<-list.files(pattern = "*.txt$") ## loading all the txt files
filenames

for (filename in filenames) {
  type<- read.delim(filename, header = TRUE,) ## read all the txt files
  type_FDR<-mutate(type, FDR = -log10(type$HypergeoPval)) ## calulates -log10 values of HypergeoPval column.
  type_FDR<-type_FDR[,c(1:2,9,3:8)] ## rearrange the columns
  type_FDR<-arrange(type_FDR, desc(FDR)) ## sort the FDR column values from high to low
  type_FDR<-filter(type_FDR, FDR>1.35) ## filter out the rows which have FDR <- 1.35
  type_FDR<-type_FDR[!(type_FDR$FDR == Inf), ]
  write.csv(type_FDR, sub("*.txt", "-FDR.txt", filename), row.names = FALSE)
}

recurrent_up<- read.delim("Recurrent_Enrichment_up.txt", header = TRUE, sep = ",")
view(recurrent_up)
recurr_up_FDR<- mutate(recurrent_up, FDR = -log10(recurrent_up$HypergeoPval))
view(recurr_up_FDR)
recurr_up_FDR<-recurr_up_FDR[,c(1:2,9,3:8)] ## rearrange the coulumn position
recurr_up_FDR<-arrange(recurr_up_FDR, desc(FDR)) ## sort the FDR column values from high to low 
view(recurr_up_FDR)
recurr_up_FDR<-filter(recurr_up_FDR, FDR>1.35)
recurr_up_FDR<-recurr_up_FDR[!(recurr_up_FDR$FDR ==Inf), ]
view(recurr_up_FDR)
write.csv(recurr_up_FDR, file = "Recurrent_Enrichment_up-FDR.txt", row.names = FALSE)

############Barplot using ggplot2############

recurrent_up_data<-read.delim("Recurrent_Enrichment_up-FDR.txt", header = TRUE, sep = ",")

plot<-ggplot(data = recurrent_up_data, aes(x = PathAnnotation, y = FDR)) +
      geom_bar(stat = "identity", width = 0.5, fill ="steelblue")
plot + labs(x = "Pathways", y = " -log10(FDR)") + coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

primary_up_data <- read.delim("Primary_Enrichment_up-FDR.txt", header = TRUE, sep = ",")
plot<-ggplot(data = primary_up_data, aes(x = PathAnnotation, y = FDR)) +
  geom_bar(stat = "identity", width = 0.5, fill ="steelblue")
plot + labs(x = "Pathways", y = " -log10(FDR)") + coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())




