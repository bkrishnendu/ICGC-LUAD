#################################################################################
## Heatmap of  interaction score of the selected edges across the subtypes ##
## The code takes a subset of output data to create the heatmaps ##
#################################################################################


library(gplots)

edge <-read.delim("Expressiondata_edges.txt", sep = "\t", header = TRUE)
View(edge)
survivingN = read.delim("EdgesSelected_minThres.txt",sep = "\t", header = TRUE)
View(survivingN)
sampleProb = read.delim("SampleClass_Probabilities.txt",sep = "\t", header = TRUE)
View(sampleProb)
subtype = read.delim("LUAD_groupinfo_example.txt",sep = "\t", header = TRUE)

normal_edge<- survivingN$Edge[which(survivingN$Normal_absdiknew>0)]
primary_edge<- survivingN$Edge[which(survivingN$Primary_absdiknew >0)]
recurrent_edge<-survivingN$Edge[which(survivingN$Recurrent_absdiknew >0)]


Xsurv<- edge[which(edge$Edge %in% survivingN$Edge),]
rownames(Xsurv) = Xsurv$Edge
Xsurv= Xsurv[, -1]
View(Xsurv)


normal_sub <- survivingN[which(survivingN$Normal_absdiknew >0), c(1:8)]
primary_sub <- survivingN[which(survivingN$Primary_absdiknew >0), c(1:4, 9:12)]
recurrent_sub <- survivingN[which(survivingN$Recurrent_absdiknew >0), c(1:4, 13:16)]

y = subtype$LUAD_Subtype[match(colnames(Xsurv), subtype$SubjectID)]
View(y)
ycol = as.numeric(factor(y))
ycol[which(ycol==1)] = "blue"
ycol[which(ycol==2)] = "green"
ycol[which(ycol==3)] = "purple"

HM = heatmap.2(as.matrix(Xsurv), trace="n", Rowv = T, Colv=T,col=bluered(20),breaks=seq(-5,5,by=0.5),ColSideColors  = ycol, labRow=F,labCol =F,mar=c(6,18), keysize=1) 


legend("right",      # location of the legend on the heatmap plot
       legend = c("Normal","Primary","Recurrent"),# category labels
       title = "Tissue Type",
       col =unique(ycol),  # color key
       lty= 1,             # line style
       lwd = 5,# line width
       cex = 0.8,
       box.lty = 0) # line type of legend box 
dev.off()

### creating the  large heatmap for selected edges based on predicted classes from iOmicsPASS output  ###

predictedClass<- sampleProb$PredictedClass[match(colnames(Xsurv), sampleProb$Subject)] 
View(predictedClass)

predictedCol <- as.numeric(factor(predictedClass)) 
predictedCol

predictedCol[which(predictedCol == 1)] = "blue"
predictedCol[which(predictedCol == 2)] = "green"
predictedCol[which(predictedCol == 3)] = "purple"
predictedCol

HM = heatmap.2(as.matrix(Xsurv), trace="n", Rowv = T, Colv=T,col=bluered(20),breaks=seq(-5,5,by=0.5),ColSideColors  = predictedCol, labRow=F,labCol =F,mar=c(6,18))

legend("right",      # location of the legend on the heatmap plot
       legend = c("Normal","Primary","Recurrent"),# category labels
       title = "Subtypes",
       col =unique(ycol),  # color key
       lty= 1,             # line style
       lwd = 5,# line width
       cex = 0.8,
       box.lty = 0)
dev.off()
### creating heatmap for selected edges separately for each subtype ###
## subsetting normal dataset from total dataset

which(predictedClass == "Normal")
NormalSurv<- Xsurv[, c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,99,126,165,180,218,233,256,288,369,371,464)]
View(NormalSurv)

NormalClass<- sampleProb$PredictedClass[match(colnames(NormalSurv), sampleProb$Subject)]
View(NormalClass)

NormalCol<-as.numeric(factor(NormalClass))
NormalCol[which(NormalCol == 1)] = "blue"
NormalCol

NormalRow<-rownames(NormalSurv)
NormalRow[survivingN$Edge[which(survivingN$Normal_absdiknew == 0)]] = "grey" 
NormalRow[survivingN$Edge[which(survivingN$Normal_absdiknew > 0)]] = "orange"
NormalRow

heatmap.2(as.matrix(NormalSurv), trace="n", Rowv = T, Colv=T,col=bluered(20),breaks=seq(-5,5,by=0.5),RowSideColors = NormalRow, ColSideColors  = NormalCol, labRow=F,labCol =F,mar=c(6,18))
dev.off()

## subsetting primary dataset from total dataset

which(predictedClass == "Primary")
PrimarySurv<-Xsurv[, c(20,21,22,23,25,26,27,28,29,31,33,34,35,36,37,38,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,61,62,63,64,65,66,68,69,71,72,73,74,75,76,77,78,80,82,83,84,85,86,87,88,89,90,91,93,94,95,96,97,98,100,101,102,104,105,108,109,110,111,112,113,114,115,117,118,119,120,121,122,123,124,125,127,128,129,130,131,132,134,136,137,138,139,140,141,142,144,145,146,147,148,149,150,151,154,155,156,157,158,159,160,161,163,164,166,167,168,169,171,172,173,174,175,176,177,178,179,181,182,183,184,185,187,189,190,191,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,219,221,222,223,224,225,226,227,228,229,231,232,234,235,238,239,240,242,243,244,245,246,247,249,250,251,254,257,258,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,277,278,279,280,282,283,284,285,286,287,289,290,291,292,294,295,296,297,298,299,300,302,303,304,305,306,307,308,310,313,314,316,317,318,319,320,321,322,323,324,325,326,327,328,329,331,332,334,335,336,337,338,339,341,342,344,345,346,347,348,349,350,351,352,353,354,355,358,359,360,361,362,363,364,365,366,367,368,370,372,373,374,375,376,377,378,379,380,381,383,384,387,388,389,390,391,392,393,394,395,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,416,417,418,421,422,423,424,425,426,427,430,431,432,434,435,436,437,438,439,440,441,443,444,445,446,447,448,450,451,453,454,455,456,457,458,459,460,461,462,463,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,483,484)]

PrimaryClass<-sampleProb$PredictedClass[match(colnames(PrimarySurv), sampleProb$Subject)]
View(PrimaryClass)

PrimaryCol<-as.numeric(factor(PrimaryClass))
PrimaryCol[which(PrimaryCol == 1)] = "green"
PrimaryCol

PrimaryRow<-rownames(PrimarySurv)
PrimaryRow[survivingN$Edge[which(survivingN$Primary_absdiknew == 0)]] = "grey" 
PrimaryRow[survivingN$Edge[which(survivingN$Primary_absdiknew > 0)]] = "orange"
PrimaryRow

heatmap.2(as.matrix(PrimarySurv), trace="n", Rowv = T, Colv=T,col=bluered(20),breaks=seq(-5,5,by=0.5),RowSideColors = PrimaryRow, ColSideColors  = PrimaryCol, labRow=F,labCol =F,mar=c(6,18))
dev.off()

## subsetting recurrent dataset from total dataset

which(predictedClass == "Recurrent")

RecurrentSurv <- Xsurv[, c(24,30,32,39,60,67,70,79,81,92,103,106,107,116,133,135,143,152,153,162,170,186,188,192,220,230,236,237,241,248,252,253,255,259,276,281,293,301,309,311,312,315,330,333,340,343,356,357,382,385,386,396,415,419,420,428,429,433,442,449,452,482,485,486)]

View(RecurrentSurv)

RecurrentClass<-sampleProb$PredictedClass[match(colnames(RecurrentSurv), sampleProb$Subject)]

RecurrentCol<-as.numeric(factor(RecurrentClass))
RecurrentCol[which(RecurrentCol == 1)] ="purple"
RecurrentCol

RecurrentRow<-rownames(RecurrentSurv)
RecurrentRow[survivingN$Edge[which(survivingN$Recurrent_absdiknew == 0)]] = "grey" 
RecurrentRow[survivingN$Edge[which(survivingN$Recurrent_absdiknew > 0)]] = "orange"
RecurrentRow

heatmap.2(as.matrix(RecurrentSurv), trace="n", Rowv = T, Colv=T,col=bluered(20),breaks=seq(-5,5,by=0.5),RowSideColors = RecurrentRow, ColSideColors  = RecurrentCol, labRow=F,labCol =F,mar=c(6,18))
dev.off()
