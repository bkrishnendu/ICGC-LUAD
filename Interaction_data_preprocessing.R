## condisering the interactions for each type with edge score over 0.1, cosidered as upregulated interactions.
## condisering the interactions for each type with edge score less than -0.1, cosidered as downregulated interactions.
## row binding for the combined  data frame.
## deleting reversed duplicates.  

interaction<-read.delim("EdgesSelected_minThres.txt", header = TRUE, sep = "\t") ## loading the total interaction data
interaction<-interaction[, c(-1,-6,-7,-8,-10,-11,-12,-14,-15,-16)] ## deleting the first coulmn

interaction[,c(1,2)]<-apply(interaction[,c(1,2)], 2, function(y) gsub("_mrna|_prot", "", y))
View(interaction)
class(interaction$InteractionType)
interaction$InteractionType<-as.character(interaction$InteractionType)


recurrent_up<-interaction[(interaction$Recurrent_dikNew >= 0.1 & interaction$Primary_dikNew <=0), ] 
View(recurrent_up)
recurrent_down<-interaction[(interaction$Recurrent_dikNew <= -0.1 & interaction$Primary_dikNew >= 0), ]
View(recurrent_down)
recurrent_up_down<-rbind(recurrent_up,recurrent_down)
View(recurrent_up_down)
recurrent_up_down<-recurrent_up_down[!duplicated(data.frame(t(apply(recurrent_up_down,1, sort)))),] ## Deleting reversed duplicates
View(recurrent_up_down)
write.csv(recurrent_up_down, file = "recurrent_up_down.txt", row.names = F)

primary_up<-interaction[(interaction$Primary_dikNew>= 0.1 & interaction$Recurrent_dikNew <= 0), ]
View(primary_up)
primary_down<-interaction[(interaction$Primary_dikNew<= -0.1 & interaction$Recurrent_dikNew >= 0), ]
View(primary_down)
primary_up_down<- rbind(primary_up, primary_down)
View(primary_up_down)
primary_up_down<-primary_up_down[!duplicated(data.frame(t(apply(primary_up_down,1, sort)))),] ## Deleting reversed duplicates
View(primary_up_down)
write.csv(primary_up_down, file = "primary_up_down.txt", row.names = F)

normal_up<-interaction[(interaction$Normal_dikNew >= 0.1 & interaction$Primary_dikNew<=0 & interaction$Recurrent_dikNew<= 0), ]
View(normal_up)
normal_down<-interaction[(interaction$Normal_dikNew <= -0.1 & interaction$Primary_dikNew>=0 & interaction$Recurrent_dikNew>= 0), ]
View(normal_down)
normal_up_down<-rbind(normal_up, normal_down)
View(normal_up_down)
normal_up_down<-normal_up_down[!duplicated(data.frame(t(apply(normal_up_down,1, sort)))),] ## Deleting reversed duplicates
View(normal_up_down)
write.csv(normal_up_down, file = "normal_up_down.txt", row.names = F)
