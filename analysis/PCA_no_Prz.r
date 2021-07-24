### script to create text file with individuals to keep for second run of PCA

breeds <-read.csv("results/ancestry/clusters.csv", fill=T, stringsAsFactors = F) #dataframe with individual id and Breeds
Benson<-c("Benson","BEN","unknown")
breeds_with_Benson <- rbind(breeds, Benson)

Prz_hyb<- noquote(subset(breeds_with_Benson, CLUS == "PRZ" | CLUS == "PRZ-HY")[,1]) #left with 160 individuals

pcs <-read.table("results/ancestry/PCA_analysis.eigenvec")
list_id<- pcs[,1]
no_prz_ind<-list_id[-c(80,81,82,83,84,110,111)]
no_prz <- list_id[-c(80,81,82,83,84,110,111)]
x<- data.frame(V1=no_prz_ind, V2=no_prz)

write.table(x, file = "results/ancestry/no_Prz.txt", sep = "\t", row.names = FALSE, col.names=F, quote = F)
