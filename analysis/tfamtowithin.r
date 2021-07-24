tfam<-read.delim("analysis/treemix.tfam", header = F) ##header set to false
breeds <-read.csv("results/ancestry/clusters.csv", fill=T, stringsAsFactors = F) #dataframe with individual id and Breeds
#fam <- read.table("results/ancestry/PCA_ready.fam")
Benson<-c("Benson",1,"unknown")
breeds_with_Benson <- rbind(breeds, Benson)

merged = cbind(tfam, breeds_with_Benson)

##remove column 3 to 8--> no need
merged<-merged[,-c(3:8)]

##change some of the names
merged$CLUSTER <- gsub("\\(>50%Quarter)","50Quarter", gsub(" ", "", merged$CLUSTER))

##write table

table <- cbind(merged$V1, merged$V2, merged$sub_group)
#table <- cbind(df_all$name,df_all$name, apprev)
write.table(table, row.names=F, sep="\t", file="results/ancestry/strat_file", quote=F, col.names=F)
