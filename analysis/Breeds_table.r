# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-11
# Last Modified: 2021-07-02

#packages load
require(dplyr)
library(gridExtra)

#load files needed
info <- read.csv("data/cleaned_data/info_all.csv")
runName = read.table("data/ancestry/identifiers.txt", fill=T, stringsAsFactors = F)
#the above file was generated from: for F in $(cat data/bam1.list); do     basename $F _dedup_overlapclipped.bam >> data/ancestry/identifiers.txt; done


runDf <- data.frame(index = seq(1,nrow(runName)), name = runName) #add an index from 1 to 166 to each element
names(runDf)[2]<- "name" #rename second column in the data frame 

# match the code with info_all
info_trim <- info[c("Run", "BioSample", "sub_group", "era")]

run_join <- left_join(runDf, info_trim, by=c("name" = "Run"))

run_join$sub_group <- as.character(run_join$sub_group)
run_join$sub_group[run_join$name == "final"] <- "BENSON"
run_join$era[run_join$name == "final"] <- "modern"


# strip na groups for alternative join on biosample
if(any(is.na(run_join$sub_group))){
  run_join_true <- run_join[which(!is.na(run_join$sub_group)),]
  no_group <- subset(run_join, is.na(run_join$sub_group))[,1:2]
  info_uniq <- unique(info_trim[,2:4]) # group by biosample due to multiple runs
  no_joined <- left_join(no_group, info_uniq, by=c("name" = "BioSample"))
  run_join_true <- run_join_true[-3]
  run_join <- rbind(run_join_true, no_joined)
  
}

# pop "BioSample" if exists
if("BioSample" %in% colnames(run_join)){
  run_join = run_join[which(colnames(run_join)!="BioSample")]
}

# check if other na's have passed
if(any(is.na(run_join$sub_group))){print("check bam.list and info file, NA found")}

# then final.bam
# if still no match, output error
len <- length(run_join$sub_group)

table <- cbind(run_join$index,rep(1,len), run_join$sub_group)
colnames(table) <- c("ID","FID","CLUSTER") ##maybe do something like if they belong to the same species they get a special id flag
df <- data.frame(table)
df$ID <- as.numeric(as.character(df$ID))
df <- df[order(df$ID),]

# write table out
write.table(df, row.names=F, sep="\t", file="results/ancestry/clusters", quote=F)
write_csv(df, file="results/ancestry/clusters.csv")
## Summary tables for report writeup
#a table of the individuals used and accession codes
df_info_all <- run_join %>% 
  filter(sub_group != "BENSON") %>%
  transmute(ascension_biosample = name,
            breed = sub_group
  ) %>%
  arrange(breed, ascension_biosample)

write.csv(df_info_all, "results/general/reference_panel.csv", row.names=F, quote=F)

##table with number of different breeds used
breed_table <- df_info_all %>%
  group_by(breed) %>%
  summarise(n=n())
names(breed_table)[1]<- "Breeds"#rename first column in the data frame 
names(breed_table)[2]<- "Number"

write.csv(breed_table, "results/general/breed_table.csv", row.names=F, quote=F)

pdf(file='results/general/breeds.pdf', width=7, height=12)

a<-grid.table(breed_table, rows=NULL)

dev.off()
