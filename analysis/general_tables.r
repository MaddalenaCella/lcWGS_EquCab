#!/usr/bin/Rscript
# Author: Maddalena Cella
# Email mc2820@ic.ac.uk
# Date: 24-07-2021
# Last Modified: 24-07-2021
# Description: script that plots summary tables for reference panel information

#packages load
require(dplyr)
library(gridExtra)

breeds <-read.csv("results/ancestry/clusters.csv", fill=T, stringsAsFactors = F)
breeds_SRcodes<- read.csv("results/general/clusters_SRcodes.csv", fill=T, stringsAsFactors = F)


breeds_noid <- subset(breeds, select = -ID) #remove ID column
a<-breeds_noid$CLUS
a ##find a way to add this vector to column

breed_table <- breeds_noid %>%
  group_by(CLUSTER, CLUS) %>%
  summarise(n=n())
names(breed_table)[1]<- "Breeds" #rename first column in the data frame 
names(breed_table)[2]<- "Breeds code" 
names(breed_table)[3]<- "Number"


pdf(file='results/general/breeds_code.pdf', width=7, height=12)
a<-grid.table(breed_table, rows=NULL)
dev.off()

breeds_SRcodes_noid <- subset(breeds_SRcodes, select = -c(index,era))

  accession_table<- breeds_SRcodes_noid %>%
  arrange(sub_group)
names(accession_table)[2]<- "Breeds" #rename first column in the data frame 
names(accession_table)[1]<- "Accession" 

pdf(file='results/general/breeds_accession.pdf', width=7, height=48)
b<-grid.table(accession_table, rows=NULL)
dev.off()
