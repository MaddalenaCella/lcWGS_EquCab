#!/usr/bin/Rscript
# Author: Maddalena Cella
# Email mc2820@ic.ac.uk
# Date: 24-07-2021
# Last Modified: 24-07-2021
# Description: Scripts that summaries the output of the f3 stat from qp3Pop (AdmixTools),


library(dplyr)
library(Hmisc) #for plot of SD
library(knitr)
library(gridExtra) ##look into this
library("grid")
library("ggplotify")

##get the poplist and outgroup txt files produced by f3 stat analysis

f3_poplist <- read.table("results/ancestry/eigenstrat/f3_poplist_filtered.txt")
f3_poplist<- subset(f3_poplist, select = -c(V1,V8))
col_names<- c("Pop A", "Pop B", "Target", "F3", "SE", "Z-score")
names(f3_poplist)<- col_names ##data frame with column names

f3_poplist_ordered<- f3_poplist[order(f3_poplist$F3, decreasing = F),]

table_f3<- f3_poplist_ordered[1:10,]
row.names(table_f3) <- NULL #to remove row numbering from table

pdf(file='results/ancestry/eigenstrat/poplist_table.pdf', width=5, height=3.3)
b<-grid.table(table_f3, rows=NULL)
dev.off()

pdf("results/ancestry/eigenstrat/errorbars_f3_poplist.pdf")
plot1<-errbar(1:10, table_f3$F3[1:10],
       (table_f3$F3+table_f3$SE)[1:10],
       (table_f3$F3-table_f3$SE)[1:10], pch=20, las=2, cex.axis=0.4, xaxt='n',
       xlab="Population A (above) and B (below)", ylab="F3")+ 
  axis(1, at=1:10, labels=table_f3$`Pop B`[1:10], las=2, cex.axis=0.6) +
  axis(3, at=1:10, labels=table_f3$`Pop A`[1:10], las=2, cex.axis=0.6)
dev.off()

table_f3$Populations <- paste(table_f3$`Pop A`, "-", table_f3$`Pop B`, "-", table_f3$Target)
write.csv(table_f3, "results/ancestry/eigenstrat/modified_f3.csv", row.names = F)

#outgroup
f3_outgroup <- read.table("results/ancestry/eigenstrat/f3_outgroup_filtered.txt")
f3_outgroup<- subset(f3_outgroup, select = -c(V1,V8))
col_names_outgroup<- c("Target", "Pop B", "Outgroup", "F3", "SE", "Z-score")
names(f3_outgroup)<- col_names_outgroup ##dataframe with column names

f3_outgroup_ordered<- f3_outgroup[order(f3_outgroup$F3, decreasing = T),] ##In this scenario, the statistic F3(A, B; C) measures the branch length from C to the common ancestor of A and B, coloured red. 
#So this statistic is simply a measure of how closely two population A and B are related with each other, as measured from a distant outgroup. 
#It is thus a similarity measure: The higher the statistic, the more genetically similar A and B are to one another.

table_f3_out<- f3_outgroup_ordered[1:10,]
row.names(table_f3_out) <- NULL #to remove row numbering from table

pdf(file='results/ancestry/eigenstrat/outgroup_table.pdf', width=5, height=3.3)
a<-grid.table(table_f3_out, rows=NULL)
dev.off()


pdf("results/ancestry/eigenstrat/errorbars_f3_outgroup.pdf")
errbar(1:10, table_f3_out$F3[1:10],
       (table_f3_out$F3+table_f3_out$SE)[1:10],
       (table_f3_out$F3-table_f3_out$SE)[1:10], pch=20, las=2, cex.axis=0.4, xaxt='n',
       xlab="Population B", ylab="F3")
axis(1, at=1:10, labels=table_f3_out$`Pop B`[1:10], las=2, cex.axis=0.6)
dev.off()

table_f3_out$Populations <- paste(table_f3_out$`Target`, "-", table_f3_out$`Pop B`, "-", table_f3_out$Outgroup)
write.csv(table_f3_out, "results/ancestry/eigenstrat/modified_f3_outgroup.csv", row.names = F)


f4_stat <- read.table("results/ancestry/eigenstrat/f4_results_filtered.txt")
f4_stat<- subset(f4_stat, select = -c(V1,V8,V9,V10))
col_names<- c("Outgroup", "Pop A", "Pop B", "Target", "F4", "Z-score")
names(f4_stat)<- col_names ##data frame with column names

f4_stat_ordered<- f4_stat[order(f4_stat$F4, decreasing = T),]

row.names(table_f3) <- NULL #to remove row numbering from table

pdf(file='results/ancestry/eigenstrat/poplist_table.pdf', width=5, height=3.3)
b<-grid.table(table_f3, rows=NULL)
dev.off()

