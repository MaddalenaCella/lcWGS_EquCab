library(dplyr)
library(Hmisc) #for plot of SD
library(formattable) ##to make nice tables
library(knitr)

# Export a Formattable as PNG, PDF, or JPEG
# check: https://github.com/renkun-ken/formattable/issues/26
export_formattable 

##tables
f3_all_combos <- read.table("results/ancestry/f3stat_poplist_filtered.txt")

f3_all_combos_unique <- 
  f3_all_combos[!duplicated(f3_all_combos$V7), ]

f3_close_breeds <- read.table("results/ancestry/eigenstrat/f3stat_filtered.txt")
f3_close_breeds<- subset(f3_close_breeds, select = -V1)
col_names<- c("Pop A", "Pop B", "Target", "F3", "std_err", "Z", "SNPs")
names(f3_close_breeds)<- col_names ##dataframe with column names

ordered_f3<- f3_close_breeds[order(f3_close_breeds$F3, decreasing = F),]

table_f3_try<- ordered_f3[1:10,]
row.names(table_f3_try) <- NULL #to remove row numbering from table
formattable(table_f3_try, align =c("c","c","c","c","c", "c", "c"))

pdf("results/ancestry/eigenstrat/errorbars_f3")
errbar(1:10, ordered_f3$F3[1:10],
       (ordered_f3$F3+ordered_f3$std_err)[1:10],
       (ordered_f3$F3-ordered_f3$std_err)[1:10], pch=20, las=2, cex.axis=0.4, xaxt='n',
       xlab="Population", ylab="F3")
#axis(1, at=1:10, labels=ordered_f3$`Source 2`[1:10], las=2, cex.axis=0.6)
dev.off()


##get the poplist and outgroup txt files produced by f3 stat analysis

f3_poplist <- read.table("results/ancestry/eigenstrat/f3_poplist_filtered.txt")
f3_poplist<- subset(f3_poplist, select = -V1)
col_names<- c("Pop A", "Pop B", "Target", "F3", "std_err", "Z", "SNPs")
names(f3_poplist)<- col_names ##dataframe with column names

f3_poplist_ordered<- f3_poplist[order(f3_poplist$F3, decreasing = F),]

table_f3<- f3_poplist_ordered[1:10,]
row.names(table_f3) <- NULL #to remove row numbering from table
formattable(table_f3, align =c("c","c","c","c","c", "c", "c"))

#outgroup
f3_outgroup <- read.table("results/ancestry/eigenstrat/f3_outgroup_filtered.txt")
f3_outgroup<- subset(f3_outgroup, select = -V1)
col_names_outgroup<- c("Target", "Pop B", "Outgroup", "F3", "std_err", "Z", "SNPs")
names(f3_outgroup)<- col_names_outgroup ##dataframe with column names

f3_outgroup_ordered<- f3_outgroup[order(f3_outgroup$F3, decreasing = T),] ##In this scenario, the statistic F3(A, B; C) measures the branch length from C to the common ancestor of A and B, coloured red. 
#So this statistic is simply a measure of how closely two population A and B are related with each other, as measured from a distant outgroup. 
#It is thus a similarity measure: The higher the statistic, the more genetically similar A and B are to one another.

table_f3_out<- f3_outgroup_ordered[1:10,]
row.names(table_f3_out) <- NULL #to remove row numbering from table
formattable(table_f3_out, align =c("c","c","c","c","c", "c", "c"))

