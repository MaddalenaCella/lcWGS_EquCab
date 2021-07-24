#!/usr/bin/Rscript
# Author: Maddalena Cella
# Email mc2820@ic.ac.uk
# Date: 15-07-2021
# Last Modified: 15-07-2021
# Desc: script that creates all possible combinations of populations for f3 statistics

##import breeds table
breeds <-read.csv("results/ancestry/clusters.csv", fill=T, stringsAsFactors = F)

## generate all possibile combinations of populations
vec_breeds<-breeds$CLUS ##get the column of interest and convert it to a vector
breeds_combo<-t(combn(vec_breeds,2)) ##combn generates all possibile combinations
Ben<- rep("BEN", nrow(breeds_combo)) ## the third column should be just Benson (the focal individual)
breeds_combo_B<- cbind(breeds_combo, Ben) ##combine into  single dataframe

##export as tsv
write.table(breeds_combo_B, "results/ancestry/eigenstrat/pop_combos_all.tsv", col.names=F, sep=" ", row.names=F, quote=F)
