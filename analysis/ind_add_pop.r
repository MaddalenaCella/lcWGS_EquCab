#!/usr/bin/Rscript
# Author: Maddalena Cella
# Email mc2820@ic.ac.uk
# Date: 15-07-2021
# Last Modified: 15-07-2021
# Desc: script that adds pop info to treemix.ind file for f3 statistics

ind_eigen<- read_delim("results/ancestry/eigenstrat/treemix.ind", col_names=F, " ")
breeds <-read.csv("results/ancestry/clusters.csv", fill=T, stringsAsFactors = F)
Benson<-c("Benson","BEN","unknown")
breeds_with_benson <- rbind(breeds, Benson)
ind_eigen_pop = cbind(ind_eigen, breeds_with_benson)

drop <- c("ID", "X3", "CLUSTER")
to_write<- ind_eigen_pop[,!(names(ind_eigen_pop) %in% drop)] ##drop unnecessary columns defined in the drop vector

write.table(to_write, "results/ancestry/eigenstrat/treemixpop.ind", col.names=F, sep=" ", row.names=F, quote=F)

