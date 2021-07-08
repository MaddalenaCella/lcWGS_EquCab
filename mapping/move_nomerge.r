#!/usr/bin/Rscript
# Author: Maddalena Cella
# Email mc2820@ic.ac.uk
# Date: 19-02-2021
# Last Modified: 19-02-2021
# Desc: script that creates a list of files that need to be moved from wgs_data/sorted to wgs_data/mapped

no_merge<- read.csv("../data/no_merge.csv", header= F, sep = ",")
list <- no_merge[,2]
write.table(list, "../data/move.list", quote=F, row.names = F, col.names=F)
