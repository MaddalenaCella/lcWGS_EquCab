#!/usr/bin/Rscript
# Author: Maddalena Cella
# Email mc2820@ic.ac.uk
# Date: 19-02-2021
# Last Modified: 19-02-2021
# Desc: script that creates a list of files that need to be moved from wgs_data/sorted to wgs_data/mapped

no_merge<- read_csv("~/Documents/lcWGS_EquCab/data/merge/no_merge.csv", col_names = F)
list<- length(no_merge)
for (i in no_merge[2]){
  list<- c(i)  
}
write(list, "~/Documents/lcWGS_EquCab/data/move.list")
