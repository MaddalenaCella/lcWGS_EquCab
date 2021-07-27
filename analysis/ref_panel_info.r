#!/usr/bin/Rscript
# Author: Maddalena Cella
# Email mc2820@ic.ac.uk
# Date: 01-07-2021
# Last Modified: 01-07-2021
# Desc: script that obtains info from individuals that make out the reference panel

df = read.csv("../data/cleaned_data/info_all.csv")
individuals_used = read.table("../data/identifiers.txt", fill=T, stringsAsFactors = F)
#the above file was generated from: for F in $(cat data/bam1.list); do     basename $F _dedup_overlapclipped.bam >> data/ancestry/identifiers.txt; done
x <- as.vector(t(as.matrix(individuals_used[1])))

newdf<-dplyr::filter(df, Run %in% x | BioSample %in% x & sub_group=='Przewalski') ##data frame with rows where Run or BioSample match to vector x of individuals used