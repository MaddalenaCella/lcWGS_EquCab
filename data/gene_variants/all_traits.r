#!/usr/bin/Rscript
# Author: Maddalena Cella
# Email mc2820@ic.ac.uk
# Date: 15-03-2021
# Last Modified: 15-07-2021
# Description: creates angsd text file of genomic positions of interest

all_traits<-read.csv("../data/gene_variants/all_traits.csv")
angsd_pos<-subset(all_traits, genomic_coordinates != 'failed in Lift genome annotation')
angsd_pos<-data.frame(angsd_pos$genomic_coordinates)
write.table(angsd_pos, "../data/gene_variants/angsd_pos.txt", col.names=F, row.names = F, quote = F, na=" ") #angsd compatible text file
