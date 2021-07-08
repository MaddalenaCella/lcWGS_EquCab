##create angsd text file of genomic positions of interest
all_traits<-read.csv("../data/gene_variants/all_traits.csv")
angsd_pos<-subset(all_traits, genomic_coordinates != 'failed in Lift genome annotation')
angsd_pos<-data.frame(angsd_pos$genomic_coordinates)
write.table(angsd_pos, "../data/gene_variants/angsd_pos.txt", col.names=F, row.names = F, quote = F, na=" ") #angsd compatible text file
