#QTLs from https://www.animalgenome.org/cgi-bin/QTLdb/EC/index
size<- read.delim("~/Documents/lcWGS_EquCab/data/gene_variants/qtls/size.txt", header=T)
fertility<- read.delim("~/Documents/lcWGS_EquCab/data/gene_variants/qtls/fertility.txt", header=T)
gaits<- read.delim("~/Documents/lcWGS_EquCab/data/gene_variants/qtls/gaits.txt", header=T)
racing<- read.delim("~/Documents/lcWGS_EquCab/data/gene_variants/qtls/racing_ability.txt", header=T)
disease_injury<- read.delim("~/Documents/lcWGS_EquCab/data/gene_variants/qtls/disease_injury.txt", header=T)

##keep only those that are significant and for which there is information on the genomic position
sign_subset <- function(x){
  subset(x, x$Significance=='Significant' & x$Chromosome != 'NA' & x$Cood_A_bp != 'NA')
}

size_s<- sign_subset(size) #612
fertility_s<- sign_subset(fertility) #50
gaits_s<- sign_subset(gaits) #220
racing_s <- sign_subset(racing) #38
disease_s<- sign_subset(disease_injury) #940

##merge all these QTLs into a signle dataframe
qtls<- merge(merge(merge(merge(size_s, fertility_s, all = T), racing_s, all=T), gaits_s, all = T), disease_s, all=T) #tot 1860
qtls<- qtls[order(qtls$Chromosome),] # sort by ascending chromosome number

##create bed file of the format chr start end snp_id
bed <- qtls[,c('Chromosome', 'Cood_A_bp', 'Coord_B_bp', 'QTL_ID')]
colnames(bed) <- c('chrom', 'chromStart', 'chromEnd', 'name')
write.table(bed, "~/Documents/lcWGS_EquCab/data/gene_variants/qtls/QTL.bed", row.names = F)

##save qtls as csv 
write.csv(qtls, "~/Documents/lcWGS_EquCab/data/gene_variants/qtls/qtls_interest.csv")

#Mendelian variants from https://onlinelibrary.wiley.com/doi/epdf/10.1111/age.12857
pigmentation<- read.csv("~/Documents/lcWGS_EquCab/data/gene_variants/mendelian/genetic_variants_pigmentation.csv", header=TRUE, skip=1, na.strings=c(""," ","NA"))
disease<- read.csv("~/Documents/lcWGS_EquCab/data/gene_variants/mendelian/disease_performance_variants.csv", header=TRUE, skip=1, na.strings=c(""," ","NA"))
read
mend_subset <- function(x){
  subset(x, x$ECA.number != 'NA' & x$coordinate != 'NA')
}
keep=c('regulatory', 'missense', 'nonsense (stop-gain)', 'splicing') #type of variants to keep
pigmentation_s<- mend_subset(pigmentation)
pigmentation_s<- pigmentation_s[pigmentation_s$type.of.variant %in% keep, ] #remove indels from list (hard to map)

disease_s_m<- subset(disease, Genomic.Coordinate != 'NA')# & type.of.variant != "regulatory" & type.of.variant != "nonsense (stop-gain)")
disease_s_m<- disease_s_m[disease_s_m$type.of.variant %in% keep, ] #remove indels from list (hard to map)

#create a new column with just coordinate position
pigmentation_s$chromStart<- pigmentation_s$coordinate

library(tidyverse)
pigmentation_s$chromPos<-str_remove_all(pigmentation_s$chromStart, "[^0-9_ins]")
#str_remove_all(pigmentation_s$chromPos, "[^0-9_ins]")
pigmentation_s$chromStart<-str_extract(pigmentation_s$chromPos, "[0-9]{8}")
pigmentation_s$chromStart[c(27,33)]<- c('122833887', '141677402')
pigmentation_s$chromeEnd<- pigmentation_s$chromStart
pigmentation_s$chromPos<-NULL
pigmentation_s$trial<-NULL
pigmentation_s$ID<- 1:nrow(pigmentation_s)

pigm_bed <- pigmentation_s[,c('ECA.number', 'chromStart', 'chromeEnd', 'ID')]
colnames(pigm_bed) <- c('chrom', 'chromStart', 'chromEnd', 'name')

disease_s_m$chromPos<-disease_s_m$Genomic.Coordinate
disease_s_m$chromPos<-str_remove_all(disease_s_m$chromPos, "[^0-9_ins]")
disease_s_m$chromStart<-str_extract(disease_s_m$chromPos, "[0-9]{8}")
disease_s_m$chromStart[c(2,9)]<- c("128056148",'4535550')
disease_s_m$chromeEnd<- disease_s_m$chromStart
disease_s_m$ID<- 1:nrow(disease_s_m)+20
disease_s_m$chromPos<-NULL

dis_bed <- disease_s_m[,c('ECA.Number', 'chromStart', 'chromeEnd', 'ID')]
colnames(dis_bed) <- c('chrom', 'chromStart', 'chromEnd', 'name')

mend_bed<- merge(pigm_bed, dis_bed, all = T)
write.table(mend_bed, "~/Documents/lcWGS_EquCab/data/gene_variants/mendelian/mend.bed", row.names = F)
write.csv(disease_s_m, "~/Documents/lcWGS_EquCab/data/gene_variants/mendelian/mend_disease.csv")
write.csv(pigmentation_s, "~/Documents/lcWGS_EquCab/data/gene_variants/mendelian/mend_pigm.csv")


#merge qtl and mendelian bed 
all_bed<- merge(mend_bed, bed, all=T)
write.table(all_bed, "~/Documents/lcWGS_EquCab/data/gene_variants/all.bed", col.names=F, row.names = F, quote = F)
