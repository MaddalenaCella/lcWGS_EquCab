#!/usr/bin/Rscript
# Author: Maddalena Cella
# Email mc2820@ic.ac.uk
# Date: 15-07-2021
# Last Modified: 15-07-2021
# Desc:

#libraries
library(data.table) #good for reading vcf files
library(dplyr)
library(gridExtra)

vcf_pheno <- fread("results/phenotype/pheno_pos_Benson.recode.vcf")
all_traits<-read.csv("data/gene_variants/all_traits.csv")

vcf_pheno$coordinates <- gsub(" ", "", paste(vcf_pheno$`#CHROM`, ":", vcf_pheno$POS, "-", vcf_pheno$POS)) # add a column to vcf_pheno with coordinates in the format og=f chr:start-end
vcf_pheno<- vcf_pheno[order(`#CHROM`),]
pos_found <-all_traits[all_traits$genomic_coordinates %in% vcf_pheno$coordinates, ]
pos_found_imp<- subset(pos_found, select = -c(type.of.variant, coordinates_paper, allele, ECA.number, Position.based.on.reference., Year.Published, Genomic.Coordinates, coding.position, protein.position, genomic.indexing, Breed, Notes, Notes2))
pos_found_v2 <- subset(pos_found, select = -c(type.of.variant, coordinates_paper, allele, ECA.number, Position.based.on.reference., Year.Published, Genomic.Coordinates, coding.position, protein.position, genomic.indexing, Breed, Notes, Notes2, Major, Minor))

all_info <- cbind(pos_found_imp, Format= vcf_pheno$FORMAT, Horse_sample= vcf_pheno$Benson)
all_info_v2 <- cbind(pos_found_v2, Format= vcf_pheno$FORMAT, Horse_sample= vcf_pheno$Benson, Major=vcf_pheno$REF, Minor=vcf_pheno$ALT)
all_info <- all_info[,c(1,2,3,4,5,6,7)]
colnames(all_info)[4] <- "EquCab3.0 coordinates"

write.csv(all_info_v2, 'results/phenotype/info_phen_2.csv')
info_modified <- read.csv("results/phenotype/info_phen_2.csv")
all_info <- cbind(all_info, `Genotype Benson`= info_modified$Benson.genotype, `Genotype probability`=info_modified$Genotype.probability)
all_info <- all_info[,c(1,2,3,4,5,6,8,9,7)]
colnames(all_info)[2] <- "Trait type"

### pdf with what was found in Benson
pdf(file='results/phenotype/postable.pdf', width=17, height=9)
a<-grid.table(all_info, rows=NULL)
dev.off()


### pdf with all Mendelian traits

all_mendelian <- subset(all_traits, all_traits$Trait.Inheritance == "Mendelian" | all_traits$Trait.Inheritance== "Mendelian recessive" |
                          all_traits$Trait.Inheritance== "Mendelian dominant" )
all_mendelian <- all_mendelian[,c(1,2,3,4,5,6,9,10,11,12,15,16,17)]
colnames(all_mendelian)[7] <- "Original coordinates"
colnames(all_mendelian)[8] <- "EquCab3.0 coordinates"
all_mendelian <- sapply(all_mendelian, as.character)
all_mendelian[is.na(all_mendelian)] <- " " ##replace blanks with NAs

pdf(file='results/phenotype/all_mendelian.pdf', width=26, height=16)
b<-grid.table(all_mendelian, rows=NULL) 
dev.off()


### pdf wih all QTLs
all_QTLs <- subset(all_traits, all_traits$Trait.Inheritance == "QTL" )
all_QTLs <- all_QTLs[,c(1,2,3,4,5,6,9,10,11,12,15,16,17)]
colnames(all_QTLs)[7] <- "Original coordinates"
colnames(all_QTLs)[8] <- "EquCab3.0 coordinates"
all_QTLs <- sapply(all_QTLs, as.character) 
all_QTLs[is.na(all_QTLs)] <- " " ##replace NAs with blanks

pdf(file='results/phenotype/all_QTLs.pdf', width=26, height=19)
c<-grid.table(all_QTLs, rows=NULL) 
dev.off()

