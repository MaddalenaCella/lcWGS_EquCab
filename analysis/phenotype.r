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

pos_found <-all_traits[all_traits$genomic_coordinates %in% vcf_pheno$coordinates, ]
pos_found_imp<- subset(pos_found, select = -c(type.of.variant, coordinates_paper, allele, ECA.number, Position.based.on.reference., Year.Published, Genomic.Coordinates, Major, Minor, coding.position, protein.position, genomic.indexing, Breed, Notes, Notes2))

all_info <- cbind(pos_found_imp, Major=vcf_pheno$REF, Minor= vcf_pheno$ALT, Format= vcf_pheno$FORMAT, Horse_sample= vcf_pheno$Benson)
all_info <- all_info[,c(1,2,3,5,6,7,8,4)]
colnames(all_info)[3] <- "Genomic.Coordinates"

pdf(file='results/phenotype/postable.pdf', width=15, height=12)
a<-grid.table(all_info, rows=NULL)
dev.off()

