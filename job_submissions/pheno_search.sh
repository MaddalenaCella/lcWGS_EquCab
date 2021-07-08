#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=20gb

echo '=================================='
echo -e "\nLoading modules\n"
module load vcftools/0.1.13

##look for the presence of SNPs in Benson

HOME_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/

echo '=================================='
echo -e "\nFind SNPs\n"
mkdir $HOME_DIR/results/phenotype
vcftools --gzvcf $HOME_DIR/results/Benson_imputation/Merged.vcf.gz --recode --out $HOME_DIR/results/phenotype/pheno_pos_Benson --positions $HOME_DIR/data/gene_variants/angsd_final.file

vcftools --gzvcf $HOME_DIR/results/phased_chr/MergedRef.vcf.gz --recode --out $HOME_DIR/results/phenotype/pheno_pos_Ref --positions $HOME_DIR/data/gene_variants/angsd_final.file