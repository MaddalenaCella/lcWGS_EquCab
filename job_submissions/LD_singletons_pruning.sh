#! /bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=5:mem=5gb

module load vcftools
module load plink/1.90
HOME_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/


echo '=================================='
echo -e "\nRemove singletons\n"
vcftools --vcf $HOME_DIR/results/PCA_ready/combinedNoFI.vcf --singletons --out $HOME_DIR/results/PCA_ready/combinedNoFI

vcftools --vcf $HOME_DIR/results/PCA_ready/combinedNoFI.vcf \
--exclude-positions $HOME_DIR/results/PCA_ready/combinedNoFI.singletons \
--recode --recode-INFO-all --out $HOME_DIR/results/PCA_ready/final_vcf

echo '=================================='
echo -e "\nRemove LD\n"
plink --vcf $HOME_DIR/results/PCA_ready/final_vcf.recode.vcf --maf 0.01 --indep-pairwise 50 5 0.2 --horse --out $HOME_DIR/results/PCA_ready/final-noLD 
plink --vcf $HOME_DIR/results/PCA_ready/final_vcf.recode.vcf --extract $HOME_DIR/results/PCA_ready/final-noLD.prune.in --make-bed --horse --out $HOME_DIR/results/PCA_ready/PCA_ready
