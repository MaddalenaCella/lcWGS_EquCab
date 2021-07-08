#! /bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=1:mem=96gb

module load bcftools/1.3.1
module load htslib/1.3.2
module load vcflib/2016-10-05
module load tabix/0.2.6
HOME_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/
 
echo '=================================='
echo -e "\nIndexing vcf.gz files\n"
#indexing of vcf.gz files
cd $HOME_DIR/results/Benson_imputation/
for F in chr*.vcf.gz ; do   bcftools index -f ${F}  ; done
for F in chr*.vcf.gz ; do   tabix -f -p vcf ${F}  ; done

echo '=================================='
echo -e "\nConcatenation vcf.gz files\n"
bcftools concat $HOME_DIR/results/Benson_imputation/chr*.out.vcf.gz -O z -o $HOME_DIR/results/Benson_imputation/Merged.vcf.gz

echo '=================================='
echo -e "\nIndexing merged file\n"
bcftools index $HOME_DIR/results/Benson_imputation/Merged.vcf.gz
tabix -f -p vcf $HOME_DIR/results/Benson_imputation/Merged.vcf.gz

echo '=================================='
echo -e "\nLeave just GT field\n"
bcftools annotate -x 'FORMAT','INFO' $HOME_DIR/results/Benson_imputation/Merged.vcf.gz -Oz -o $HOME_DIR/results/Benson_imputation/MergedNoFI.vcf.gz

echo '=================================='
echo -e "\nIndexing file\n"
bcftools index $HOME_DIR/results/Benson_imputation/MergedNoFI.vcf.gz
tabix -f -p vcf $HOME_DIR/results/Benson_imputation/MergedNoFI.vcf.gz

echo '=================================='
echo -e "\nIMerge Sample and Reference\n"
mkdir $HOME_DIR/results/PCA_ready
vcfcombine $HOME_DIR/results/phased_chr/RefPanelNoFI.vcf.gz $HOME_DIR/results/Benson_imputation/MergedNoFI.vcf.gz > $HOME_DIR/results/PCA_ready/combinedNoFI.vcf

##after this run the remove singletons script
