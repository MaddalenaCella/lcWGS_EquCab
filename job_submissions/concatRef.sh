#! /bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=1:mem=96gb

echo '=================================='
echo -e "\nLoad modules\n"

module load bcftools/1.3.1
module load htslib/1.3.2
module load vcflib
module load tabix

HOME_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/
PHASED_CHR=$HOME_DIR/results/phased_chr/

##REFERENCE PANEL PROCESSING
#bgzip $HOME_DIR/data/gene_variants/phased_chr/chr*.out.gt.vcf #if vcf file is not .gz, then uncomment this line of code 

echo '=================================='
echo -e "\nIndexing vcf.gz files\n"
cd $PHASED_CHR
for F in chr*.out.gt.vcf.gz ; do   bcftools index -f ${F}  ; done
for F in chr*.out.gt.vcf.gz ; do   tabix -f -p vcf ${F}  ; done

echo '=================================='
echo -e "\nConcatenating vcf.gz files\n"
bcftools concat $PHASED_CHR/chr*.out.gt.vcf.gz -O z -o $PHASED_CHR/RefPanel.vcf.gz

echo '=================================='
echo -e "\nIndexing concatenated file\n"
cd $PHASED_CHR
bcftools index -f $PHASED_CHR/RefPanel.vcf.gz
tabix -f -p vcf  $PHASED_CHR/RefPanel.vcf.gz

echo '=================================='
echo -e "\nKeep only gt field in vcf file\n"
bcftools annotate -x 'FORMAT','INFO' $PHASED_CHR/RefPanel.vcf.gz -o $PHASED_CHR/RefPanelNoFI.vcf

echo '=================================='
echo -e "\nUnzip\n"
bgzip $PHASED_CHR/RefPanelNoFI.vcf

echo '=================================='
echo -e "\nIndexing concatenated file\n"
bcftools index -f $PHASED_CHR/RefPanelNoFI.vcf.gz
tabix -f -p vcf $PHASED_CHR/RefPanelNoFI.vcf.gz