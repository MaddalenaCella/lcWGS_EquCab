#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=32:mem=62gb

HOME_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq
module load java
module load vcftools
module load htslib
module load bcftools

ANGSD=/rds/general/project/human-popgen-datasets/live/Maddalena/angsd/angsd
NUM=`echo $PBS_ARRAY_INDEX`
REF=$HOME_DIR/results/ref_genome/EquCab3.fna

echo '=================================='
echo -e "\nExtract variants for each chromosome in the reference panel\n" 
vcftools --gzvcf $HOME_DIR/results/phased_chr/chr$NUM.out.gt.vcf.gz --out $HOME_DIR/results/variants_chr/chr$NUM.filtered --site-quality

##generate SNP sites list in a format compatible with ANGSD
awk 'NF{NF-=1};1' <$HOME_DIR/results/variants_chr/chr$NUM.filtered.lqual >$HOME_DIR/results/variants_chr/chr$NUM.sites
sed '1d' $HOME_DIR/results/variants_chr/chr$NUM.sites > $HOME_DIR/results/variants_chr/chr$NUM.txt

echo '=================================='
echo -e "\nCall SNPs in Benson\n" 

$ANGSD sites index $HOME_DIR/results/variants_chr/chr$NUM.txt

$ANGSD -i $HOME_DIR/results/novel_data/merged/merged_dedup_overlapclipped.bam -ref $REF -P 7 \
        -out $HOME_DIR/results/beagle_anc_input/gl.benson.chr$NUM -sites $HOME_DIR/results/variants_chr/chr$NUM.txt -r chr$NUM: \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
        -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
        -checkBamHeaders 0 \
        -dovcf 1 -GL 1 -doGlf 2 -doPost 1 -doMajorMinor 1 -SNP_pval 1e-1 -doMaf 1 \
        -doGeno 8 -dumpCounts 2 -doDepth 1 -doCounts 1

echo '=================================='
echo -e "\nUnzip chr file\n" 
bgzip -d $HOME_DIR/results/phased_chr/chr$FILE_LINE.out.gt.vcf.gz

echo '=================================='
echo -e "\nRename sample file\n" 
bcftools index -f $HOME_DIR/results/beagle_anc_input/gl.benson.chr$NUM.vcf.gz

mkdir $HOME_DIR/data/sample_name
nano $HOME_DIR/data/sample_name/name.txt #write preferred name for the sample, in my case BENSON
bcftools reheader --samples $HOME_DIR/data/sample_name/name.txt $HOME_DIR/results/beagle_anc_input/gl.benson.chr$NUM.vcf.gz -o $HOME_DIR/results/beagle_anc_input/gl.benson.reheaded.chr$NUM.vcf.gz

echo '=================================='
echo -e "\nImputation with Beagle\n" 
mkdir $HOME_DIR/results/Benson_imputation/
java -Xmx96g -jar $HOME_DIR/beagle.r1399.jar ref=$HOME_DIR/results/phased_chr/chr$NUM.out.gt.vcf gt=$HOME_DIR/results/beagle_anc_input/gl.benson.reheaded.chr$NUM.vcf.gz out=$HOME_DIR/results/Benson_imputation/chr$NUM.out window=100000000

echo '=================================='
echo -e "\nZip back chr file\n" 
bgzip $HOME_DIR/results/phased_chr/chr$FILE_LINE.out.gt.vcf
