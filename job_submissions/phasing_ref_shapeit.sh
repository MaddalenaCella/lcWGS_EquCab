#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=16:ncpus=32:mem=124gb

#qsub -J 1-31 phasing_ref_shapeit.sh


echo '=================================='
echo -e "\nLoad modules\n"
module load java
module load shapeit/2.778

CODE_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/


if [ $PBS_ARRAY_INDEX == 1 ]
then 
	cd $CODE_DIR/results/
	mkdir phased_chr #phased chromosomes
fi

FILE_LINE=`echo $PBS_ARRAY_INDEX`

echo '=================================='
echo -e "\nPhasing variants in each chromosome\n"

shapeit -V $CODE_DIR/results/SNP_calls_ref/chr$FILE_LINE:.refpanel.anc..vcf.gz -O $CODE_DIR/results/phased_chr/chr$FILE_LINE.phased

echo '=================================='
echo -e "\nConverting haps into vcf format\n"
shapeit -convert --input-haps $CODE_DIR/results/phased_chr/chr$FILE_LINE.phased \
--output-vcf $CODE_DIR/results/phased_chr/chr$FILE_LINE.out.gt.vcf

#remove file in haps format
rm $CODE_DIR/results/phased_chr/chr$FILE_LINE.phased