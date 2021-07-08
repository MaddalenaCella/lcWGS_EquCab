#!/bin/bash
#PBS -l walltime=22:00:00
#PBS -l select=1:ncpus=32:mem=62gb

#create directory where to store merged files
CODE_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/ #change it to prefererred code directory

cd $CODE_DIR/results/novel_data/
mkdir merged

cd $PBS_O_WORKDIR

DIR=$CODE_DIR/results/novel_data/
DATA=($DIR/sorted/*.bam)
DATA_FILTERED=($DIR/sorted/*_minq10_sorted.bam)
FILE_OUT=$DIR/merged/merged.bam

#----- load modules ----#
echo '=================================='
echo -e "\nLoad samtools\n"
module load samtools/1.3.1

#----- filter poorly mapped reads----#
echo '=================================='
echo -e "\nRemove reads with quality lower than 10\n"

for file in $DIR/sorted/*.bam; do
    samtools view -h -q 10 $file | samtools view -buS - | samtools sort -o $file'_minq10_sorted.bam'
done

#----- merge ----#
echo '=================================='
echo -e "\nMerge Files\n"
samtools merge -f --threads 31 $FILE_OUT $DATA_FILTERED


echo '=================================='
echo -e "\nIndex\n"

samtools index $FILE_OUT $FILE_OUT.bai


echo '=================================='
echo -e "\nflagstat\n"

samtools flagstat $FILE_OUT > $DIR/merged'.stat1.txt'

## clean up memory space
cd $DIR
rm -r converted
