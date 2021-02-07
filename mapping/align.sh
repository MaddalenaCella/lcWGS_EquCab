#!/bin/bash
# align pair-ended reads, convert to bam, index, summary stats


REF_GEN=$PBS_O_WORKDIR/ref_genome/EquCab3.fna

# create names and paths
FILE_1=$1
FILE_2=$2
BASE_NAME=$3
RES=$4

echo '-----------------------'
echo -e "\nAligning, converting bam\n"


bwa mem -t 29 $REF_GEN $FILE_1 $FILE_2 | samtools view -bS --threads 29 - > \
		$RES/converted/$BASE_NAME'.bam'


echo '-----------------------'
echo -e "\nSorting\n"

samtools sort -m 300GiB --threads 29 $RES/converted/$BASE_NAME'.bam' -o  \
		$RES/sorted/$BASE_NAME'.sorted.bam'


echo '-----------------------'
echo -e "\nIndex\n"

samtools index $RES/sorted/$BASE_NAME'.sorted.bam'

echo '-----------------------'
echo -e "\nFlagstat\n"

# should be high due to 
samtools flagstat $RES/sorted/$BASE_NAME'.sorted.bam' > \
		$RES/stats/$BASE_NAME'.stat.txt'


