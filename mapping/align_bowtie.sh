#!/bin/bash
# Author: Maddalena Cella
# Email m.cella20@imperial.ac.uk
# Date:   20-02-2021
# Last Modified: 02-07-2021
# Description: align pair-ended reads, convert to bam, index, summary stats

CODE_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/
REF_GEN=CODE_DIR/results/ref_genome/EquCab3.fna
REFBASENAME="${REF_GEN%.*}"

# create names and paths
FILE_1=$1
FILE_2=$2
BASE_NAME=$3
RES=$4

echo '-----------------------'
echo -e "\nAligning, converting bam\n"

bowtie2 -q --phred33 --very-sensitive -p 1 -I 0 -X 1500 --fr --rg-id $BASE_NAME --rg SM:$BASE_NAME --rg LB:$BASE_NAME -x $REFBASENAME -1 $FILE_1 -2 $FILE_2 -S $RES/converted/$BASE_NAME'.sam'

samtools view -bS -F 4 $RES/converted/$BASE_NAME'.sam' > $RES/converted/$BASE_NAME'.bam'
 rm -f $RES/converted/$BASE_NAME'.sam' #to save space
    

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