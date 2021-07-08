#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=124gb

CODE_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/ #change paths as appropriate
DIR=$CODE_DIR/results/novel_data/merged/ #change path as appropriate


# import unix functions
source $HOME/general/unix_functions.sh

echo '=================================='
echo -e "\nLoading modules\n"

module load samtools/1.3.1 # general
module load java/jdk-8u144 # picard associated
module load picard/2.6.0 # cleaning
PICARD=$PICARD_HOME/picard.jar #location of Picard jar

module load anaconda3/personal
source activate /rds/general/user/mc2820/home/anaconda3/envs/horse_project
#source activate /anoconda3/envs/horse_project

cd $PBS_O_WORKDIR

echo '=================================='
echo -e "\nRemove duplicates\n"

java -Xmx60g -jar $PICARD MarkDuplicates I=$DIR/merged.bam O=$DIR/merged_dedup.bam M=$DIR/merged_dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true


echo '=================================='
echo -e "\nClip overlapping pair end reads\n"

bam clipOverlap --in $DIR/merged_dedup.bam --out $DIR/merged_dedup_overlapclipped.bam --stats

echo '=================================='
echo -e "\nIndex\n"

samtools index $DIR/merged_dedup_overlapclipped.bam $DIR/merged_dedup_overlapclipped.bam.bai

echo '=================================='
echo -e "\nAverage depth\n"

samtools depth -aa $DIR/merged_dedup_overlapclipped.bam | cut -f 3 | gzip > $DIR/merged_dedup_overlapclipped.bam.depth.gz