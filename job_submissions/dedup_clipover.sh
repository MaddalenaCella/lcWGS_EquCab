#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=124gb

HOME=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/ #change paths as appropriate
DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/job_submissions/novel_data/merged/ #change path as appropriate


# import unix functions
source $HOME/general/unix_functions.sh

echo '=================================='
echo -e "\nLoading modules\n"

module load samtools/1.3.1 # general
module load java/jdk-8u144 # picard associated
module load picard/2.6.0 # cleaning
module load bamutil

PICARD=$PICARD_HOME/picard.jar

timer 

echo '=================================='
echo -e "\nFilter poorly mapped reads\n"

samtools view -h -q 10 $DIR/merged.bam | samtools view -buS - | samtools sort -o $DIR/merged_minq10.bam #filter out MAPQ smaller than 10

echo '=================================='
echo -e "\nRemove duplicates\n"

java -Xmx60g -jar $PICARD MarkDuplicates I=$DIR/merged.bam O=$DIR/merged_dedup.bam M=$DIR/merged_dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

timer 

echo '=================================='
echo -e "\nClip overlapping pair end reads\n"

bam clipOverlap --in $DIR/merged_dedup.bam --out $DIR/merged_dedup_overlapclipped.bam --stats

timer



