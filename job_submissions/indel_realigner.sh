#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=32:mem=124gb

echo '=================================='
echo -e "\nLoading modules\n"

module load samtools/1.3.1 # general

BASEDIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/job_submissions/


echo '=================================='
echo -e "\nCalculating average depth of bam filen"

#samtools depth -aa $BASEDIR/novel_data/merged/merged_dedup_overlapclipped.bam | cut -f 3 | gzip > $BASEDIR/novel_data/merged/merged_dedup_overlapclipped.bam.depth.gz
samtools depth -a $BASEDIR/novel_data/merged/merged_dedup_overlapclipped.bam | awk '{c++;s+=$3}END{print s/c}' > $BASEDIR/novel_data/merged/merged_dedup_overlapclipped_depth.txt