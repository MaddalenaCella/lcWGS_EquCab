#!/bin/bash
#PBS -l walltime=22:00:00
#PBS -l select=1:ncpus=32:mem=62gb

# import unix functions
COD_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/
source $COD_DIR/general/unix_functions.sh

timer
#----- load modules ----#
echo '=================================='
echo -e "\nLoad samtools\n"
module load samtools/1.3.1

#----- load modules ----#
echo '=================================='
echo -e "\nCreate loimpute imput file\n"
samtools mpileup -l $COD_DIR/data/gene_variants/all.bed $COD_DIR/job_submissions/novel_data/merged/merged_dedup_overlapclipped.bam | gzip -c > $COD_DIR/job_submissions/novel_data/loimpute_input.gz

timer