#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=32:mem=124gb


##define variables
#regenerated bam list
HOME=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/
LIST=$HOME/data/bam.list #check this
BASE_DIR=$HOME/job_submissions/
BAM_FILE=${PBS_ARRAY_INDEX} #bam file indexing in bam.list is the same as the current job array number
PATH_TO_FILE=$(sed -n "${BAM_FILE}p" $LIST)
name=$(echo "$PATH_TO_FILE" | cut -f 13 -d '/' | cut -f 1 -d '.')

# import unix functions
source $HOME/general/unix_functions.sh

echo '=================================='
echo -e "\nLoading modules\n"

module load samtools/1.3.1 # general
module load gatk

timer
java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator \
    -R $BASE_DIR/ref_genome/ \
    -I $PATH_TO_FILE \
    -o $BASE_DIR/wgs_data/realigned/$name"_realigner.intervals"
timer

java -jar GenomeAnalysisTK.jar -T IndelRealigner \
    -R $BASE_DIR/ref_genome/ \
    -I $PATH_TO_FILE \
    -targetIntervals $BASE_DIR/wgs_data/realigned/$name"realigner.intervals"
    -o $BASE_DIR/wgs_data/realigned/$name"_realigned.bam"

timer