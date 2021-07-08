#!/bin/bash
#PBS -l walltime=35:00:00
#PBS -l select=1:ncpus=32:mem=62gb

##qsub -J 1-172 ## wc -l number of lines in the newly generted bam1.list

CODE_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/
WGS_DATA=$CODE_DIR/results/wgs_data/ #change paths as appropriate
DIR=$WGS_DATA/processed/ #change path as appropriate
LIST=$CODE_DIR/data/bam1.list #to create this list remember to run the makebam.sh file!
BAM_FILE=${PBS_ARRAY_INDEX} #bam file indexing in bam.list is the same as the current job array number

if [ $PBS_ARRAY_INDEX == 0 ]
then
    	cd $WGS_DATA                      
        mkdir processed #merged files
        mkdir depth
fi


echo '=================================='
echo -e "\nLoading modules\n"

module load samtools/1.3.1 # general
module load java/jdk-8u144 # picard associated
module load picard/2.6.0 # cleaning
module load anaconda3/personal
source activate /rds/general/user/mc2820/home/anaconda3/envs/horse_project #bamutil is downloaded here
PICARD=$PICARD_HOME/picard.jar #where Picard jarfile is located

#path to file to perform filtering, marking duplicates, clip overlapping pairs, indexing and calculating average depth

PATH_TO_FILE=$(sed -n "${BAM_FILE}p" $LIST)
name=$(echo "$PATH_TO_FILE" | cut -f 13 -d '/' | cut -f 1 -d '.') ##probably need to change this

echo '=================================='
echo -e "\nFilter poorly mapped reads\n"

samtools view -h -q 10 $PATH_TO_FILE | samtools view -buS - | samtools sort -o $DIR$name'_minq10.bam'

echo '=================================='
echo -e "\nRemove duplicates\n"

java -Xmx60g -jar $PICARD MarkDuplicates I=$DIR$name'_minq10.bam' O=$DIR$name'_dedup.bam' M=$WGS_DATA/stats/$name'_merged_dupstat.txt' VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

##remove minq10 files to free memory
rm $DIR$name'_minq10.bam'

echo '=================================='
echo -e "\nClip overlapping pair end reads\n"

bam clipOverlap --in $DIR$name'_dedup.bam' --out $DIR$name'_dedup_overlapclipped.bam' --stats

timer
##remove dedup.bam files to free memory
rm $DIR$name'_dedup.bam'

echo '=================================='
echo -e "\nIndexing\n" ##should i do this????

samtools index $DIR$name'_dedup_overlapclipped.bam' $DIR$name'_dedup_overlapclipped.bam.bai'


echo '=================================='
echo -e "\nCalculating average depth of bam filen"

#samtools depth -aa $BASEDIR/novel_data/merged/merged_dedup_overlapclipped.bam | cut -f 3 | gzip > $BASEDIR/novel_data/merged/merged_dedup_overlapclipped.bam.depth.gz
samtools depth -a $DIR$name'_dedup_overlapclipped.bam' | awk '{c++;s+=$3}END{print s/c}' > $WGS_DATA/depth/$name'_dedup_overlapclipped_depth.txt'
