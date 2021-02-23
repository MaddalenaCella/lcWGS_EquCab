#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=124gb

# if $PBS_ARRAY_INDEX==0
# do
# cd $PBS_O_WORKDIR/wga_data/
# mkdir processed
# fi

HOME=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/
WGS_DATA=$HOME/job_submissions/wgs_data/ #change paths as appropriate
DIR=$WGS_DATA/processed/ #change path as appropriate
LIST=$HOME/data/bam1.list #to create this list remember to run the makebam.sh file!
BAM_FILE=${PBS_ARRAY_INDEX} #bam file indexing in bam.list is the same as the current job array number

# import unix functions
source $HOME/general/unix_functions.sh

echo '=================================='
echo -e "\nLoading modules\n"

module load samtools/1.3.1 # general
module load java/jdk-8u144 # picard associated
module load picard/2.6.0 # cleaning
module load anaconda3/personal
source activate /rds/general/user/mc2820/home/anaconda3/envs/horse_project #bamutil is downloaded here

PICARD=$PICARD_HOME/picard.jar #not sure about this:comapre to the one used for novel_data

timer 

#path to file to perform filtering, marking duplicates, clip overlapping pairs, indexing and calculating average depth

PATH=$(sed -n "${BAM_FILE}p" $LIST)
name=$(echo "$PATH" | cut -f 13 -d '/' | cut -f 1 -d '.')
echo '=================================='
echo -e "\nFilter poorly mapped reads\n"

samtools view -h -q 10 $PATH | samtools view -buS - | samtools sort -o $DIR$name'_minq10.bam' #filter out MAPQ smaller than 10

echo '=================================='
echo -e "\nRemove duplicates\n"
java -Xmx60g -jar $PICARD MarkDuplicates I=$DIR$name'_minq10.bam' O=$DIR$name'_dedup.bam' M=$WGS_DATA$name'_merged_dupstat.txt' VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

echo '=================================='
echo -e "\nClip overlapping pair end reads\n"
bam clipOverlap --in $DIR$name'_dedup.bam' --out $DIR$name'_dedup_overlapclipped.bam' --stats

echo '=================================='
echo -e "\nIndexing\n" ##should i do this????
samtools index $DIR$name'_dedup_overlapclipped.bam' $DIR$name'_dedup_overlapclipped.bam.bai'

echo '=================================='
echo -e "\nCalculating average depth of bam filen"

#samtools depth -aa $BASEDIR/novel_data/merged/merged_dedup_overlapclipped.bam | cut -f 3 | gzip > $BASEDIR/novel_data/merged/merged_dedup_overlapclipped.bam.depth.gz
samtools depth -a $DIR$name'_dedup_overlapclipped.bam' | awk '{c++;s+=$3}END{print s/c}' > $HOME/job_submissions/wgs_data/merged/depth/$name'_dedup_overlapclipped_depth.txt'


# while IFS= read -r line
# do
#     name=$(echo "$line" | cut -f 13 -d '/' | cut -f 1 -d '.') #have a look at this

#     echo '=================================='
#     echo -e "\nFilter poorly mapped reads\n"

#     samtools view -h -q 10 $line | samtools view -buS - | samtools sort -o $DIR$name'_minq10.bam' #filter out MAPQ smaller than 10

#     echo '=================================='
#     echo -e "\nRemove duplicates\n"
#     java -Xmx60g -jar $PICARD MarkDuplicates I=$DIR$name'_minq10.bam' O=$DIR$name'_dedup.bam' M=$DIR$name'merged_dupstat.txt' VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

#     echo '=================================='
#     echo -e "\nClip overlapping pair end reads\n"
#     bam clipOverlap --in $DIR$name'_dedup.bam' --out $DIR$name'_dedup_overlapclipped.bam' --stats

#     echo '=================================='
#     echo -e "\nIndexing\n" ##should i do thisd????
#     samtools index $DIR$name'_dedup_overlapclipped.bam' $DIR$name'_dedup_overlapclipped.bam.bai'

#     echo '=================================='
#     echo -e "\nCalculating average depth of bam filen"

# #samtools depth -aa $BASEDIR/novel_data/merged/merged_dedup_overlapclipped.bam | cut -f 3 | gzip > $BASEDIR/novel_data/merged/merged_dedup_overlapclipped.bam.depth.gz
#     samtools depth -a $DIR$name'_dedup_overlapclipped.bam' | awk '{c++;s+=$3}END{print s/c}' > $HOME/job_submissions/wgs_data/merged/depth/$name'_dedup_overlapclipped_depth.txt'

# done < "$LIST"

timer



