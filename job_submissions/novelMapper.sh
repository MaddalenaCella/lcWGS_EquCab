#!/bin/bash
#PBS -l walltime=40:00:00
#PBS -lselect=1:ncpus=32:mem=62gb


# qsub -J 0-24 novelMapper.sh

##create directories where results of script will be stored, just for the first PBS ARRAY INDEX
CODE_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/ #change it to prefererred code directory

if [ $PBS_ARRAY_INDEX == 0 ]
then 
	cd $CODE_DIR/results
	mkdir novel_data 
	cd novel_data
	mkdir raw_files #renamed raw files
	mkdir trimmed #trimmed files
	mkdir converted #converted files from SAM to BAM
	mkdir sorted #sorted and indexed files
	mkdir stats #summary stat
fi

cd $PBS_O_WORKDIR

DATA_DIR=$CODE_DIR/data/Benson/
RES=$CODE_DIR/results/novel_data/

echo '=================================='
echo -e "\nLoading modules\n"
module load fastx/0.0.14 # trimming 
module load bwa/0.7.8 # alignment
module load samtools/1.3.1 # general
module load java/jdk-8u144 # picard associated
module load picard/2.6.0 # cleaning
module load anaconda3/personal # python
module load htslib


# file names based on job number 
FILE=($(python3 $CODE_DIR/mapping/names.py | tr -d "[''],"))

# pair ended reads
READ1=${FILE[0]} # 1st pair ended read
READ2=${FILE[1]} # 2nd pair ended read

# new shorter file name
FILE_PREFIX=${FILE[2]}

echo "read1: "$READ1
echo "read2: "$READ2
echo "prefix: "$FILE_PREFIX


echo '=================================='
echo -e "\nCopy reads from DATA_DIR\n"

# copy read files to where results will be stored
# read 1
cp $DATA_DIR/Clean/F19FTSEUHT1854-swab-horse-1A/$READ1 \
        $RES/raw_files/$FILE_PREFIX'_1.fq.gz'

# read 2
cp $DATA_DIR/Clean/F19FTSEUHT1854-swab-horse-1A/$READ2 \
        $RES/raw_files/$FILE_PREFIX'_2.fq.gz'


timer 
echo '=================================='
echo -e "\nUnzipping\n"

gunzip $RES/raw_files/$FILE_PREFIX'_1.fq.gz'
gunzip $RES/raw_files/$FILE_PREFIX'_2.fq.gz'

timer 
echo '=================================='
echo -e "\nTrimming\n"

FILE_1=$RES/trimmed/$FILE_PREFIX'_1.trim.fq'
FILE_2=$RES/trimmed/$FILE_PREFIX'_2.trim.fq'

echo 'read 1'
fastx_trimmer -l 90 -i $RES/raw_files/$FILE_PREFIX'_1.fq' -o $FILE_1
		# ~ 8 mins per read
echo 'read 2'
fastx_trimmer -l 90 -i $RES/raw_files/$FILE_PREFIX'_2.fq' -o $FILE_2

echo '=================================='
echo -e "\nMapping\n"


sh $CODE_DIR/mapping/align.sh $FILE_1 $FILE_2 $FILE_PREFIX $RES
