#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -lselect=1:ncpus=30:mem=360gb

# ERR868003, ERR868004 - memory maxed when aligning

# qsub -J 0-626 wgsMapper.sh 
	#qsub -J 0-3 wgsMapper.sh 
	#qsub -J 3-629 wgsMapper.sh

##create directories where results of script will be stored, just for the first PBS ARRAY INDEX
if [ $PBS_ARRAY_INDEX == 0 ]
then 
	cd results
	mkdir wgs_data 
	cd wgs_data
	mkdir converted #converted files from SAM to BAM
	mkdir sorted #sorted and indexed files
	mkdir stats #summary stat
fi

cd $PBS_O_WORKDIR

CODE_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/  ##change it with your code directory
RES=$CODE_DIR/results/wgs_data/sorted/
RUNS_LIST=$CODE_DIR/data/cleaned_data/sra_runs_to_do.txt
DATA=$CODE_DIR/results/raw_files/

echo '=================================='
echo -e "\nLoading modules\n"

module load fastx/0.0.14 # trimming 
module load bwa/0.7.8 # alignment
module load samtools/1.3.1 # general
module load java/jdk-8u144 # picard associated
module load picard/2.6.0 # cleaning
module load bowtie2/2.2.9 #aligning
module load anaconda3/personal # python
source activate horse_project #for fastq-tools


# file names based on job number 
FILE=($(python3 $CODE_DIR/mapping/wgs_mapping.py $RUNS_LIST $DATA | tr -d "[''],"))

# pair ended reads
FILE_1=${FILE[0]} # 1st pair ended read
FILE_2=${FILE[1]} # 2nd pair ended read

# new shorter file name
FILE_PREFIX=${FILE[2]}

echo "read1: "$FILE_1
echo "read2: "$FILE_2
echo "prefix: "$FILE_PREFIX

timer 


echo '=================================='
echo -e "\nAligning\n"

FILE_1=$DATA/$FILE_1
FILE_2=$DATA/$FILE_2

FILE_RES=$RES/$FILE_PREFIX.sorted.bam

if [ -f "$FILE_RES" ]; then
    echo "$FILE_RES exists."
else 
    sh $CODE_DIR/mapping/align_bowtie.sh $FILE_1 $FILE_2 $FILE_PREFIX $CODE_DIR/results/wgs_data/
fi

