#!/bin/bash
#PBS -l walltime=35:00:00
#PBS -l select=1:ncpus=32:mem=62gb


# qsub -J 0-36 wgs_merge.sh
	# qsub -J 0-3 wgs_merge.sh
	# qsub -J 3-36 wgs_merge.sh
	# number of lines in to_merge.csv

# import unix functions
CODE_DIR=$HOME/genomics/EquSeq/

source $CODE_DIR/scripts/unix_functions.sh


#----- load modules ----#
echo '=================================='
echo -e "\nLoad samtools\n"
module load samtools/1.3.1
module load anaconda3/personal # python

#----- merge ----#
echo '=================================='
echo -e "\nCreate samtools script to run\n"

python $HOME/mapping/samtools_scripts.py $HOME/data/to_merge.csv $HOME/job_submissions/wgs_data/sorted/ $HOME/job_submissions/wgs_data/merged/ > output.txt

#python $CODE_DIR/mapping/wgs_merge.py

echo '=================================='
echo -e "\nMerge files\n"
merge=$(head -n 1 output.txt) #save the first line of the output text file as a variable 

$merge #run the variable

echo '=================================='
echo -e "\nIndex files\n"
index=$(head -n 2 output.txt) #save the second line of the output text file as a variable

$index #run the variable

rm output.txt #delete output.txt file

timer
