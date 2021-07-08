#!/bin/bash
#PBS -l walltime=35:00:00
#PBS -l select=1:ncpus=32:mem=62gb


# qsub -J 0-36 wgs_merge_2.sh
	# qsub -J 0-3 wgs_merge_2.sh
	# qsub -J 3-36 wgs_merge_2.sh
	# number of lines in to_merge.csv


CODE_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/


if [ $PBS_ARRAY_INDEX == 1 ]
then 
	cd $CODE_DIR/results/wgs_data/
	mkdir merged #merged files
fi

cd $PBS_O_WORKDIR

#----- load modules ----#
echo '=================================='
echo -e "\nLoad samtools\n"
module load samtools/1.3.1
module load anaconda3/personal # python

#----- merge ----#
echo '=================================='
echo -e "\nCreate samtools script to run\n"

python $CODE_DIR/mapping/samtools_scripts.py $CODE_DIR/data/to_merge.csv $CODE_DIR/results/wgs_data/sorted/ $CODE_DIR/results/wgs_data/merged/ > output.txt


echo '=================================='
echo -e "\nMerging and indexing files\n"

sh output.txt

rm output.txt #delete output.txt file

