#!/bin/bash
#PBS -lwalltime=48:00:00
#PBS -lselect=1:ncpus=1:mem=1gb

# Desc: Get files from ncbi sra
	# wc -l 
	# qsub -J 1-630 fastqDownload.sh
		# qsub -J 0-3 fastqDownload.sh
		# qsub -J 3-629 fastqDownload.sh

CODE_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/ #insert your code directory path here (pwd command)
RES_DIR=$CODE_DIR/results/raw_files #change it to desired results directory

echo '=================================='
echo -e "\nLoad modules\n"
module load anaconda3/personal

echo '=================================='
echo -e "\nGet Fasta files\n"

python3 $CODE_DIR/ref/getFastq.py $CODE_DIR/data/cleaned_data/sra_runs_to_do.txt $RES_DIR