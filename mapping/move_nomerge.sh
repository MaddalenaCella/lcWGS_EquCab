#! /bin/bash
# Author: Maddalena Cella
# Email mc2820@ic.ac.uk
# Date: 19-02-2021
# Last Modified: 19-02-2021
# Desc: move files that did not require merging to merged folder

CODE_DIR=rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/
OUT_DIR=$CODE_DIR/job_submissions/wgs_data/merged/

echo '=================================='
echo -e "\nCreate move list\n"
module load anaconda3/personal #make sure R is installed 
FILES=($(cat $CODE_DIR/data/move.list))

echo '=================================='
echo -e "\nCreate move list\n"
Rscript --vanilla CODE_DIR/mapping/move_nomerge.r 


for file in "${FILES[*]}"
do
	#echo "moving " $file " to " $OUT_DIR "\n"
	mv $CODE_DIR/job_submissions/wgs_data/sorted/$file.bam $OUT_DIR
    mv $CODE_DIR/job_submissions/wgs_data/sorted/$file.bam.bai $OUT_DIR
done