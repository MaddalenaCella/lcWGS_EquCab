#! /bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=40:mem=480gb

## run as qsub -J 1-8 evalAdmix.sh --> one job per number of .Q files you want to evaluate
EVAL_ADMIX=/rds/general/project/human-popgen-datasets/live/Maddalena/evalAdmix ##where the executable is located
CODE_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq ## change as appropriate


if [ $PBS_ARRAY_INDEX == 1 ]
then
	NUM=5
elif [ $PBS_ARRAY_INDEX == 2 ]
then
	NUM=6
elif [ $PBS_ARRAY_INDEX == 3 ]
then
	NUM=7
elif [ $PBS_ARRAY_INDEX == 4 ]
then
    	NUM=8
elif [ $PBS_ARRAY_INDEX == 5 ]
then
    	NUM=9
elif [ $PBS_ARRAY_INDEX == 6 ]
then
    	NUM=15
elif [ $PBS_ARRAY_INDEX == 7 ]
then
    	NUM=25
elif [ $PBS_ARRAY_INDEX == 8 ]
then
    	NUM=30
fi

$EVAL_ADMIX/evalAdmix -plink $CODE_DIR/data/gene_variants/PCA_ready/PCA_ready \
-fname $CODE_DIR/data/gene_variants/PCA_ready/PCA_ready.$NUM.P -qname $CODE_DIR/data/gene_variants/PCA_ready/PCA_ready.$NUM.Q \
-o $CODE_DIR/data/gene_variants/PCA_ready/evaladmix$NUM.corr -autosomeMax 31 -P 30 
