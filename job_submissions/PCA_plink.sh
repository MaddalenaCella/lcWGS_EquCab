#! /bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=5:mem=5gb

echo '=================================='
echo -e "\nLoading modules\n"
module load plink/1.90

HOME_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq

echo '=================================='
echo -e "\nRun PCA with plink\n"

plink --bfile $HOME_DIR/results/PCA_ready/PCA_ready --pca --horse --mind --out $HOME_DIR/results/PCA_ready/PCA_analysis