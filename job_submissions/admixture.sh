#! /bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=8:mem=72gb

echo '=================================='
echo -e "\nLoad modules\n"
module load plink/1.90
#installing admixture in a conda env
module load anaconda3/personal
#conda create -n admixture
#conda activate admixture
#conda install -c bioconda admixture #version admixture-1.3.0
#conda deactivate

HOME_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/

#plink --vcf $HOME_DIR/data/gene_variants/PCA_ready/final_vcf.recode.vcf --recode12 --horse --out $HOME_DIR/data/gene_variants/PCA_ready/admixture12

source activate admixture #activate conda environment

echo '=================================='
echo -e "\nPerform cross validation\n"
##perform cross validation for different values of K: A good value of K will exhibit a low cross-validation error compared to other K values before you run this
# for K in {1..31}
# do
#   admixture --cv $HOME_DIR/results/PCA_ready/PCA_ready.bed $K
# done

echo '=================================='
echo -e "\nRun admixture\n"
for K in 5, 6, 7, 8, 9, 15, 20, 25, 30
do
	cd $HOME_DIR/results/PCA_ready/; admixture $HOME_DIR/results/PCA_ready/PCA_ready.bed $K -j8
done
