#! /bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=20:mem=240gb

module load eigensoft/6.1.4

CODEDIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq
DIR=$CODEDIR/data/gene_variants/treemix

convertf -p $DIR/PAR.PED.eigenstrat