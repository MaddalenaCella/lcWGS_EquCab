#! /bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=8:mem=96gb

#qsub -J 1-31 #one run per chromosome
echo '=================================='
echo -e "\nLoading modules\n"
ANGSD=/rds/general/project/human-popgen-datasets/live/Maddalena/angsd/angsd #directory where angsd is installed

CODE_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/

if [ $PBS_ARRAY_INDEX == 1 ]
then 
	cd $CODE_DIR/results/
	mkdir snp_calls_ref #SNP calls per chromosome
fi

OUT_DIR=$CODE_DIR/results/snp_calls_ref/

FILE_LINE=`echo $PBS_ARRAY_INDEX`
RF_INPUT=`sed -n "$FILE_LINE"p $CHROM` #get chromosome number to call SNPs by chromosome

# reference genome
REF=$CODE_DIR/results/ref_genome/EquCab3.fna

echo '=================================='
echo -e "\nSNP calling\n"

$ANGSD -b $CODE_DIR/data/bam1.list -ref $REF -P 7 \
        -out $OUT_DIR/$RF_INPUT.refpanel.anc. -r $RF_INPUT \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
        -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
        -checkBamHeaders 0 \
        -dovcf 1 -GL 1 -doGlf 2 -doPost 1 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 \
        -doGeno 8 -dumpCounts 2 -doDepth 1 -doCounts 1
