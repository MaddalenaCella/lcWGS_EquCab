#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=32:mem=124gb

#define variables
HOME=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/
LIST=$HOME/data/bam.list #new bam.list with new path to file
BASE_DIR=$HOME/job_submissions/
BAM_FILE=${PBS_ARRAY_INDEX} #bam file indexing in bam.list is the same as the current job array number
PATH_TO_FILE=$(sed -n "${BAM_FILE}p" $LIST)
name=$(echo "$PATH_TO_FILE" | cut -f 13 -d '/' | cut -f 1 -d '.')

echo '=================================='
echo -e "\nLoading modules\n"
module load samtools/1.3.1
load module gatk/4.0

echo '=================================='
echo -e "\nConvert input bam to fasta for HaplotypeCaller\n"

samtools bam2fq $PATH | seqtk seq -A - > $BASE_DIR/wgs_data/realigned/$name.fa

echo '=================================='
echo -e "\nCalling Haplotypes\n"
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller \
    -R $BASE_DIR/ref_genome/EquCab.fna \
    -I $BASE_DIR/wgs_data/realigned/$name.fa \
    -o $BASE_DIR/wgs_data/calling/$name.g.vcf \
    -ERC GVFC \ ##SNELLING et al.
    --variant_index_type LINEAR \
    --variant_index_parameter 128000

## should I do joint genotyping???? probs yes 
#you need to probably do the data aggregation step and then the joint genotyping
