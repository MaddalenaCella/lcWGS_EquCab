#!/bin/bash
#PBS -l walltime=05:00:00
#PBS -l select=1:ncpus=8:mem=10gb

# Desc: # retrieve and index reference genome

CODE_DIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/ #change it to prefererred code directory
DIR_REF=$CODE_DIR/results/ref_genome/ #change it to preferred location to store reference genome sequence

REF_TAB=assembly_report.txt 
HEADER=header.bam

echo '=================================='
echo -e "\nLoading Modules\n"
module load bwa/0.7.8
module load samtools/1.3.1 
module load anaconda3/personal
source activate horse_project


echo '=================================='
echo -e "\nDownload reference genome\n"


#-------------- EquCab3 --------------#
# specific file 
rsync --copy-links --times --verbose \
	rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_genomic.fna.gz $DIR_REF


echo '=================================='
echo -e "\nUnzipping reference genome\n"

gunzip $DIR_REF/*.gz 


echo '=================================='
echo -e "\nRename\n"

mv  $DIR_REF/GCF_002863925.1_EquCab3.0_genomic.fna $DIR_REF/EquCab3.fna


echo '=================================='
echo -e "\nGet assembly report\n"

wget -O $REF_TAB ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/all_assembly_versions/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_assembly_report.txt

echo '=================================='
echo -e "\nUpdate fasta header\n"

#echo -e "update"
python3 $CODE_DIR/scripts/updateHeader.py 'fasta' $REF_TAB $DIR_REF/EquCab3.fna

echo '=================================='
echo -e "\nClean up the environment\n"

rm $REF_TAB
mv $DIR_REF/EquCab3.fna $DIR_REF/old_EquCab3.fna 
mv $DIR_REF/EquCab3.fna.corrected.fna $DIR_REF/EquCab3.fna 

echo '=================================='
echo -e "\nIndex ref genome\n"

bwa index $DIR_REF/EquCab3.fna

samtools faidx $DIR_REF/EquCab3.fna

