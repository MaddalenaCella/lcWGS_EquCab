#! /bin/bash
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=5:mem=5gb
# Desc: prepare treemix files on the hpc

module load anaconda3/personal
module load vcftools
module load plink
module load bcftools

CODEDIR=/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq
DIR=$CODEDIR/data/gene_variants/treemix
FILE=$CODEDIR/data/gene_variants/PCA_ready/final_vcf.recode.vcf
FILE_NEW=$DIR/treemix


################################################
### ALL OF THE 	ABOVE CAN BE RUN IN THE TERMINAL, SUB THE REST BELOW AS A JOB
################################################
# below is modified of: https://speciationgenomics.github.io/Treemix/


echo '=================================='
echo "\nmaking treemix file\n"

# snps don't appear in map ped with plink, first use vcftools 
	# and apply some filtering

# make map and ped
	# allows for reformatting
echo "vcftools"

vcftools --vcf $FILE --plink-tped --mac 2 \
	--remove-indels --max-alleles 2 --out $FILE_NEW #--positions snp.list

plink --tped $FILE_NEW.tped --tfam $FILE_NEW.tfam --horse --recode --out $FILE_NEW

echo "correcting map"
# snps cause issues if no name is present - rename them
awk -F"\t" '{
        split($2,chr,":")
	$1="1"
	$2="1:"chr[2]
        print $0
}' ${FILE_NEW}.map > better.map
mv better.map ${FILE_NEW}.map

echo "generate stratified freq file"
# generate stratified freq file 
	# and filter 
	# --mind 0.1 # important for SE but benson is absent, instead remove sites manually later
plink --file $FILE_NEW --snps-only --geno 0.1 \
	--maf 0.02 --freq --missing \
	--family --out $FILE_NEW #--mind 0.1 

echo "zipping"
gzip -f $FILE_NEW".frq.strat"

#zcat < $file".frq.strat.gz" | head

echo "plink2treemix.py"
# convert using treemix python script
python2 $EPHEMERAL/dependencies/plink2treemix.py $FILE_NEW".frq.strat.gz" $DIR/treemix.frq.gz
#zcat treemix.frq.gz | head
#zcat treemix.frq.gz | grep BENSON
#zcat treemix.frq.gz | head -n1 | tr " " "\n" | wc -l

### remove SNPs that aren't present with our sample 
	# issues arise with treemix algorithms if lots of our sample is absent

#### KEEP FOR NOW - remove snps not associated with BENSON
	# run as parallel analysis
#python $CODEDIR/ancestry/bensonSNPs.py

echo "snp count:"
zcat $DIR/treemix.frq.gz | wc -l

echo '=================================='
echo "\nprepare cluster file for plotting\n"
# create poporder file for residual plot
cat $DIR/ALL.clst | cut -f3 | sort | uniq > $DIR/poporder # if benson is present
#cat $DIR/ALL.clst | cut -f3  | sed 's/BENSON//g' | sed '/^$/d' | sort | uniq > poporder.benson