# lcGWS_EquCab
> Tool to obtain ancestry and genotype information from a mouth swab sample of a horse sequenced at low-coverage 

## Table of contents
* [General info](#general-info)
* [Requirements](#requirements)
* [Setup](#setup)
* [Retrieving the raw data](#Retrievingtherawdata)
* [Mapping](#Mapping)
* [Merging](#Merging)
* [Deduplicate and clip overlapping read pairs](#Deduplicateandclipoverlappingreadpairs)
* [SNP calling](#SNPcalling)
* [Phasing of the reference panel](#Phasingofthereferencepanel)
* [SNP calling and imputation](#SNPcallingandimputation)
* [Status](#status)
* [Contact](#contact)

## General info
This project aims to develop a pipeline that allows to obtain ancestry and phenotype information from a horse sample sequenced at low-coverage. It was developed as part of my Masters Thesis in Computational Methods in Ecology and Evolution. This repo contains all the scripts used to achieve the Thesis' aim and instructions to obtain the raw data and reproduce the analysis. In order to automate the process I additionally developed two nextflow pipelines: one for mapping of the horse sample and one for imputation of the missing genotype positions to allow ancestry and phenotype analysis available at [nf-pipelines](https://github.com/MaddalenaCella/nf-lcWGS-mapping-and-imputation). 

## Requirements
In order to utilise this tool and reproduce my analysis, make sure have access to an HPC or HTC. 
The computing languages and bioinformatics tools used are the following:

### Languages
* bash (4.2.46)
* Python (3.6.1)
* R (3.6.1)

### Bioinformatic Tools
* BamUtils (1.0.14)
* BWA (0.7.8)
* Bowtie2 (bowtie/2.2.6)
* fastq-tools (0.8.3)
* FastQC (0.11.9)
* FastX (0.0.14)
* Java (jdk-8u144)
* Picard (2.6.0)
* Samtools (1.3.1)
* ShapeIt (2.778)
* Beagle (4.0)
* evalAdmix (0.95)
* Eigensoft (3.0)
* AdmixTools (7.0.2)

### Python modules 
* argparse
* copy
* csv
* numpy
* os
* pandas
* re
* matplotlib

### R packages
* tidyverse
* pophelper
* dplyr
* ggplot2
* gridExtra
* tools
* RColorBrewer
* ggforce
* gghighlight
* ggpubr
* scales
* grid

## Setup
Before reproducing this analysis, make sure you have the required tools in the correct version installed on your HPC machine.

## Retrieving the raw data
* The horse sample analysed is available upon request from Vincent Savolainen HPC environment. To run the analysis using the script provided in this repo, move the reads to `data/Benson`.

* The reference panel raw sequences (FASTQ format) can be downloaded from ENA (European Nucleotide Archive) by submitting the job `job_submissions/fastqDownload.sh`

* The Horse reference sequence is available to download from NCBI by submitting the job `job_submissions/reference_genome.sh`. By running this script you will also update the header to contain information on the chromosome position and index the reference genome (bwa index).

## Mapping
* To map the raw sample reads to the Horse reference genomes submit the job `job_submissions/novelMapper.sh`.
* To map the reference panel sequences submit the job `job_submissions/reference_mapping.sh`.

## Merging
* To remove reads with quality lower than 10 and merge the mapped reads of the horse sample, submit the job `job_submissions/merge.sh`.
* to merge the mapped reads of the reference panel you need to:
1. generate a bam list by running `sh general/makebam.sh -m` on the remote cluster;
2. `scp` the bam list (data/bam.list) to your local machine and run `python3 mapping/toMerge.py -b data/bam.list -o data/ -i data/cleaned_data/info_all.csv -r Run -g BioSample` to generate two csv files with files that require and do not require merging;
3. run `Rscript --vanilla mapping/move_nomerge.r` to generate a move.list of files that do not require merging;
4. `scp` the to_merge.csv, no_merge.csv and move.list files to the HPC;
5. now you can submit the job `job_submissions/merge_ref.sh` as a job array with the same number of jobs as the to_merge.csv file (`wc -l data/to_merge.csv`);
6. move the files that do not require merging to the same folder as the merged files by running `sh mapping/move_nomerge.sh`;
7. regenerate the bam.list with the merged reference panel `sh general/makebam.sh -m` (comment the line used to generate the first list and un-comment the following line).

## Deduplicate and clip overlapping read pairs
* To deduplicate and clip overlapping read pairs on the merged horse bam file submit the job `job_submissions/dedup_clipover.sh`;
* for the reference panel run `job_submissions/dedup_clipover_wgs.sh`;
* to calculate the average depth after this step in the reference panel run: 
```
for i in results/wgs_data/depth/*.txt; do cat $i>> results/wgs_data/depth/depths.txt; done

awk '{ total += $1 } END { print total/NR }' results/wgs_data/depth/depths.txt > results/wgs_data/depth/avg_depth.txt
```

## SNP calling
ANGSD was used to call genotypes (p-value to call SNPs set to -SNP_pval 1e-6).
* To perform variant calling on the reference panel run `job_submissions/SNP_calling_ref.sh`.

## Phasing of the reference panel
ShapeIt was used to phase the reference panel. The phasing algorithm was performed on each individual chromosome.
* To phase the reference panel per chromosome run `job_submissions/phasing_ref_shapeit.sh`;
* to concatenate the phased vcf files run `job_submissions/concatRef.sh`.

## SNP calling and Imputation
Angsd was used for a preliminary variants calling (p-value to call SNPs set to -SNP_pval 1e-1). To speed up the SNP calling algorithm, the search of SNPs in the horse sample file was limited to the variants identified in the reference panel. Variants calling was followed by imputation of the missing geneotypes with Beagle 4.0. 
* To call variants and impute missing genotypes run the job `job_submissions/SNP_caling_imputation.sh`

## Ancestry analysis
* First, remove singletons and sites in LD with vcftools and Plink respectively by submitting `job_submissions/LD_singletons_pruning.sh`;
* to create a --keep file to run the PCA also without Przewalski and Przewalski hybrids run `analysis/PCA_no_Prz.r`;
* to conduct the PCA analysis with and without Przewalski and Przewalski hybrids submit `job_submissions/PCA_plink.sh`;
* to conduct the ancestry analysis submit `job_submissions/ancestry.sh`;
* to evaluate the results of admixture submit `job_submissions/evalAdmix.sh` as a job array with one job per number of ancestral K populations tried;
* scp the PCA results, the .Q files from the admixture analysis and the .txt files from the evalAdmix runs locally in the results/ancestry directory;
* to extract the basename from the bam file list `for F in $(cat data/bam1.list); do     basename $F _dedup_overlapclipped.bam >> data/ancestry/identifiers.txt; done`;
* run the scripts `analysis/Breeds_table.r`, `analysis/PCA_plotting.r` and `analysis/admixture_plotting.r` to produce the plots of interest;
* create a parameter file to convert ped to eigenstrat format and submit the job `job_submissions/ped_to_eigenstrat.sh`. The parameter file should look like this:

```
genotypename:    treemix.ped
snpname:         treemix.map # or example.map, either works 
indivname:       treemix.ped # or example.ped, either works
outputformat:    EIGENSTRAT
genotypeoutname: treemix.eigenstratgeno
snpoutname:      treemix.snp
indivoutname:    treemix.ind
familynames:     NO
```
* scp locally the .ind file produced from the above analysis and run `analysis/ind_add_pop.r` to add the population individuals belong to the file and move the file back to the HPC; 
* generate all possible combinations of source populations running the R script `analysis/generate_pops.r`. Move the .tsv file obtained to HPC to run the f3 analysis or alternatively, create a customised tsv file with the combinations of populations of interest;
* create the parameter file for the f3 analysis including the eigenstrat files (.snp, .ind, .geno) and the populations file and submit `job_submissions/treemix.sh`. The parameter file should look like this:

```
genotypename: treemix.eigenstratgeno
snpname: treemix.snp
indivname: treemixpop.ind ##modified .ind file obtained from ind_add_pop.r
popfilename: pop_combos_all.tsv ##tab-delimited file with population combinations
```
* to filter out the unnecessary information from the output of the f3 analysis and keep just the rows starting with 'result:' use `grep 'result:' f3stat_poplist.txt > f3stat_poplist_filtered.txt`;
* to visualise and summarise the results of the f3 analysis in plots and tables run `analysis/f3_analysis.r` and `analysis/f3_plots.py`;


## Plotting
The plots used in the write-up were created with the following scripts:
* The table containing the number of individuals in the reference panel and their breeds was produced with the R script `analysis/Breeds_table.r`;
* PCA results were plotted with the R script `analysis/PCA_plotting.r`;
* cross-validation, Admixture and evalAdmix plots were produced by running `analysis/admixture_plotting.r`;


## Phenotype inference
To find whether some SNPs associated with phenotypes of interest were present on the horse sample you need to:
1. generate a list of genomic positions for the QTLs and Mendelian traits of interest run `Rscript data/gene_variants/all_traits.r`. The csv file (all_traits.csv) from which the genomic positions were obtained was manually curated by me based on the article: Ten years of the horse reference genome: insights into equine biology, domestication and population dynamics in the post-genome era.
2. covert between a list of genomic positions in the format chromosome:start position-end position run `cat data/gene_variants/angsd_pos.txt | sed 's/:/\t/' | sed 's/-/\t/' > data/gene_variants/angsd_final.file`.
3. finally, run the script `job_submissions/pheno_search.sh`.

To produce tables with all Mendelian traits and QTLs included in the anlysis, as well as the phenotypic variants found in the horse sample `analysis/phenotype.r`.
 
## Status
Project is: _in progress_

## Contact
Created by [@MaddalenaCella](mc2820@ic.ac.uk) - feel free to contact me!
