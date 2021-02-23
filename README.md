# lcGWS_EquCab
> Tool to obtain ancestry and genotype information from a mouth swab sample of a horse sequenced at low-coverage 

## Table of contents
* [General info](#general-info)
* [Requirements](#requirements)
* [Setup](#setup)
* [Retrieving the raw data](#Retrieving the raw data)
* [Mapping](#Mapping)
* [Merging](#Merging)
* [Deduplicate and clip overlapping read pairs](#Deduplicate and clip overlapping read pairs)
* [Loimpute input](#loimpute input)
* [Loimpute phased reference](#loimpute phased reference)
* [Status](#status)
* [Contact](#contact)

## General info
Add more general information about project. What the purpose of the project is? Motivation?

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
* fastq-tools (0.8.3)
* FastQC (0.11.9)
* FastX (0.0.14)
* gatk (4.0)
* Java (jdk-8u144)
* Picard (2.6.0)
* Samtools (1.3.1)

### Python modules 
* argparse
* copy
* csv
* numpy
* os
* pandas
* re

### R packages
* tidyverse

## Set-up


## Retrieving the raw data
* The horse sample analysed is available upon request from Vincent Savolainen HPC environmnet. 

* The reference panel raw sequences (FASTQ format) can be downloaded from ENA (European Nucleotide Archive) by submitting the job `job_submissions/fastqDownload.sh`

* The Horse reference sequence is available to download from NCBI by submitting the job `job_submissions/reference_genome.sh`. By running this script you will also update the header to contain information on the chromosome position and index the reference genome (bwa index).

## Mapping
* To map the raw sample reads to the Horse reference genomes submit the job `job_submissions/novelMapper.sh`.
* To map the reference panel sequences submit the job `job_submissions/wgsMapper.sh`.

## Merging
* To merge the mapped reads of the horse sample submit the job `job_submissions/merge.sh`.
* to merge the mapped reads of the reference panel you need to:
1. generate a bam list by running `sh general/makebam.sh -m` on the remote cluster.
2. `scp` the bam list (data/bam.list) to your local machine and run `python 3 mapping/toMerge.py -b data/bam.list -o data/ -i data/cleaned_data/info_all.csv -r Run -g BioSample` to generate two csv files with files that require and do not require merging.
3. run `Rscript --vanilla mapping/move_nomerge.sh` to generate a move.list of files that do not require merging.
4. `scp` the to_merge.csv, no_merge.csv and move.list files to the HPC.
5. now you can submit the job `job_submissions/wgs_merge_2.sh` as a job array with the same number of jobs as the to_merge.csv file (`wc -l data/to_merge.csv`).
6. move the files that do not require merging to the same folder as the merged files by running `sh mapping/move_nomerge.sh`.
7. regenerate the bam.list with the merged reference panel `sh general/makebam.sh -m` (comment the line used to generate the first list and un-comment the following line).

## Deduplicate and clip overlapping read pairs
* to deduplicate and clip overlapping read pairs on the merged horse bam file submit the job `job_submissions/dedup_clipover.sh`.
* for the reference panel run `job_submissions/dedup_clipover_wgs.sh`.

## LowImpute Input 
To prepare the horse sample for imputation using Loimpute submit the job `job_submissions/loimpute_input.sh`.

## LowImpute phased reference
To make the reference panel ready to be used for imputation you need to run to:
* realign indels by submitting the job `job_submissions/indel_realign.sh`.
* calling haplotypes by submitting the job `job_submissions/haplotype_caller.sh`.


 
## Status
Project is: _in progress_

## Inspiration
Add here credits. Project inspired by..., based on...

## Contact
Created by [@MaddalenaCella](mc2820@ic.ac.uk) - feel free to contact me!
