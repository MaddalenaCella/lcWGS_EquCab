#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-16
# Last Modified: 2020-06-24



""" determine which fastq files need to be merged """



# python3 scripts/toMerge.py -b data/bam.list -o data/ -i data/cleaned_data/info_all.csv -r Run -g BioSample

###########################################
################# Modules #################
###########################################

import csv
import os
import argparse
import numpy as np

###########################################
################# Options #################
###########################################

parser = argparse.ArgumentParser(description='Which files need merging.')

# input list
parser.add_argument("-b", "--bam.list", dest="bamlist", type=str,
					required=True, metavar="IN_LIST",
					help="IN_LIST.csv of unique RUN ID and grouping, eg. individual")

# input file
parser.add_argument("-i", "--in.file", dest="infile", type=str,
                 	 help="IN_FILE with information on RUN_ID and GROUP_ID",
                 	 required=True, metavar="IN_FILE")

# output list prefix
parser.add_argument("-o", "--out.list", dest="outfile", type=str,
                 	 help="write OUT_LIST of RUN IDs to merge and the grouping",
                 	 required=True, metavar="OUT_LIST")
# unique id
parser.add_argument("-r", "--runID", dest="runid", type=str,
                  required=True, help="RUN_ID name, should be unique.",
                  metavar="RUN_ID")
# group by key
parser.add_argument("-g", "--groupby", dest="grp", type=str,
                  required=True, help="group by GROUP_ID variable name.",metavar="GROUP_ID")


# define args
args = parser.parse_args()

###########################################
############## Function(s) ################
###########################################


def open_csv(file):

	""" open a csv into a list format """

	tmp = [] # initialise the list
	with open(file, 'r') as f:
		reader = csv.reader(f)
		for row in reader:
			tmp.append(row) # add row to list

	return tmp

def write_csv(list_file, path):
	
	""" Write list to csv """

	with open(path, 'w') as f:
		writer = csv.writer(f, delimiter=',')
		for i in list_file:
			if isinstance(i, list):
				writer.writerow(i) # multi-column
			else:
				writer.writerow([i]) # single column

def stripPE(string):
	""" strip path and extension from 
		file string """
	noext = string.split(os.extsep, 1)[0]
	elem = os.path.split(noext)[1]
	return elem

###########################################
######### Input(s) and Parameters #########
###########################################

#bamlist = args.bamlist
#args.runid  # Run
#args.grp # BioSample
#args.outfile = 'data/cleaned_data/to_merge.csv'
#args.infile #'data/cleaned_data/info_all.csv'

#info_all = open_csv('data/cleaned_data/info_all.csv')
#bam_txt = open('data/bam.list', "r")

info_all = open_csv(args.infile)
bam_txt = open(args.bamlist, "r")
bams = [i.split()[0] for i in bam_txt.readlines()]

###########################################
############### Wraggling #################
###########################################

bamlist_no_ext = [stripPE(i) for i in bams]

i_run = info_all[0].index(args.runid)
i_grp = info_all[0].index(args.grp)

# only files present in bamlist
present_only = [i for i in info_all if i[i_run] in bamlist_no_ext]

# unique codes - i.e. don't group or merge
groups = [i[i_grp] for i in present_only]

# return non unique grouping codes
not_unique_vals = [i for i in {*groups} if groups.count(i) > 1]

store = []
for val in not_unique_vals:
	tmp = [i[i_run] for i in present_only if val == i[i_grp]]
	tmp.insert(0, val)
	store.append(tmp)

# and unique ones
unique_vals = [i for i in {*groups} if groups.count(i) == 1]

store_uniq = []
for val in unique_vals:
	tmp = [i[i_run] for i in present_only if val == i[i_grp]]
	tmp.insert(0, val)
	store_uniq.append(tmp)


# append a header
store.insert(0, [args.grp, args.runid + '*'])

# write to file
path_merge = args.outfile + "to_merge.csv"
path_not = args.outfile + "no_merge.csv"

write_csv(store, path_merge) # those to merge 
write_csv(store_uniq, path_not) # those not to merge




