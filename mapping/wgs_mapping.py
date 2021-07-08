#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-09
# Last Modified: 2020-05-23



""" Which sra files need aligning together? """

###########################################
################# Modules #################
###########################################

import csv
import os
import sys

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


def pairFiles(fileList, runList):

	""" pair files in list with run list """

	store_files = []
	for i in fileList:
		name = ''
		position = ''
		
		i.split('.fastq.gz')
		filename, file_extension = i.split('.', 1)
		file_extension = '.' + file_extension

		# paired reads
		if '_' in filename:
			name, position = filename.split('_')
		else:
			name = filename


		store_files.append([name, position, i])



	store = []
	for ru in runList:
		tmp = [i[2] for i in store_files if i[0] == ru and i[1] != '']

		# if longer than one, paired and unmapped are included
		if len(tmp) > 1:
			tmp_sorted = sorted(tmp)
			tmp_sorted.append(ru)
			store.append(tmp_sorted)

	return store
		


###########################################
################# main ####################
###########################################


def main(argv):

	""" Run pairs functions return 
		job number specific pair """ 
	# catch the job number
	runs_file = argv[1] 
	
	with open(runs_file) as f:
		runs = f.read().splitlines()

	files = os.listdir(argv[2])
	
	pairs = pairFiles(files, runs)
	job_num = int(os.environ['PBS_ARRAY_INDEX'])
	job_pair = pairs[job_num]
	
	print(job_pair) # return the pair list

	return 0


if __name__ == "__main__":

    """Makes sure the "main" function is called from the command line"""
    
    status = main(sys.argv) 
    sys.exit(status)
