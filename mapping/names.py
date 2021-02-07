#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-03-23
# Last Modified: 2020-06-23



""" Generate names to be used for novel data """

###########################################
################# Modules #################
###########################################

import os 
import re
import sys

###########################################
############## Function(s) ################
###########################################

def pairs(job_number, path='/rds/general/user/jb1919/home/genomics/sequences/cdts-hk.genomics.cn/Clean/F19FTSEUHT1854-swab-horse-1A/'):

	""" Determine the pair-ended read names """

	# list of files
	files = os.listdir(path)  
	files.sort()
	
	# remove unwanted .file
	if '.listing' in files:
		files.remove('.listing')

	### new name to create pairs form ###
	# regex patterns
	pat_1 = r"L\d" # first distinguishing aspect
	pat_2 = r"\-\d\d\d" # and end

	store = []

	for i in files:
		
		# extract all info - read and pair
		match_read1 = re.search(pat_1, i)
		match_read2 = re.search(pat_2, i)
		read1 = match_read1.group()
		read2 = match_read2.group()
		read12 = read1 + read2
		# update store
		store.append([i, read12])

	### match the pairs ###
	pairs = []

	for c in store:

		# match the pairs
		p = [i[0] for i in store if i[1] == c[1]]
		
		# append common name
		p.append(c[1]) # to list 
		
		# only append when not in th elist already
			# matches will return doubles
		if p not in pairs:
			pairs.append(p) # to storage

	return pairs[job_number]




###########################################
################# main ####################
###########################################


def main(argv):

	""" Run pairs functions return 
		job number specific pair """ 
	# catch the job number
	job_num = int(os.environ['PBS_ARRAY_INDEX'])
	job_pair = pairs(job_num)
	
	print(job_pair) # return the pair list

	return 0


if __name__ == "__main__":

    """Makes sure the "main" function is called from the command line"""
    
    status = main(sys.argv) 
    sys.exit(status)











