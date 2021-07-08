#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-01
# Last Modified: 2020-05-25


"""  SQ header in bam files notes the NCBI Assecsion number 
	rather than the chromosome number. Update them."""

###########################################
################# Modules #################
###########################################

import re
import sys
from Bio import SeqIO

###########################################
############## Function(s) ################
###########################################

def assembly_table(tab):

	""" Handle assembly table txt file """

	# open reference table file
	lines = open(tab,"r").readlines()

	### list of names and NCBI assecsion numbers ###
	store = []

	# create list of ncbi codes and chrX values 
	for n, elem in enumerate(lines):
		if n > 30 and "#" not in elem:
			# ignore the first 30 entries (table metadata) 
			tmp = elem.split()
			#store.append(tmp)
			if tmp[1] == 'assembled-molecule':
				if tmp[0] == 'MT': # mitochondrial entry has a shorter row
					val = tmp[9]
				else:
					val = tmp[0]
			else:
				val = tmp[0] # contigs

			store.append([tmp[6], val])

	return store


def bam_head(ref_table, bam_file):

	""" pass a bam file and retrun a new header """ 

	assembly = assembly_table(ref_table)

	# open bam header 
	# open file
	bam = open(bam_file,"r").readlines()

	# regex pattern
	pattern = "(?<=@SQ\\tSN:)(.*)(?=\\tLN:)"

	### update the bam list ###
	for n, elem in enumerate(bam):

		### find ncbi num and replace with chrX ###
		try: 
			# search the element for pattern
			ncbi = re.search(pattern, elem).group() 

			# pull the matching name
			new_name = [i[1] for i in assembly if i[0] == ncbi][0]

			# substitute in the new name
			new_line = re.sub(pattern, new_name, elem)

			bam[n] = new_line # update the line

		except AttributeError:
			# skip if pattern doesn't match
			pass 

	# convert to single string 
	bam_string = ''.join(bam)

	return bam_string


def fasta_head(ref_table, fasta_file):

	""" pass a fasta file and retrun a new file """ 

	fasta_file_new = fasta_file + '.corrected.fna'

	assembly = assembly_table(ref_table)
	fasta = open(fasta_file,"r").readlines()

	# regex pattern
	pattern = "(?<=\>)(.*)(?= Equus)"

	### read write fasta file ###
	with open(fasta_file) as original, open(fasta_file_new, 'w') as corrected:
	    records = SeqIO.parse(fasta_file, 'fasta')
	    
	    # scaffold counter - only print once
	    counter = 0

	    print("ID changes:")

	    # iterate through records
	    for record in records:
	    	
	    	# pull the matching name and update ID
	    	record.id = [i[1] for i in assembly if i[0] == record.id][0]
	    	
	    	# keep track of what is happening 
	    	if 'chr' in record.id:
	    		print(record.id)
	    	elif 'scaffold' in record.id and counter == 0:
	    		print('un-assigning scaffolds. (Printing once as many)')
	    		counter += 1
	    	
	    	# write to new fasta
	    	SeqIO.write(record, corrected, 'fasta')


def main(argv):

	""" Main entry point.
		Generate new string header """ 

	if argv[1] == 'bam':
		bam_head(argv[2], argv[3]) # ref table and BAM header
	elif argv[1] == 'fasta':
		fasta_head(argv[2], argv[3]) # ref table and fasta file
		#print(argv)
	else:
		return "Error: file type not selected"
	return 0


if __name__ == "__main__":

    """Makes sure the "main" function is called from the command line"""
    
    status = main(sys.argv) 
    sys.exit(status)

