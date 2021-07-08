# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-05
# Last Modified: 2020-06-25



""" merge wgs_all bam files """

###########################################
################# Modules #################
###########################################

import csv
import os

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


# timer
import time
start = time.time() # start the timer from import

def timer():
	
	end = time.time()
	duration = end-start
	duration = round(duration)
	string = "\n..........................\n"\
			"   Time elapsed: {} sec"\
			"\n..........................\n"\
			.format(duration)

	print(string)



###########################################
######### Input(s) and Parameters #########
###########################################

#make it so that it takes arguments from command
to_merge = open_csv("/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/data/to_merge.csv")
directory = "/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/job_submissions/wgs_data/sorted/"
out_directory = "/rds/general/project/human-popgen-datasets/live/Maddalena/EquSeq/job_submissions/wgs_data/merged/"

###########################################
############### Wraggling #################
###########################################
del to_merge[0]

job_num = int(os.environ["PBS_ARRAY_INDEX"])

job = to_merge[job_num]

out_name = out_directory + job[0] + ".bam"
files = [directory + i + ".sorted.bam" for i in job[1:]]

files_string = " ".join(files)


# merge
print("merging")
cmd_merge = "samtools merge --threads 31 {} {}".format(out_name, files_string)

#print(cmd_merge)
os.system(cmd_merge)

timer()

# index
print("indexing")
cmd_index = "samtools index {} {}.bai ".format(out_name, out_name)
#print(cmd_index)
os.system(cmd_index)

timer()