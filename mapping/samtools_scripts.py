#!/usr/bin/env python
# Author: Maddalena Cella
# Email mc2820@ic.ac.uk
# Date: 15-02-2021
# Last Modified: 02-07-2021

""" Produce merge and indexing scripts for each line of the to_merge.csv file """

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

###########################################
######### Input(s) and Parameters #########
###########################################

#it takes arguments from command line

def main(argv):

    to_merge = open_csv(argv[1])
    directory = argv[2]
    out_directory = argv[3]

    ###########################################
    ############### Wraggling #################
    ###########################################
    del to_merge[0]

    job_num = int(os.environ["PBS_ARRAY_INDEX"])
    job = to_merge[job_num]

    out_name = out_directory + job[0] + ".bam"

    files = [directory + i + ".sorted.bam" for i in job[1:]]

    files_string = " ".join(files)



    cmd_merge = "samtools merge {} {}".format(out_name, files_string)

    print(cmd_merge)
    #os.system(cmd_merge)

    cmd_index = "samtools index {} {}.bai ".format(out_name, out_name)
    print(cmd_index)


if __name__ == "__main__":

    """Makes sure the "main" function is called from the command line"""
    
    status = main(sys.argv) 
    sys.exit(status)


