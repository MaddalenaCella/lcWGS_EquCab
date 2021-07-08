#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-26
# Last Modified: 2020-06-26



""" Determine which bam files belong to which population.
	Lists will be used to generate summary stats """

###########################################
################# Modules #################
###########################################

import csv
import os
import pandas as pd
import re
import numpy as np

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

### save a text file without a new line at the end
def saveTxt(dirfile, listToSave, sep='\n'):
	with open(dirfile, 'w') as f:
		for num, val in enumerate(listToSave):
			if num == len(listToSave) - 1:
				f.write(val)
			else:
				f.write(val + sep)


# strip a filepath to its basename
stripToBase = lambda elem : os.path.basename(elem).split("_")[0]

###########################################
######### Input(s) and Parameters #########
###########################################

fopen = open("data/bam1.list", "r")
list_raw = fopen.readlines()
list_full = [i.replace("\n", "") for i in list_raw]
info = pd.read_csv("data/cleaned_data/info_all.csv")

###########################################
############### Wraggling #################
###########################################

print("joining groups")
# full paths and respective basenames
list_df = pd.DataFrame(list_full)
list_df.columns = ["fullpath"]
list_df["basename"] = [stripToBase(i) for i in list_df["fullpath"]]

# info table trimmed to required columns
info_trim = info[['Run', 'BioSample', 'sub_group']]

# join on run
run_paths = pd.merge(list_df, info_trim, how='inner', left_on="basename", right_on="Run")

# join NAs on BioSample 
	# make a unique df 
bio_paths = pd.merge(list_df, info_trim[["BioSample", "sub_group"]], \
	how='inner', left_on="basename", right_on="BioSample").drop_duplicates()

# add benson
#benson_path = list_df[list_df["basename"] == "final"]
#benson_path["sub_group"] = "BENSON"
# concat the dfs
con_elems = {"runs": run_paths, "bioS": bio_paths }
df_all = pd.concat(con_elems, sort=True)

# convert to list 
df_all_trim = df_all[["fullpath", "sub_group"]]
list_data = df_all_trim.values.tolist()

# remove spaces and special characters in grp name
for val, elem in enumerate(list_data):
	tmp = elem[1].replace(" ", "_") # remove spaces
	list_data[val][1] = re.sub('[^A-Za-z0-9]+', '', tmp) # remove all special characters

# unique groups - prep for saving lists
uniq_grps = np.unique([i[1] for i in list_data])

print("saving lists")
# save lists based on matches of groups
for grp in uniq_grps:
	path_list = [i[0] for i in list_data if i[1] == grp]
	path_name = "data/ancestry/bam_list_grps/" + grp.lower() + ".list"
	saveTxt(path_name, path_list)