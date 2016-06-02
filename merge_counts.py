#!/usr/bin/env python

###############################################################################
# Use this script to merge all count.txt files in a count directory into a count table. The resulting table can be loaded straight into DESeq2 for analysis
###############################################################################

import pandas as pd
import os
import argparse

# Get target directory

parser = argparse.ArgumentParser(description="Merge multiple count files into a single count table")

parser.add_argument("directory", help="directory containing count files (each file must end with '.txt'")

args = parser.parse_args()

directory = args.directory 
os.chdir(directory)

# Merge files

count_files = []

for f in os.listdir("."):
	if f.endswith("counts.txt"):
		count_files.append(f)

merged = pd.read_csv(count_files[0], sep="\s", header=None, names=[count_files[0], 'gene'])
merged = merged[['gene', count_files[0]]]

for i in range(1, len(count_files)):
	count_file = pd.read_csv(count_files[i], sep="\s", header=None, names=[count_files[i], 'gene'])
	merged = merged.merge(count_file, how='left', left_on='gene', right_on='gene')

# Get rid of the row containing unmapped reads

merged = merged[merged['gene'] != '*']

# If a particular transcript was detected in some samples but not others, the samples in which the transcript was absent will report an "NA" value after merging of the count files. This needs to be changed to "0" for downstream analysis

merged = merged.fillna(0)

# Save the count table

merged.to_csv('count_table.txt', sep='\t', index=False)

