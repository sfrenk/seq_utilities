#!/usr/bin/env python

###############################################################################
# Merge multiple bed files containing count data
###############################################################################

import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description="Merge multiple count files into a single count table")
parser.add_argument("directory", help="directory containing count files (each file must end with '.bed'")
parser.add_argument("-o", "--output", help = "output file name (default: count_table.txt", default = "bed_count_table.txt")
parser.add_argument("-t", "--total", help = "total mapped reads file for Normalization (optional)")

args = parser.parse_args()

count_files = []

for f in os.listdir(args.directory):
	if f.endswith(".bed"):
		count_files.append(f)

# Read in first count file. All further counts will be appended to this data frame
bed_name = count_files[0].replace(".bed", "").replace("_sorted", "")
merged = pd.read_csv(count_files[0], sep="\t", header=None, names=['chr', 'start', 'end', bed_name])

# Remove the first count file from the list so that it doesn't get counted twice

count_files = count_files[1:]

# Merge the other count files to the first one

for i in count_files:
	bed_name = i.replace(".bed", "").replace("_sorted", "")
	count_file = pd.read_csv(i, sep="\t", header=None, names=['chr', 'start', 'end', bed_name])

	merged = merged.merge(count_file, how='left', on=['chr', 'start', 'end'])

merged = merged.fillna(0)

# Merge the chromosome, start and end variables into one column

def make_locus(chr, start, end):
	return(str(chr)+":"+str(start)+".."+str(end))

merged['locus'] = merged.apply(lambda x: make_locus(x['chr'], x['start'], x['end']), axis=1)

# Rearrange columns

indexes = [x for x in range(3,len(merged.columns)-1)]
indexes = [len(merged.columns)-1] + indexes


merged = merged[indexes]

# Save the count table

merged.to_csv(args.output, sep='\t', index=False)

# Normalization

if args.total is not None:
	merged_norm = merged
	total_mapped = pd.read_csv(args.total, sep = "\t", header = None, names = ["sample", "count"])

	for i in merged_norm:
		if i != "locus":
			mapped_count = total_mapped.loc[total_mapped["sample"] == i, "count"]
			merged_norm[i] = merged_norm[i].apply(lambda x: x/mapped_count)

	merged_norm.to_csv(args.output+".norm", sep='\t', index=False)
	