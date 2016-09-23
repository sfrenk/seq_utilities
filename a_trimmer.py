#!/usr/bin/env python

###############################################################################
# Trim 3' A nucleotides from raw format reads
###############################################################################

import argparse
import os

parser = argparse.ArgumentParser("Description = trim 3' A nucleotides from raw format reads")

parser.add_argument("file", help = "Input file (must be raw format - one read per line")
parser.add_argument("-o", "--output", help = "Output file name", default = "reads_a_trimmed.txt")

args = parser.parse_args()

# Remove output file if it already exists
if os.exists(args.output):
	os.remove(args.output)


in_file = open(args.file, "r")
out_file = open(args.file, "a")

print("removing 3' A nucleotides")
count_a = 0
for line in in_file:
	read = line.split()
	while read.endswith('A'):
		# Sequentially remove all 3' A nucleotides
		read = read[:-1]
		count_a += 1
	out_file.write(read+"\n")

print('3 prime A nucleotides removed: ', str(count_a))
