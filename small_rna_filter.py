#!/usr/bin/env python

###############################################################################
# small_rna_filter.py version 2.0
#
#
# Filter reads from a fastq, fasta or raw file. Output is a .txt file containing raw sequences, which can be mapped in bowtie using the -r flag
###############################################################################

import os
import argparse
import gzip

parser = argparse.ArgumentParser(description="Extract 22g and/or 21u rnas from fastq, fasta or raw text file")

parser.add_argument("infile", help = "Input file. File type will be detected based on the first line of the file. File may be in gzipped format.")
parser.add_argument("-o", "--output", help = "Output file.", required = True)
parser.add_argument("-f", "--filter", help = "Filter reads based on the first nucleotide (eg. to select 22g RNAs, use options -f g -s 22,22", default = "A,T,G,C")
parser.add_argument("-s", "--size", help = "Specify size range of reads to keep by providing min and max length seperated by a comma (eg. to keep reads between 19 and 24 nucleotides (inclusive), use '-s 19,24' (Default size range: 18 to 30 nucleotides)", default = "18,30")
parser.add_argument("-a", "--trim_a", help = "trim 3' A nucleotides", action = "store_true", default = False)

args = parser.parse_args()

# Define starting nucleotide
args.filter = args.filter.upper()
first_base = args.filter.split(",")

# Define minimum and maximum read size

min_size = int(args.size.split(",")[0])
max_size = int(args.size.split(",")[1])

###############################################################################

print("##### small_rna_filter.py #####")

# Input file may or may not be gzipped

if args.infile.endswith(".gz"):
	f = gzip.open(args.infile, "rb")
else:
	f = open(args.infile, "r")

seqs = f.readlines()

reads = []

if seqs[0].startswith("@"):
	print("extracting sequence information from fastq file")
	# Each fastq record consists of four lines. The sequence info is on the second line.
	for i, line in enumerate(seqs):
		if i % 4 == 1:
			reads.append(line)
	reads = [line.strip() for line in reads]

elif seqs[0].startswith(">"):
	print("extracting sequence information from fasta file")
	reads = [line.strip() for line in seqs if line[0] != '>']

elif seqs[0][0].isalpha():
	print("extracting sequence information from raw text file")
	reads = [line.strip() for line in seqs]

else:
	raise IOError("Invalid input file type")

f.close()

print("input reads: " + str(len(reads)))

# Trim 3' A nucleotides

if args.trim_a:
	print("removing 3' A nucleotides")
	trimmed_reads = []
	count_a = 0
	for read in reads:
		while read.endswith('A'):
			# Sequentially remove all 3' A nucleotides
			read = read[:-1]
			count_a += 1
		trimmed_reads.append(read)

	reads = trimmed_reads
	print("3 prime A nucleotides removed: " + str(count_a))

# Extract reads based on filtering parameters

print("extracting reads in size range "+ str(min_size) + " to " + str(max_size) + " nucleotides with " + args.filter + " as the first nucleotide")

reads = [x for x in reads if len(x) >= min_size and len(x) <= max_size and x[0].upper() in first_base]

# Write out all (filtered) reads to raw txt file

print("output reads: " + str(len(reads)))

print('writing filtered sequences to: ' + args.output)

with open(args.output, 'w') as f:
	for line in reads:
		f.write(line+'\n')

print("#####	Finished	#####")
