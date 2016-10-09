#!/usr/bin/env python

###############################################################################
# Extract all 22G sRNAs fastq, fasta or raw file. This gives a .txt file containing raw sequences, which can be mapped in bowtie using the -r flag
###############################################################################

import os
import argparse
import gzip

parser = argparse.ArgumentParser(description="Extract 22g and/or 21u rnas from fastq, fasta or raw text file")

parser.add_argument("I", help = "Input file. File type will be detected based on the first line of the file. File may be in gzipped format.")
parser.add_argument("-o", help = "Output file.", required = True)
parser.add_argument("-g", help = "Retrieve 22g RNAs.", action = "store_true")
parser.add_argument("-u", help = "Retrieve 21u RNAs.", action = "store_true")
parser.add_argument("-s", help = "If not filtering for 22g or 21u RNAs, you can specify size range of reads to keep by providing min and max length seperated by a comma (eg. to keep reads between 19 and 24 nucleotides (inclusive), use '-s 19,24' (Default size range: 18 to 30 nucleotides)", default = "18,30")

args = parser.parse_args()

infile = args.I
outfile = args.o

# Define minimum and maximum read size

min_size = int(args.s.split(",")[0])
max_size = int(args.s.split(",")[1])

###############################################################################

print '##### small_rna_filter.py #####'

# Input file may or may not be gzipped

if infile.endswith(".gz"):
	f = gzip.open(infile, "rb")
else:
	f = open(infile, "r")

seqs = f.readlines()

reads = []

if seqs[0].startswith("@"):
	print 'extracting sequence information from fastq file'
	# Each fastq record consists of four lines. The sequence info is on the second line.
	for i, line in enumerate(seqs):
		if i % 4 == 1:
			reads.append(line)
	reads = [line.strip() for line in reads]

elif seqs[0].startswith(">"):
	print 'extracting sequence information from fasta file'
	reads = [line.strip() for line in seqs if line[0] != '>']

elif seqs[0][0].isalpha():
	print 'extracting sequence information from raw text file'
	reads = [line.strip() for line in seqs]

else:
	raise IOError("Invalid input file type")

f.close()

print 'input reads:'
print(len(reads))

G22 = []
U21 = []

# Extract 22G RNAs

if args.g:
	G22 = [x for x in reads if len(x) == 22 and x[0] == 'G']
	print '22G RNAs extracted:'
	print(len(G22))

# Extract 21U RNAs

if args.u:
	U21 = [x for x in reads if len(x) == 21 and x[0] == 'T']
	print '21U RNAs extracted:'
	print(len(U21))

# Combine all filtered reads

reads_out = G22 + U21

# If neither -g nor -u were specified, just convert input reads to raw format

if args.g == False and args.u == False:
	print('no sRNA type selected. Selecting for reads that are ' + str(min_size) + ' to ' + str(max_size) + ' nucleotides in length and converting all reads to raw format')
	reads_out = reads

	# Select reads that fall in the expected size range (18-30 nucleotides)

	reads_out = [x for x in reads_out if len(x) >= min_size and len(x) <= max_size]

# Write out all (filtered) reads to raw txt file

print 'writing filtered sequences to:'
print(outfile)

with open(outfile, 'w') as f:
	for line in reads_out:
		f.write(line+'\n')

print('#####	Finished	#####')
