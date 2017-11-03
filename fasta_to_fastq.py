#!/usr/bin/env python3

import argparse
import sys

parser = argparse.ArgumentParser(description = "convert fasta or raw sequence to fastq with fake quality scores")

parser.add_argument("file", help = "input filename (use - for stdin)")
parser.add_argument("-o", "--output", help = "output filename")
parser.add_argument("-n", "--name", help = "Name prefix for read names. This will be used when there is no fasta header (default: seq)", default = "name")

args = parser.parse_args()

if args.file == "-":
	infile = sys.stdin
else:
	infile = open(args.file, "r")

outfile = open(args.output, "w")

counter = 1

for line in infile:
	name = None
	
	if line.startswith(">"):
		# fasta header
		name = line(line.strip()[1:])
		f.next()
	seq = line.strip()
	
	if name == None:
		# Create name for raw sequence
		name = args.name + "_" + str(counter)
		counter += 1

	# Construct fake quality line
	qual = "I" * len(seq)

	outfile.write("@" + name + "\n" + seq + "\n+\n" + qual + "\n")

outfile.close()
