#!/usr/bin/python3

# Reads in output from sra_lookup.sh and downloads SRA files then extracts adn gzips the fastq file(s). 

import argparse
import sys
import subprocess
import os
import shutil

###############################################################################
# ARGUMENTS
###############################################################################

parser = argparse.ArgumentParser(description = "Script for automation of SRA file downlading and fastq file extraction. NOTE: make sure sratoolkit is loaded! Also, make sure that there is only one SRA file per sample")

parser.add_argument("file", help = "Tab delimitted file containing sample name in first column and SRA name in second column (use - for stdin)")
parser.add_argument("-d", "--directory", help = "Destination directory for fastq files", default = ".")

args = parser.parse_args()

# Deal with stdin vs file input
if args.file == "-":
	input_file = sys.stdin
else:
	input_file = open(args.file, "r")

# Remove trailing / in drectory name
if args.directory.endswith("/"):
	args.directory = args.directory[:-1]

# Throw an error if the target directory doesnt exist
if not os.path.exists(args.directory):
		sys.exit("ERROR: Invalid target directory")

###############################################################################
# VARIABLES
###############################################################################

# The directory in which SRA files downloaded with prefetch are stored
sra_dir = "/nas02/home/s/f/sfrenk/ncbi/public/sra"

###############################################################################

# Process input
for line in input_file:
	sample = line.split("\t")[0]
	print("Getting sequencing data for", sample)
	
	sra = line.strip().split("\t")[1]
	
	# prefetch command downloads .sra file matching the query name and stores it in the local sra drectory 
	subprocess.call(["prefetch", sra])
	
	# move sra file to destination directory
	# need to use shutil rather than move, as it's likely that files will be moving between disks
	shutil.move(sra_dir+"/"+sra+".sra", args.directory+"/"+sra+".sra")
	
	subprocess.call(["fastq-dump", "--outdir", args.directory, "--split-files", args.directory+"/"+sra+".sra"])

	# split-files will either create one file with a _1 suffix for single-end, or two files with _1 and _2 suffixes respictively for paired-end.

	if os.path.exists(args.directory+"/"+sra+"_2.fastq"):

		# paired-end
		shutil.move(args.directory+"/"+sra+"_1.fastq", args.directory+"/"+sample+"_1.fastq")
		shutil.move(args.directory+"/"+sra+"_2.fastq", args.directory+"/"+sample+"_2.fastq")
	else:

		# single end
		shutil.move(args.directory+"/"+sra+"_1.fastq", args.directory+"/"+sample+".fastq")

	os.remove(args.directory+"/"+sra+".sra")
	# gzip all the fastq files to save disk space
	print("gzipping fastq files for ", sample)
	for f in os.listdir(args.directory):
		if f.endswith(".fastq"):
			subprocess.call(["gzip", args.directory+"/"+f])
