#!/usr/bin/env python

import argparse
import os
import sys

parser = argparse.ArgumentParser()

parser.add_argument("file", help = "input file ('-' for stdin)")
parser.add_argument("-o", "--output", help = "output filename", default = "uncollapsed.txt")

args = parser.parse_args()

if args.file == "-":
	f = sys.stdin
else:
	f = open(args.file, "r")

if os.path.exists(args.output):
	os.remove(args.output)
outfile = open(args.output, "a")

for line in f:
	if not line.startswith(">"):
		for i in range(int(line.strip().split("\t")[1])):
			outfile.write(line.strip().split("\t")[0] + "\n")

f.close()
outfile.close()
