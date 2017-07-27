#!/usr/bin/bash

module add fastqc
module list

# Defaults
output="."

usage="
	USAGE:
		This script maps fastq files and converts sam files to bam. Make sure the following modules are loaded:
		bwa
		samtools

		Trims adapter sequence from fastq files.

		-d/--directory: directory containing fastq files
		-o/--output: output directory (default: current directory)
"

if [ -z "$1" ]; then
    echo "$usage"
    exit
fi

while [[ $# > 0 ]]
do
    key="$1"
    case $key in
        -d|--directory)
		directory="$2"
		shift
		;;
		-o|--output)
		output="$2"
		shift
		;;
    esac
shift
done

# Check directories
if [ ! -d $directory ]; then
	echo "ERROR: Invalid fastq directory path"
	exit 1
fi

if [ ! -d $output ]; then
	echo "ERROR: Invalid output directory path"
	exit 1
fi

# Get adapter sequences for library type
#if [ srna = false ]; then
#	adapter=/nas02/apps/bbmap-36.64/bbmap/resources/truseq.fa.gz
#else
#	adapter=/nas02/apps/bbmap-36.64/bbmap/resources/truseq_rna.fa.gz
#fi

for file in ${directory}/*.fastq*; do
	fastqc -o $output $file
done
