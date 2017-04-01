#!/usr/bin/env bash

###############################################################################
# Collapse reads in raw/fasta format to unique sequences. Output consists of raw reads.
###############################################################################

usage="
    Collapse reads in raw/fasta format to unique sequences. Output consists of raw reads.

    ARGUMENTS
        -i/--input		input directory (default: current directory)
        -o/--output		output directory for collapsed reads (default: current directory)
        -c/--count	output directory for collapsed reads with counts (default: same as output directory)       
    "

# Defaults:
input="."
output="."
count="same as output"

# Parse arguments

while [[ $# > 0 ]]; do
	key="$1"
	case $key in
		-i|--input)
		input="$2"
		shift
		;;
		-o|--output)
		output="$2"
		shift
		;;
		-c|--count)
		count="$2"
		shift
		;;
		-h|--help)
		echo "$usage"
		exit
		;;
	esac
	shift
done

# If no count output directory is selected, use same directory as main output

if [[ $count == "same as output" ]]; then
	count=${output}
fi
# Remove trailing "/" from directories

if [[ ${input:(-1)} == "/" ]]; then
	input=${input::${#input}-1}
fi

if [[ ${output:(-1)} == "/" ]]; then
	output=${output::{#output}-1}
fi

if [[ ${count:(-1)} == "/" ]]; then
	count=${count::{#count}-1}
fi

# Make output file if it doesn't exist already

if [ ! -d $output ];then
	mkdir $output
fi

# Collapse read files

for file in ${input}/*; do
	if [[ -f $file ]]; then
	
		if [[ ${file:(-3)} == ".gz" ]]; then
			gunzip -c $file | grep -v ">" - | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' > ${count}/$file"_collapsed_counts.txt"
		else
			grep -v ">" $file | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' > ${count}/$file"_collapsed_counts.txt"
		fi

		cut -f2 ${count}/$file"_collapsed_counts.txt" > ${output}/$file"_collapsed.txt"
	fi
done

gzip ${output}/*
gzip ${count}/*
