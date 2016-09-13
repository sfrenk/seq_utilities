#!/usr/bin/env bash

###############################################################################
# Find SRR file names from sample/GSM ids.
###############################################################################

###############################################################################
# ARGUMENTS
###############################################################################

usage="
	USAGE:
		Find SRR file names from sample/GSM ids. output file consists of two columns containing name and SRR ID respectively for each sample.

		sra_lookup -o <output_file> -i <input_file>

		-i/--input: two column tab-separated table with ID (eg. published sample name/GSM ID) in the first column and sample name in the second
		-o/--output: output filename
"

if [ -z "$1" ]; then
    echo "$usage"
    exit
fi

while [[ $# > 0 ]]
do
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
    esac
shift
done

# Make sure the required options have been specified

if [ -z $input ]; then
	echo "ERROR: please specify input file with -i option"
	exit 1
fi

if [ -z $output ]; then
	echo "ERROR: please specify output file name with -o option"
	exit 1
fi

# Remove any pre-existing file with the same name as the output file

if [ -e $output ]; then
	rm $output
fi

##############################################################################
##############################################################################

# Make a list of samples, whith each sample assigned to its ID

samples=()

while read line; do
	# Extract id and name from the file and remove any illegal characters
	id=`echo "$line" | cut -f 1 | sed s/" "/"_"/g`
	name=`echo "$line" | cut -f 2 | sed s/" "/"_"/g | sed s/[\(\),\.]//g`
	samples+=`echo "${id}","${name}"`" "
done < $input

# Get the srr file names for each sample

for sample in ${samples[@]}; do
	
	id=`echo "$sample" | cut -d"," -f 1`
	name=`echo "$sample" | cut -d"," -f 2`

	echo $name

	# Use the Entrez command line utility to search for srr files based on the sample ID
	
	srr_files=(`esearch -db sra -query ${id} | efetch --format runinfo | cut -d ',' -f 1 | grep SRR`)

	# If more than one srr files exist for a single sample, make a comma-separated list of the file names

	srr_files=$(printf ",%s" "${srr_files[@]}")
	srr_files=${srr_files:1}
	echo $srr_files

	# Output the results as a tsv file with sample name in the first column and srr file name in the second

	printf ${name}"\t"${srr_files}"\n" >> ${output}

done
