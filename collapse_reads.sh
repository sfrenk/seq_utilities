#!/usr/bin/env bash
#SBATCH -t 1-0
#SBATCH --mem 4G

###############################################################################
# Collapse/uncollapse reads in raw/fasta format to unique sequences.
###############################################################################

usage="
    Collapse reads in raw/fasta format to unique sequences. Output consists of raw reads.

    ARGUMENTS
        -i/--input		input file
        -o/--output		output file
        -u/--uncollapse uncollapse reads     
    "

# Defaults:
uncollapse=false

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
		-u|--uncollapse)
		uncollapse=true
		;;
		-h|--help)
		echo "$usage"
		exit
		;;
	esac
	shift
done

# Check input file

if [[ ! -f $input ]]; then
	echo "ERROR: Invalid input file"
	exit 1
fi

###############################################################################
###############################################################################

if [[ $uncollapse = false ]]; then

	# Collapse read file

	if [[ ${file:(-3)} == ".gz" ]]; then
		
		zgrep -v ">" $input | sort | uniq -c | awk '{ print ">read_" NR "_count="$1"\n"$2 }' > ${output}"_temp.txt"
		gzip -c ${output}"_temp.txt" > ${output}
		rm ${output}"_temp.txt"

	else

		grep -v ">" $input | sort | uniq -c | awk '{ print ">read_" NR "_count="$1"\n"$2 }' > ${output}

	fi

else

	# Uncollapse read files

	if [[ ${file:(-3)} == ".gz" ]]; then

		gunzip -c $input > ${output}"_temp.txt"
		input_file=${output}"_temp.txt"

	else

		input_file=${input}
	fi

	# Remove output file if it already exists

	if [[ -f $output ]]; then
		rm $output
	fi

	while read line; do

		if [[ ${line:0:1} == ">" ]]; then

			header="$line"
			count=$(echo "$line" | sed -r 's/.*count=([0-9]+).*/\1/g')
			read line
			seq="$line"

			while [[ $count > 0 ]]; do

				printf "%s\n" "$header" >> ${output}
				printf "%s\n" "$seq" >> ${output}
				
				let count-=1

			done

		fi

	done < $input_file

	if [[ -f ${output}"_temp.txt" ]]; then
		rm ${output}"_temp.txt"
	fi

fi
