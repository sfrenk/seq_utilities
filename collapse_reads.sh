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

	if [[ ${input:(-3)} == ".gz" ]]; then
		
		zgrep -v ">" $input | sort | uniq -c | awk '{ print ">read_" NR "_count="$1"\n"$2 }' > ${output}"_temp.txt"
		gzip -c ${output}"_temp.txt" > ${output}
		rm ${output}"_temp.txt"

	else

		grep -v ">" $input | sort | uniq -c | awk '{ print ">read_" NR "_count="$1"\n"$2 }' > ${output}

	fi

else

	# Uncollapse read files

	echo "$input"

	if [[ ${input:(-3)} == ".gz" ]]; then
		
		gunzip -c $input > ${output}"_temp.txt"
		input_file=${output}"_temp.txt"

	else

		input_file=${input}
	fi

	awk 'BEGIN{FS="_count="}{if(/^>/){x=$1;y=$2;z=getline}{for(i=1;i<=y;i++){print x"_"i"\n"$z}}}' $input_file > $output

	if [[ -f ${output}"_temp.txt" ]]; then
		rm ${output}"_temp.txt"
	fi

fi
