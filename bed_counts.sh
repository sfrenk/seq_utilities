#!/usr/bin/bash
#SBATCH -t 1-0
#SBATCH --mem 10G

module load python bedtools

# For each bam file in a directory, count the number of reads mapping to regions defined in a bed file using bedtools coverage. Requires bed_counts_merge.py

###############################################################################
# VARIABLES
###############################################################################

# path to merge_counts.py script

mc_path="/proj/ahmedlab/steve/seq/util/bed_counts_merge.py"


usage="USAGE:
			For each bam file in a directory, count the number of reads mapping to regions defined in a bed file. Make sure python is loaded!

			bed_counts.sh -r <regions bed file> -d <bam directory>

			Other options:

			-t --total: Specify for total_mapped_reads.txt file for normalization (optional)

"

# Set defaults
total="not_specified"

# Print help if no arguments given

if [[ $# == 0 ]]; then
	echo "$usage"
	exit
fi

# Parse arguments

while [[ $# > 0 ]]; do
	key="$1"
	case $key in
		-d|--dir)
		dir="$2"
		shift
		;;
		-r|--regions)
		regions="$2"
		shift
		;;
		-t|--total)
		total="$2"
		shift
		;;
	esac
	shift
done

# Display error message if one of the arguments is missing

if [ -z $dir ]; then
	echo "ERROR: please specify bam directory with -d option"
	exit 1
fi
if [ -z $regions ]; then
	echo "ERROR: please specify bed file containing regions with -r option"
	exit 1
fi

# Remove trailing slash from directory if it exists

if [[ ${dir:(-1)} == "/" ]]; then
	dir=${dir::${#dir}-1}
fi

# Process all bam files in bam directory

for file in ${dir}/*; do
	if [[ ${file:(-4)} == ".bam" ]]; then
		base=$(basename $file .bam)
		bedtools coverage -counts -b $file -a $regions > ${base}.bed
	fi
done

if [[ $total == "not_specified" ]]; then
	python $mc_path .
else
	python $mc_path -t $total .
fi
