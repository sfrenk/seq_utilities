#!/usr/bin/env bash
#SBATCH -t 5-0
#SBATCH --mem 30G

module add picard samtools

###############################################################################
# ARGUMENTS
###############################################################################

# Defaults
output="."

usage="
	USAGE:
		Remove PCR duplicates from MMP bam files.

		Make sure the following modules are loaded:
		
		samtools

		-b/--bamfile: input bam file
		-o/--output: destination directory for duplicate-removed bam files

"
while [[ $# > 0 ]]
do
    key="$1"
    case $key in
        -h|--help)
        echo "$usage"
        exit
        ;;
		-o|--output)
		output="$2"
		shift
		;;
		-b|--bamfile)
		bamfile="$2"
		shift
		;;
    esac
shift
done


# Check input file
if [ ! -f $bamfile ]; then
	echo "ERROR: invalid bam file"
	exit 1
fi

# Check output directory
if [[ ! -d $output ]]; then
	echo "ERROR: Invalid output directory path"
	exit 1
fi

# Remove trailing "/" from directories if present
if [[ ${output:(-1)} == "/" ]]; then
    output=${output::${#output}-1}
fi

###############################################################################

echo "Processing "$bamfile

base=$(basename $bamfile)

java -jar /nas/longleaf/apps/picard/2.2.4/picard-tools-2.2.4/picard.jar MarkDuplicates INPUT=$bamfile OUTPUT=${output}/$base METRICS_FILE=${output}/picard_metrics.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true

samtools index ${output}/$base
