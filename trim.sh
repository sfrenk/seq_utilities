#!/usr/bin/env bash

# Trims adapter sequence from fastq files


###############################################################################
# ARGUMENTS
###############################################################################

# Defaults
min_length=20
paired=false
adapter=""

usage="
	USAGE:
		This script maps fastq files and converts sam files to bam. Make sure the following modules are loaded:
		bwa
		samtools

		Trims adapter sequence from fastq files.

		-d/--directory: directory containing fastq files
		-m/--min_length: minimum read length (default: 20)
		-p/--paired: paired end reads
		-a/--adapter: specify adapter sequence (optional)
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
		dir="$2"
		shift
		;;
		-m|--min_length)
		min_length="$2"
		shift
		;;
		-p|--paired)
		paired=true
		;;
		-a|--adapter)
		adapter="-a "$2" "
		shift
		;;
    esac
shift
done

# Check directory
if [ ! -d $dir ]; then
	echo "ERROR: Invalid fastq directory path"
	exit 1
fi

# Remove trailing "/" from fastq directory if present
if [[ ${dir:(-1)} == "/" ]]; then
    dir=${dir::${#dir}-1}
fi

# Get adapter sequences for library type
#if [ srna = false ]; then
#	adapter=/nas02/apps/bbmap-36.64/bbmap/resources/truseq.fa.gz
#else
#	adapter=/nas02/apps/bbmap-36.64/bbmap/resources/truseq_rna.fa.gz
#fi


###############################################################################
# VARIABLES
###############################################################################

#  Set the command below to display software versions used during the run
modules=$(/nas02/apps/Modules/bin/modulecmd tcsh list 2>&1)

###############################################################################
# MODULE TEST
###############################################################################

req_modules=("trim_galore")

for i in ${req_modules[@]}; do
	if [[ $modules != *${i}* ]]; then
		echo "ERROR: Please load ${i}"
		exit 1
	fi
done

###############################################################################

echo "$modules"

if [ ! -d "trimmed" ]; then
	mkdir "trimmed"
fi

echo "$(date +"%m-%d-%Y_%H:%M") Starting pipeline"

for file in ${dir}/*.fastq.gz; do

    if [[ $paired = true ]]; then
            
        # paired end

        if [[ ${file:(-11)} == "_1.fastq.gz" ]]; then
        
            Fbase=$(basename $file .fastq.gz)
            base=${Fbase%_1}

            echo $(date +"%m-%d-%Y_%H:%M")" Trimming ${base} as paired end..."

            trim_galore ${adapter}--dont_gzip --length $min_length -o ./trimmed --paired ${dir}/${base}_1.fastq.gz ${dir}/${base}_2.fastq.gz
        fi
    else
    	
    	# single end
		base=$(basename $file .fastq.gz)


        echo $(date +"%m-%d-%Y_%H:%M")" Trimming ${base} as single end..."

        trim_galore ${adapter}--dont_gzip --length $min_length -o ./trimmed ${dir}/${base}.fastq.gz

	fi
done
