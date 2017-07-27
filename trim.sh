#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -n 8

# Trims adapter sequence from fastq files

module add bbmap
module list
###############################################################################
# ARGUMENTS
###############################################################################

# Defaults
min_length=20
paired=false
adapter_file=/nas/longleaf/apps/bbmap/37.25/bbmap/resources/adapters.fa
adapter=""
right=0
left=0

usage="
	USAGE:
		This script maps fastq files and converts sam files to bam. Make sure the following modules are loaded:
		bwa
		samtools

		Trims adapter sequence from fastq files.

		-d/--directory: directory containing fastq files
		-m/--min_length: minimum read length (default: 20)
		-p/--paired: paired end reads
		-a/--adapter: specify adapter sequence (Optional. By default, bbduk uses the adapters.fa file provided with the software)
		-l/--left: bases to remove from the left (5')
		-r/--right: bases to remove from the right (3')
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
		-l|--left)
		left="$2"
		shift
		;;
		-r|--right)
		right="$2"
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

            if [[ $adapter != "" ]]; then

            	
            	bbduk.sh in1=${dir}/${base}_1.fastq.gz in2=${dir}/${base}_2.fastq.gz out1=./trimmed/${base}_1.fastq.gz out2=./trimmed/${base}_2.fastq.gz literal=${adapter} ktrim=r overwrite=true k=23 mink=11 hdist=1 ftl=${left} ftr2=${right} maq=20 minlength=${min_length} tpe tbo
            else

            	bbduk.sh in1=${dir}/${base}_1.fastq.gz in2=${dir}/${base}_2.fastq.gz out1=./trimmed/${base}_1.fastq.gz out2=./trimmed/${base}_2.fastq.gz ref=${adapter_file} ktrim=r overwrite=true k=23 mink=11 hdist=1 maq=20 ftl=${left} ftr2=${right} minlength=${min_length} tpe tbo

            fi
        fi
    else
    	
    	# single end
		base=$(basename $file .fastq.gz)

        echo $(date +"%m-%d-%Y_%H:%M")" Trimming ${base} as single end..."

        if [[ $adapter != "" ]]; then

        	bbduk.sh in=${file} out=./trimmed/${base}.fastq.gz literal=${adapter} ktrim=r overwrite=true k=23 maq=20 mink=11 hdist=1 ftl=${left} minlength=${min_length} ftr2=${right}
        else

        	bbduk.sh in=${file} out=./trimmed/${base}.fastq.gz ref=${adapter_file} ktrim=r overwrite=true k=23 maq=20 mink=11 hdist=1 ftl=${left} minlength=${min_length} ftr2=${right}
        fi
	fi
done
