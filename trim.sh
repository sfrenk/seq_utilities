#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -n 8

# Trims adapter sequence from fastq files

module add bbmap

###############################################################################
# ARGUMENTS
###############################################################################

# Defaults
min_length=20
paired=false
adapter_file=/nas/longleaf/apps/bbmap/37.25/bbmap/resources/adapters.fa
adapter="ref=${adapter_file}"
right=0
left=0
n_mink=11
hdist=1
qual=0
k_option=23

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
		-n/--n_mink: value for mink (minimum kmer length to trim at 3' end. Default: 11)
		-h/--hdist: hamming distance (number of mismatches allowed for adapters. Default: 1)
		-q/--qual: minimum read quality to be kept (default:0)
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
		adapter="literal=$2"
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
		-n|--n_mink)
		n_mink="$2"
		shift
		;;
		-h|--hdist)
		hdist="$2"
		shift
		;;
		-q|--qual)
		qual="$2"
		shift
		;;
    esac
	shift
done

# Check directory
if [[ ! -d $dir ]]; then
	echo "ERROR: Invalid fastq directory path"
	exit 1
fi

# If using a short adapter sequence, adjust kmer size
if [[ ${adapter:0:7} == "literal" ]]; then
	adapter_seq=${adapter#literal=}
	
	if [[ ${#adapter_seq} < $k_option ]]; then
		k_option="${#adapter_seq}"
	fi
fi


# Get adapter sequences for library type
#if [ srna = false ]; then
#	adapter=/nas02/apps/bbmap-36.64/bbmap/resources/truseq.fa.gz
#else
#	adapter=/nas02/apps/bbmap-36.64/bbmap/resources/truseq_rna.fa.gz
#fi

params_file="trim_params.${SLURM_JOB_ID}"

module list &> $params_file

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

            trim_cmd="bbduk.sh in1=${dir}/${base}_1.fastq.gz in2=${dir}/${base}_2.fastq.gz out1=./trimmed/${base}_1.fastq.gz out2=./trimmed/${base}_2.fastq.gz $adapter ktrim=r overwrite=true k=$k_option mink=${n_mink} hdist=$hdist ftl=${left} ftr2=${right} minlength=${min_length} tpe tbo minavgquality=$qual"

            printf "\n${trim_cmd}\n" >> $params_file
            $trim_cmd
       
        fi

    else

    	# single end
		base=$(basename $file .fastq.gz)

        echo $(date +"%m-%d-%Y_%H:%M")" Trimming ${base} as single end..."


        trim_cmd="bbduk.sh in=${file} out=./trimmed/${base}.fastq.gz $adapter ktrim=r overwrite=true k=$k_option mink=$n_mink hdist=$hdist ftl=$left minlength=$min_length ftr2=$right minavgquality=$qual"

        printf "\n${trim_cmd}\n" >> $params_file
        $trim_cmd

	fi
done
