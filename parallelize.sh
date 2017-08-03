#!/usr/bin/bash

# Parallelization script from sequencing pipelines

pipeline_dir="/nas/longleaf/home/sfrenk/pipelines"

usage="Creates a modified copy of a pipeline script in the current directory and sets up the copy for parallelization

	USAGE:

	bash parallelize.sh -p <pipline> -d <fastq directory>

	Currently available pipelines:

	srna
	hisat2
	chip
	"

# Parse command line parameters

if [ -z "$1" ]; then
    echo "$usage"
    exit
fi

while [[ $# > 0 ]]
do
    key="$1"
    case $key in
        -d|--dir)
        dir="$2"
        shift
        ;;
        -p|--pipeline)
        pipeline="$2"
        shift
        ;;
        -h|--help)
        echo "$usage"
    	exit
    	;;
    esac
    shift
done

# Check arguments
args=(dir, pipeline)

for i in ${args[@]}; do

	if [[ -z $i ]]; then
		printf "ERROR: Please supply $i argument"
		exit
	fi
done

# Select pipeline
case $pipeline in
	srna)
	pipeline_file="${pipeline_dir}/bowtie_srna.sh"
	;;
	chip)
	pipeline_file="${pipeline_dir}/chip_bowtie.sh"
	;;
	hisat2)
	pipeline_file="${pipeline_dir}/hisat2.sh"
	;;
esac

# Count number of samples
shopt -s nullglob
fastq=(${dir}/*.fastq.gz)
fastq_number=${#fastq[@]}

# Create directory for logs
if [[ ! -d logs ]]; then
	mkdir logs
fi

# Modify the pipeline to include:
#	1. SBATCH --array, -o and -e options
#	2. modify the for loop so that only one file is processed per run

sed "2i #SBATCH --array 0-$(($fastq_number-1))\n#SBATCH -o ./logs/slurm_%a.out\n#SBATCH -e ./logs/slurm_%a.err" $pipeline_file | sed 's/for file in \${files\[@\]}/for file in \${files\[\$SLURM_ARRAY_TASK_ID\]}/' > pipeline.sh

printf "Created pipeline.sh for running $pipeline_file with $fastq_number samples\n"
