#!/usr/bin/bash

# Parallelization script from sequencing pipelines and utilities

pipeline_dir="/nas/longleaf/home/sfrenk/pipelines"
util_dir="/nas/longleaf/home/sfrenk/scripts/util"

usage="Creates a modified copy of a pipeline script in the current directory and sets up the copy for parallelization

	USAGE:

	bash parallelize.sh -p <pipline> -d <fastq directory>

	Currently available pipelines:

	srna
	hisat2
	chip
	telo_cov
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
	telo_cov)
	pipeline_file="${util_dir}/telo_cov.sh"
	;;
esac

# Count number of samples
shopt -s nullglob

if ! [[ $pipeline == "telo_cov" ]]; then
	files=(${dir}/*.fastq.gz)
	pipeline_name="pipeline.sh"
	log_dir="logs"

else
	files=(${dir}/*.bam)
	log_dir="telo_cov_logs"
	pipeline_name="telo_cov_pipeline.sh"

fi

file_number=${#files[@]}

# Create directory for logs
if [[ ! -d $log_dir ]]; then
	mkdir $log_dir
fi

# Modify the pipeline to include:
#	1. SBATCH --array, -o and -e options
#	2. modify the for loop so that only one file is processed per run

sed "2i #SBATCH --array 0-$(($file_number-1))\n#SBATCH -o ./${log_dir}/slurm_%a.out\n#SBATCH -e ./${log_dir}/slurm_%a.err" $pipeline_file | sed 's/for file in \${files\[@\]}/for file in \${files\[\$SLURM_ARRAY_TASK_ID\]}/' > $pipeline_name

printf "Created $pipeline_name for running $pipeline_file with $file_number samples\n"
