#!/usr/bin/bash
#SBATCH -t 5-0
#SBATCH --mem 4G

module add sratoolkit edirect

###############################################################################
# Download SRA data
###############################################################################


###############################################################################
# VARIABLES
###############################################################################

# Specify directory where SRA files are downloaded to 

sra_dir="/nas/longleaf/home/sfrenk/ncbi/public/sra"


###############################################################################
# CLEANUP
###############################################################################


# Make sure any partly-downloaded SRA files get removed if the script exists. This prevents the "lock" error.

function remove_sra {
	echo "Cleaning up SRA directory"
	rm -r ${sra_dir}/*
} 

trap remove_sra EXIT


###############################################################################
# ARGUMENTS
###############################################################################

# Defaults
output="."

usage="
	USAGE:
		Download SRA data.

		sra_download -o <output_directory> -i <input_file>

		-i/--input: two column tab-separated table with ID (eg. published sample name/GSM ID) in the first column and sample name in the second. Example:

			GSM336052	hermaphrodite_embryo

		-o/--output: output directory name (current directory by default)
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

if [[ -z $input ]]; then
	echo "ERROR: please specify input file with -i option"
	exit 1
fi

if [[ -z $output ]]; then
	echo "ERROR: please specify output file name with -o option"
	exit 1
fi

# Remove log file if it already exists

log_file=${output}/sra_log.txt

if [[ -f $log_file ]]; then
	rm $log_file
fi

##############################################################################
##############################################################################

# Remove any illegal/problematic characters from sample names
# Have to convert to csv format because array will split on tab 

sed 's/[ -]/_/g' $input | sed 's/[\(\),\.:;#]//g' | sed 's/[\t]/,/g' > ${input}.temp
readarray samples < ${input}.temp

for line in ${samples[@]}; do

	# Extract id and name

	id=$(echo "$line" | cut -d"," -f 1)
	name=$(echo "$line" | cut -d"," -f 2)

	# Get SRR file names from SRA

	readarray -t srr_files <<<"$(esearch -db sra -query $id | efetch --format runinfo | cut -d ',' -f 1 | grep SRR)"

	# Output sample and SRR details to log file

	echo -e "${name}\t${srr_files[@]}\n" >> $log_file
	
	# Download and srr file(s) and convert to fastq.gz

	if [[ "${#srr_files[@]}" > 1 ]]; then

		# Multiple SRR files for one sample

		mkdir ${output}/${name}
		for srr in "${srr_files[@]}"; do

			prefetch $srr
			fastq-dump --outdir ${output}/${name} --split-files --gzip ${sra_dir}/${srr}.sra
			rm ${sra_dir}/${srr}.sra

		done
	else

		# One SRR file for one sample
		srr="$srr_files"
		prefetch $srr
		fastq-dump --outdir ${output} --split-files --gzip ${sra_dir}/${srr}.sra
		rm ${sra_dir}/${srr}.sra

		# Rename file

		if [[ -f "${output}/${srr}_2.fastq.gz" ]]; then
				
			# paired end

			mv ${output}/${srr}_1.fastq.gz ${output}/${name}_1.fastq.gz
			mv ${output}/${srr}_2.fastq.gz ${output}/${name}_2.fastq.gz

		else

			# single end

			mv ${output}/${srr}_1.fastq.gz ${output}/${name}.fastq.gz

		fi
	fi

done

rm ${input}.temp
