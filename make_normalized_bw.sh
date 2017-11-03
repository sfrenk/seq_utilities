#!/usr/bin/bash

# Makes normalized BigWig files from a directory of bam files. Tracks are normalized to the mean total mapped reads per library.

module load bedtools ucsctools
module list

usage="
    Make bedgraph files of coverage at specific genome coordinates.

    ARGUMENTS
        -d/--dir
        Directory containing bam files

        -t/--total_mapped
        total_mapped_reads.txt file (tab delimitted file with sample name and total mapped reads)

    "

usage="
    USAGE

    Run this script from the working directory (containing total_mapped_reads.txt and bam/ directory)

    sbatch make_normalized_bg.sh
    "

while [[ $# > 0 ]]
do
    key="$1"
    case $key in
    	-h|--help)
		echo "$usage"
		exit
		;;
        -d|--dir)
        dir="$2"
        shift
        ;;
        -t|--total_mapped)
		total_mapped="$2"
		shift
		;;
    esac
	shift
done

if [[ ! -f $total_mapped ]]; then
	echo "ERROR: please select total_mapped_reads file with the -t/--total_mapped option"
	exit
fi

if [[ ! -d $dir ]]; then
	echo "ERROR: please select bam directory with the -d/--dir option"
	exit
fi

# Get bam files
bamfiles=(${dir}/*.bam)

# Get read counts for each sample

samples=($(cut -f1 ./total_mapped_reads.txt))
counts=($(cut -f2 ./total_mapped_reads.txt))

# Calculate total reads

total=0
for i in ${counts[@]}; do
	total=$((total + i))
done

# Calculate average read count

sample_number=${#counts[@]}
average_count=$((total/sample_number))

printf "number of samples: $sample_number\n"
printf "average count: $average_count\n\n"

# Make normalized bg and bw files

if [ ! -d "bg" ]; then
    mkdir bg
fi

if [ ! -d "bw" ]; then
    mkdir bw
fi


for i in $(seq 0 $((sample_number-1))); do

	sample=${samples[$i]}
	count=${counts[$i]}
	scale_factor=$(bc <<< "scale=3; 10000000/$count")

	printf "processing $sample\n"
	printf "Count: $count Scale factor: $scale_factor\n\n"

	bamfile=$(echo "${bamfiles[@]}" | grep "$sample")

	if [[ ! -f $bamfile ]]; then
		echo "ERROR: could not find bam file for ${sample}"
		exit
	fi

	bedtools genomecov -ibam $bamfile -bg -scale $scale_factor -g ~/proj/seq/WS251/genome/genome.bedtools_genome_file > ./bg/${sample}.bg

	wigToBigWig -clip ./bg/${sample}.bg ~/proj/seq/WS251/genome/genome.bedtools_genome_file ./bw/${sample}.bw

done
