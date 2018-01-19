#!/usr/bin/bash
#SBATCH --mem 10G
#SBATCH -t 1-0

module add samtools bedtools

bedtools_genome_file="/nas/longleaf/home/sfrenk/proj/seq/WS251/genome/genome.bedtools_genome_file"

usage="
    Make bedgraph files of coverage at specific genome coordinates.

    ARGUMENTS
        -d/--dir
        Directory containing bam files

        -b/--bed
        Bed file containing genome coordinates at which coverage is to be calculated

        -p/--paired_end
        Calculate coverage of paired-end fragments (bedgraph mode only)

        -e/--extend
        Extend reads by this amount (bedgraph mode only)

        -r/--remove_rdna
        Remove reads mapping to rDNA (recommended for most RNA-seq data)

        -g/--genomecov
        Report coverage at each base within regions defined by the bed file.

        -o/--output
        Output filename (default: combined.bg)

    "

# Defaults
dir="."
bed="/nas/longleaf/home/sfrenk/proj/seq/telomere/bed/chrom_ends_2kb.bed"
paired_flag=""
extend_flag=""
remove_rdna=false
output="combined.bg"
genomecov=false

# Arguments

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
		-p|--paired_end)
		paired_flag="-pc"
		;;
		-e|--extend)
		extend_flag="-fs $2"
		shift
		;;
		-r|--remove_rdna)
		remove_rdna=true
		;;
		-b|--bed)
		bedfile="$2"
		shift
		;;
		-g|--genomecov)
		genomecov=true
		;;
		-o|--output)
		output="$2"
		shift
		;;
    esac
	shift
done

if [[ ! -f $bed ]]; then
	echo "ERROR: please select bed file with the -b/--bed option"
	exit
fi

# Remove pre-existing run_parameters_make_bg.txt/total_mapped_reads_rdna_removed.txt file
# Remove the lines below if using this script in parallel mode
if [[ -f run_parameters_make_bg.txt ]]; then rm run_parameters_make_bg.txt; fi
if [[ -f total_mapped_reads_rdna_removed.txt ]]; then rm total_mapped_reads_rdna_removed.txt; fi
if [[ -f total_mapped_reads.txt ]]; then rm total_mapped_reads.txt; fi

# Export run parameters
if [[ ! -f run_parameters_make_bg.txt ]]; then
    printf "$(date +"%m-%d-%Y_%H:%M")\n\nPARAMETERS\n\nPipeline: chip_seq_bowtie\n\nsample directory: ${dir}\nreads_bed: ${reads_bed}\npaired_flag: ${paired_flag}\nextend_flag: ${extend_flag}\nremove rDNA: ${remove_rdna}\n\n" > run_parameters_make_bg.txt

    module list &>> run_parameters_make_bg.txt
    printf "\n" >> run_parameters_make_bg.txt
fi

# Make the necessary directories
if [[ $remove_rdna = true ]] && [[ ! -d bam_no_rdna ]]; then
	mkdir bam_no_rdna
fi

if [[ ! -d bam_subset ]]; then
	mkdir bam_subset
fi

if [[ ! -d bg ]]; then
	mkdir bg
fi

shopt -s nullglob

files=(${dir}/*.bam)

for file in ${files[@]}; do

	base=$(basename $file .bam)
	echo "$base"

	# Remove rDNA

	if [[ $remove_rdna = true ]]; then
		
		printf "Making rDNA- bam for ${base}\n"
		samtools view -b -h -L ~/proj/seq/rdna/rdna.bed -U bam_no_rdna/${base}.bam $file > ${base}_dump.bam
		rm ${base}_dump.bam

		# Calculate total mapped reads for rDNA- bam file

		total_mapped="$(samtools view -c ./bam_no_rdna/${base}.bam)"
		printf "${base}\t${total_mapped}\n" >> total_mapped_reads_rdna_removed.txt

		input_file=bam_no_rdna/${base}.bam

	else
		# Calculate total mapped reads for unfiltered bam file

		total_mapped="$(samtools view -c $file)"
		printf "${base}\t${total_mapped}\n" >> total_mapped_reads.txt

		input_file="$file"
	fi

	# Get reads at genomic regions
	samtools view -bh $input_file -L $bedfile > ./bam_subset/${base}_subset.bam
	samtools index ./bam_subset/${base}_subset.bam

	printf "Making coverage bg file for ${base}\n"

	if [[ $genomecov = true ]]; then

		bedtools genomecov -ibam bam_subset/${base}_subset.bam -dz $paired_flag $extend_flag > bg/${base}.bg
	else

		bedtools coverage -b bam_subset/${base}_subset.bam -a $bedfile -counts -sorted -g $bedtools_genome_file > bg/${base}.bg
	fi

done

# MERGING

printf "Merging files...\n"
array=(bg/*.bg)

if [[ $genomecov = true ]]; then

	# Process genomecov files into temp files so they can be processed

	for file in ${array[@]}; do
		base="$(basename $file)" 

		# Edit position columns to indicate genomic position
		awk '{print $1"\t"$2"\t"($2+1)"\t"$3}' $file > bg/${base}.temp

	done
	array=(bg/*.temp)
fi

bedtools unionbedg -i ${array[@]} -header -names ${array[@]} > $output

if [[ -f bg/*.temp ]]; then
	rm bg/*.temp
fi
