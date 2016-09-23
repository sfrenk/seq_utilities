#!/usr/bin/env bash

###############################################################################
# Trim fastq reads with trim galore then convert to raw text format with an optional filtering step
###############################################################################

###############################################################################
# BEFORE RUNNING SCRIPT DO THE FOLLOWING:
###############################################################################

# 1. Make sure modules are loaded:
#       trim_galore (if files are in fastq format)
#       python (default version)

# The location of the following files may have to be modified in this script:
#       small_rna_filter.py

###############################################################################
###############################################################################

usage="
    Trim fastq reads with trim galore then convert to raw text format with an optional filtering step

    USAGE
       step1:   load the following modules: trim_galore python (default version).
       step2:   bash trim_and_filter.sh [options]  

    ARGUMENTS
        -d/--dir
        directory containing read files in fastq.gz format

        -f/--filter 
        filter for 22G RNAs ('g') 21U RNAs ('u') or both (gu) (default = no filtering)

        -k/--keep_fastq
        keep the trimmed fastq files (these are removed by default)

    "
# Set default parameters

filter=""
keep_fastq=false

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
            -f|filter)
            filter="$2"
            shift
            ;;
            -k|--keep_fastq)
            keep_fastq=true
            ;;
    esac
shift
done

# Remove trailing "/" from dir directory if present

if [[ ${dir:(-1)} == "/" ]]; then
    dir=${dir::${#dir}-1}
fi

# parse filter option
case $filter in
    "g")
    filter_opt="-g"
    ;;
    "u")
    filter_opt="-u"
    ;;
    "gu"|"ug")
    filter_opt="-u -g"
    ;;
esac

# Print out loaded modules to keep a record of which software versions were used in this run

modules=$(/nas02/apps/Modules/bin/modulecmd tcsh list 2>&1)
echo "$modules"

# module test

req_modules=("trim_galore" "python")


for i in ${req_modules[@]}; do
    if [[ $modules != *${i}* ]]; then
        echo "ERROR: Please load ${i}"
        exit 1
    fi
done

###############################################################################
###############################################################################

# Make directories
if [ ! -d "trimmed" ]; then
    mkdir trimmed
fi
if [ ! -d "filtered" ]; then
    mkdir filtered
fi

# Start pipeline

echo $(date +"%m-%d-%Y_%H:%M")" Starting pipeline..."

for file in ${dir}/*.fastq.gz; do

    base=$(basename $file .fastq.gz)

    echo $(date +"%m-%d-%Y_%H:%M")" ################ processing ${base} ################"
    
    # Quality/Adaptor trim reads
    trim_galore -o ./trimmed $file

    # Extract 22G and or 21U RNAs or convert all reads to raw format
    python /proj/ahmedlab/steve/seq/util/small_rna_filter.py ${filter_opt} -o ./filtered/${base}.txt ./trimmed/${base}_trimmed.fq.gz

    gzip ./filtered/${base}.txt

done

# Remove the trimmed fastq files unless the -k option is set

if [[ keep_fastq = false ]]; then
    rm -r trimmed
fi
