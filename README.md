# Useful scripts to be used in combination with NGS analysis pipelines

## srna_filter

Optionally filters for 22G and/or 21U RNA sequences from fastq/fasta/raw file and outputs a raw sequence file

## sra_lookup.sh

Find SRR file names from sample IDs (eg. GSM IDs from GEO)

## batch_sra_download.py

Reads in output from sra_lookup.sh and downloads SRA files then extracts adn gzips the fastq file(s)

## merge_counts.py

Merges read count files made from individual bam files into a single file that can be loaded straight into DESeq

## remove_rdna.R

Remove reads mapping to rDNA loci from count table

## DESeq2.R

Basic script for differential gene expression analysis with DESeq2

## trim_and_filter.sh

Trims reads with trim_galore then converts reads to raw format with an optional filtering step

## a_trimmer

Remove 3' A nucleotides from raw format reads (required for some small RNA libraries)
