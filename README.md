# Useful scripts to be used in combination with NGS analysis pipelines

## srna_filter

Optionally filters for 22G and/or 21U RNA sequences from fastq/fasta/raw file and outputs a raw sequence file

## sra_download.sh

Download SRR files from sample IDs (eg. GSM IDs from GEO) and extract fastq files

## merge_counts.py

Merges read count files made from individual bam files into a single file that can be loaded straight into DESeq

## remove_rdna.R

Removes reads mapping to rDNA loci from count table

## DESeq2.R

Basic script for differential gene expression analysis with DESeq2

## trim.sh

Basic script for trimming reads

## bed_count.sh

Takes a directory of bam files and for each file counts reads mapping to regions in a specified bed file. If a total_mapped_reads.txt file exists, this can be used to calculate normalized read counts

## bed_counts_merge.py

Dependency for bed_count.sh 

## collapse_reads.sh

Used to collapse identical reads (often used in sRNA data to counteract the effect of "jackpotting"). Contains an option to uncollapsed a collapsed file.

## run_picard_markdup.sh

Run the picard duplicate removal utility on a batch of bam files

## rnaseq.R

Contains some basic processing steps for RNA-seq data that don't involve DESeq (eg. normalizing to total mapped reads)
