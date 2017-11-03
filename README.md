# Useful scripts to be used in combination with NGS analysis pipelines

## DESeq2.R

Basic script for differential gene expression analysis with DESeq2

## bed_count.sh

Takes a directory of bam files and for each file counts reads mapping to regions in a specified bed file. If a total_mapped_reads.txt file exists, this can be used to calculate normalized read counts

## bed_counts_merge.py

Dependency for bed_count.sh 

## collapse_reads.sh

Used to collapse identical reads (often used in sRNA data to counteract the effect of "jackpotting"). Contains an option to uncollapsed a collapsed file.

## fastqc.sh

Run fastqc on a directory of samples

## go_terms.R

Get GO terms from DESeq2 output file.

## make_bg.sh

Makes bedgraph file showing coverage for a group of bam files at specified loci

## make_normalized_bw.sh

Makes normalized BigWig files from a directory of bam files. Tracks are normalized to the mean total mapped reads per library.

## merge_counts.py

Merges read count files made from individual bam files into a single file that can be loaded straight into DESeq

## rnaseq.R

Contains some basic processing steps for RNA-seq data that don't involve DESeq (eg. normalizing to total mapped reads)

## run_picard_markdup.sh

Run the picard duplicate removal utility on a batch of bam files

## srna_filter

Converts between fasta/fastq/raw file types. Optionally filters for reads based on sequence features (eg. 22 nucleotides long, beginning with a "G") and allows for 3' trimming (eg. to get rid of A-tails). 

## sra_download.sh

Download SRR files from sample IDs (eg. GSM IDs from GEO) and extract fastq files

## trim.sh

Basic script for trimming reads with bbduk
