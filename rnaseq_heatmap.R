#!/usr/bin/env Rscript

# Make heatmap based on read count normalized by number of genomic mapped reads

suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(argparse))

################################################################################
# VARIABLES #

# Transposons to skip
skip <- "CELE45"

# Minimum number of counts required for a locus to be kept
min_counts <- 20
################################################################################

parser <- ArgumentParser(description = "Make heatmap based on read count normalized by number of genomic mapped reads")

parser$add_argument("input", help = "count table")
parser$add_argument("-t", "--total", help = "total_mapped_reads.txt file", required = TRUE)
parser$add_argument("-o", "--output", help = "output basename (default: 'rnaseq')", default = "rnaseq")
parser$add_argument("-m", "--max", help = "maximum number of loci to show", default = 0)

args <- parser$parse_args()

# Read in count data
count_table <- read.table(args$input, header = TRUE, sep = "\t")
counts <- count_table[,-1]
rownames(counts) <- count_table[,1]
counts <- counts[!rownames(counts) %in% skip,]

# Remove "_counts.txt" suffix from sample names
colnames(counts) <- sapply(colnames(counts), function(x) gsub("_counts.txt", "", x))

# Change any symbols that may cause problems
colnames(counts) <- sapply(colnames(counts), function(x) gsub("[.]| |-", "_", x))

# Remove any low abundance transcripts (<10 transcripts in all samples)
counts$max <- apply(counts, 1 , max)
counts <- counts[counts$max >= min_counts, ]
counts <- counts[order(counts$max), ]
counts$max <- NULL
counts <- counts + 1

# Keep only the loci with the highest counts
if (args$max > 0) {
    counts <- head(counts, args$max)
}

# Normalize each sample to total number of mapped reads when alligning to genome with hisat
mapped_reads <- read.table(args$total)
colnames(mapped_reads) <- c("sample", "count")

# Change any symbols that may cause problems
mapped_reads$sample <- sapply(mapped_reads$sample, function(x) gsub("[.]| |-", "_", x))

for (i in 1:ncol(counts)){
    mapped <- mapped_reads[mapped_reads$sample == colnames(counts)[i], "count"]
    counts[, i] <- counts[, i]/mapped
}

# Output normalized counts
write.table(counts, paste0(args$output, "_normalized_counts.txt"), sep = "\t", quote = FALSE)

# Plot heatmap
counts <- as.matrix(counts)
counts <- scale(counts)
pheatmap(counts, filename = paste0(args$output, "_heatmap.png"))
