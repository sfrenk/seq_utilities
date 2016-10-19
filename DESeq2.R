#!/usr/bin/env Rscript

###############################################################################
# R script for simple ifferential gene expression analysis with DESeq2
###############################################################################

library("optparse")
suppressPackageStartupMessages(library("DESeq2"))

option_list <- list(
    make_option(c("-o", "--output"),
        help = "output filename (default: DESeq2_results.txt)",
        type = "character",
        default = "DESeq2_results.txt"),
    make_option(c("-c", "--control"),
        help = "column index (or indices) for control sample(s) as a comma-separated list",
        type = "character"),
    make_option(c("-t", "--treatment"),
        help = "column index (or indices) for treatment sample(s) as a comma-separated list",
        type = "character"),
    make_option(c("-f", "--total_counts"),
        help = "if normalizing to total library size (eg. for small RNA libraries) use this option followed by the total_mapped_reads.txt file",
        default = "None")
    )

parser <- OptionParser(option_list = option_list,
                       description = "basic differential expression with DESeq2. Input file should be a count table created using subread featureCounts. Note: the sample columns start at index 1")
arguments <- parse_args(parser, positional_arguments = 1)

opts <- arguments$options
count_file <- arguments$args

output_file <- opts$output
control_cols <- as.numeric(unlist(strsplit(opts$control, ",")))
treatment_cols <- as.numeric(unlist(strsplit(opts$treatment, ",")))

# experimental design
# 1. choose which columns to analyze. Note that the sample columns start at column number 6 from a featureCounts ouput file or column 2 from a merge_counts.py file
sample_cols <- c(control_cols, treatment_cols)

# 2. specify the condition for each sample (eg. "treatment" or "control")
conditions <- c(rep("control", length(control_cols)), rep("treatment", length(treatment_cols)))

################################################################################

# Read in count data.
rawcounts <- read.table(count_file, header = T, row.names = 1)

# Define experimental design:
rawcounts <- rawcounts[,sample_cols]

names <- colnames(rawcounts)
sample_info <- data.frame(condition = conditions)
row.names(sample_info) <- names
print(sample_info)

# Quantifiy differential expression with DESeq
dds_count_table <- DESeqDataSetFromMatrix(countData = rawcounts,
                                          colData = sample_info,
                                          design = ~condition)
if (opts$total_counts != "None"){
    # Create custom size factors
    
    totals <- read.table(opts$total_counts)
    colnames(totals) <- c("sample", "count")
    mapped <- numeric()
    for (i in colnames(rawcounts)){
        # Get total number of mapped reads for each library
        mapped[i] <- totals[totals$sample == i, "count"]
    }
    # Get size factor for each library by dividing the total number of mapped reads by the average for the selected libraries.
    mapped_av <- mean(mapped)
    sizeFactors(dds_count_table) <- sapply(mapped, function(x) x/mapped_av)
}

dds <- DESeq(dds_count_table)
res <- results(dds)
res <- res[order(res$log2FoldChange),]
res <- na.omit(res)

# Output results
write.table(res, file = output_file, sep ="\t", col.names = NA, quote = F)

png(paste0(output_file, "_plot.png"))
plotMA(res)
dev.off()

png(paste0(output_file, "_dispersion.png"))
plotDispEsts(dds)
dev.off()
print(summary(res))