#!/usr/bin/env Rscript

# Determine log2fold change in reads mapping to transposon sequences.

library(argparse)
library(pheatmap)
library(reshape2)
library(ggplot2)

parser <- ArgumentParser(description = "Determine log2fold change in reads mapping to transposon sequences")

parser$add_argument("input", help = "Normalized count table")
parser$add_argument("-t", "--treatment", help = "column indices of treatment samples")
parser$add_argument("-c", "--control", help = "column indices of control columns")
parser$add_argument("-l", "--labels", help = "labels for treatment and control (default: 'treatment' and 'control')", default = c("treatment", "control"))
parser$add_argument("-o", "--output", help = "output filename", default = "rnaseq_results.txt")

args <- parser$parse_args()

###############################################################################

# Read in count data

counts <- read.table(args$input, header = TRUE, row.names = 1)

# Calculate read-count change

comparison_label <- paste0(unlist(strsplit(args$labels, ","))[1], "_VS_", unlist(strsplit(args$labels, ","))[2]) 
    
treatment_cols <- as.numeric(unlist(strsplit(args$treatment, ",")))
control_cols <- as.numeric(unlist(strsplit(args$control, ",")))

# Compare average of treatments with average of controls
        
print("Treatment: Average of:")
print(colnames(counts)[treatment_cols])
print("Control: Average of:")
print(colnames(counts)[control_cols])
        
if (length(treatment_cols) > 1) {
    counts$average_treatment <- rowMeans(counts[, treatment_cols])
} else {
    counts$average_treatment <- counts[, treatment_cols]
}
if (length(control_cols) > 1) {
    counts$average_control <- rowMeans(counts[, control_cols])
} else {
    counts$average_control <- counts[, control_cols]
}

counts[,comparison_label] <- log2(counts$average_treatment/counts$average_control)
counts <- counts[,comparison_label]

write.table(counts, args$output, sep = "\t", quote = FALSE, row.names = FALSE)
