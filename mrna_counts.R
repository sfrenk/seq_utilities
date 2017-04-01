#!/usr/bin/env Rscript

# Determine log2fold change in reads mapping to mRNA sequences.

library(optparse)

option_list <- list(
    make_option(c("-o", "--output"),
                help = "output filename",
                default = "log2fold_expression_change"),
    make_option(c("-t", "--treatment"),
                help = "column indices of treatment samples"),
    make_option(c("-c", "--control"),
                help = "column indices of control columns"),
    make_option(c("-n", "--no_comparison"),
                help = "report normalized read counts for each sample instead of comparing samples",
                type = "logical",
                action = "store_true",
                default = FALSE),
    make_option(c("-f", "--total_count_file"),
                help = "tab-separated text file with sample names in first column and total mapped reads for each sample in the second")
)

parser <- OptionParser(option_list = option_list, description = c("\n", "Determine log2fold change in reads mapping to mRNA sequences. Input file should be the output of merge_counts.py", "\n", "\n", "NOTE: make sure treatment and control are in the same order to allow comparison of appropriate sample pairs"))

arguments <- parse_args(parser, positional_arguments = 1)

opts <- arguments$options

###############################################################################

# Read in count data

counts <- read.table(arguments$args, header = TRUE)
colnames(counts)[1] <- "locus"

# Remove "_counts.txt" suffix from sample names

colnames(counts) <- sapply(colnames(counts), function(x) gsub("_counts.txt", "", x))

# Add 1 to counts to avoid divide by zero errors

counts[,2:ncol(counts)] <- apply(counts[,2:ncol(counts)], 1, function(x) x + 1)

# Normalize to total number of mapped reads

mapped_reads <- read.table(opts$total_count_file)
colnames(mapped_reads) <- c("sample", "count")

for (i in 2:ncol(counts)){
    mapped <- mapped_reads[mapped_reads$sample == colnames(counts)[i], "count"]
    counts[, i] <- counts[, i]/mapped
}

# Calculate read-count change

if (!opts$no_comparison) {
    
    # Print description of experiment (sanity check)
    
    treatment_cols <- as.numeric(unlist(strsplit(opts$treatment, ",")))
    control_cols <- as.numeric(unlist(strsplit(opts$control, ",")))
    
    print("Treatment:")
    print(colnames(counts)[treatment_cols])
    print("Control: ")
    print(colnames(counts)[control_cols])
    
    comparison_cols <- vector()
    for (i in 1:length(treatment_cols)){
        treatment_col <- treatment_cols[i]
        treatment_name <- colnames(counts)[treatment_col]
        control_col <- control_cols[i]
        control_name <- colnames(counts)[control_col]
        comparison_name <- paste0(treatment_name, "/", control_name)
        # Create a new column for each treatment/control pair
        counts[, comparison_name] <- log2(counts[,treatment_col]/counts[,control_col])
        comparison_cols <- c(comparison_cols, comparison_name)
    }

    write.table(counts[, c("locus", comparison_cols)], opts$output, sep = "\t", quote = FALSE, row.names = FALSE)

} else {
    
    print("No-comparison mode selected - outputting normalized count data")
    write.table(counts, opts$output, sep = "\t", quote = FALSE, row.names = FALSE)
}