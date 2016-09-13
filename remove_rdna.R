#!/usr/bin/env Rscript

library("optparse")

option_list <- list(
    make_option(c("-o", "--output"),
                help = "output filename",
                default = "counts_rdna_removed.txt",
                type = "character")
)

parser <- OptionParser(option_list = option_list, description = "remove rDNA genes from counts.txt file")
arguments <- parse_args(parser, positional_arguments = 1)
opts <- arguments$options

data <- read.table(arguments$args, sep = "\t")

# rdna_genes.txt contains a list of rDNA genes on chr I and V which need to be removed from the count file
rdna_genes <- read.table("/home/sfrenk/Documents/Resources/Seq/rdna_genes.txt")$V1
data <- data[!(data$V1 %in% rdna_genes),]

write.table(data, opts$output, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)