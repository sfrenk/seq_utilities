#!/usr/bin/env Rscript

### Filter gtf for protein-coding genes ###
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = "Extract protein coding genes from genes.gtf file")

parser$add_argument("file", help = "Input gtf")
parser$add_argument("-o", "--output", help = "Output gtf filename (default: protein_coding_genes.gtf", default = "protein_coding_genes.gtf")

args <- parser$parse_args()

################################################################################

# Get protein coding genes from biomart
print("Getting protein-coding genes from biomart...")
mart <- useMart(biomart="ensembl", dataset="celegans_gene_ensembl")
protein_coding_genes <- getBM(attributes=c("external_gene_name","ensembl_transcript_id"),filters = c("biotype"), values=list("protein_coding"), mart=mart)

# Load gtf
print("Loading GTF file...")
gtf <- read.table(args$file, sep = "\t")
#gtf <- read.table("genes.gtf", sep = "\t")
colnames(gtf) <- c("chrom", "source", "type", "start", "end", "col6", "strand", "col8", "att")

# Filter GTF
# Gene name attribute can cause problems due to multiple names for the same gene
# Therfore, get both gene name and transcript id
print("Filtering GTF file")
gtf$gene_name <- sapply(gtf$att, function(x) gsub(".*gene_name[ :]*([^;]+);.*", "\\1", x))
gtf$transcript_id <- sapply(gtf$att, function(x) gsub(".*transcript_id[ :]*([^;]+);.*", "\\1", x))
gtf <- filter(gtf, gene_name %in% protein_coding_genes$external_gene_name | transcript_id %in% protein_coding_genes$ensembl_transcript_id)

# Output filtered gtf file
write.table(gtf[, 1:9], args$output, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

print("Finished")
