#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))
library(rtracklayer)

### RUN ON CLUSTER ###

parser <- ArgumentParser(description = "Visualize reads in bam files mapping to a specific region.")

parser$add_argument("files",
                    help = "List of bam files to visualize",
                    nargs = "+")

parser$add_argument("-c", "--chrom",
                    help = "Chromosome",
                    required = TRUE,
                    type = "character")

parser$add_argument("-s", "--start",
                    help = "Start coordinate",
                    required = TRUE,
                    type = "integer")

parser$add_argument("-e", "--end",
                    help = "End coordinate",
                    required = TRUE,
                    type = "integer")

parser$add_argument("-o", "--output",
                    help = "Output file name (Default: region coordinates)",
                    type = "character",
                    default = "")

parser$add_argument("-y", "--ylim",
                    help = "Set upper limit for data track y axes (Default: auto)",
                    type = "integer",
                    default = 0)

parser$add_argument("-t", "--type",
                    help = "Type of plot for data tracks (coverage (default) or pileup)",
                    type = "character",
                    default = "coverage")

parser$add_argument("-r", "--regex",
                    help = "regex to convert file basenames to sample names. The regex should include one section in brackets, which will be extracted to obtain the sample name. Default: (.+)\\.bam",
                    type = "character",
                    default = "(.+)\\.bam")

args <- parser$parse_args()

# Set up output name
if (args$output == ""){
    outfile <- paste0(args$chrom, "_", args$start, "_", args$end, ".png")
} else {
    outfile <- args$output
}

# Annotation track
options(ucscChromosomeNames=FALSE)

# There's sometimes a problem with connecting to biomart, so need to keep trying if there is an error
backup_grtrack <- FALSE
try_counter <- 0
while(!exists("mart") & try_counter < 10){
    tryCatch(mart <- useMart(biomart="ensembl", dataset="celegans_gene_ensembl"), error = function(x) {Sys.sleep(10)})
    try_counter <- try_counter + 1
}

if (!exists("mart")){
    print("Could not connect to Biomart, defaulting to gtf backup")
    backup_grtrack <- TRUE
}

gtrack <- GenomeAxisTrack(fontsize = 10, labelPos = "above")

## Plot tracks
backup_grtrack <- TRUE
make_panel <- function(chrom, start, end){

    # Anotation track (gene models)
    if (!backup_grtrack){
        grtrack <- BiomartGeneRegionTrack(start, end, mart, chrom)
    } else{
        g = as.data.frame(import("/nas/longleaf/home/sfrenk/proj/seq/WS251/genes.gtf"))
        g <- g %>% filter(seqnames == chrom)
        grtrack <- GeneRegionTrack(start = start, end = end, chromosome = chrom, rstarts = g$start, rends = g$end, transcript = g$transcript_id, strand = g$strand, name = "", symbol = g$gene_name, showId = TRUE, geneSymbol = TRUE, shape = "smallArrow")
    }
    
    # Assemble all tracks into one list
    tracks <- list(gtrack, grtrack)
    
    # Data tracks
    for (i in 1:length(args$files)){
        if (args$ylim > 0){
            tracks[2+i] <- AlignmentsTrack(args$files[i], name = gsub(args$regex, "\\1", basename(args$files[i])), ylim = c(0, args$ylim), type = args$type)
        } else {
            tracks[2+i] <- AlignmentsTrack(args$files[i], name = gsub(args$regex, "\\1", basename(args$files[i])), type = args$type)
        }
    }

    # Plot
    png(args$output, width = 6, height = 6, units = "in", res = 300)
    plotTracks(tracks, from = start, to = end, chromosome = chrom, baseline = 0, sizes = c(1, 1, rep(2, length(args$files))))
    dev.off()
}

make_panel(args$chrom, args$start, args$end)
