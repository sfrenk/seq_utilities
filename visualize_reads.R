#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))
library(rtracklayer)

### RUN ON CLUSTER ###

parser <- ArgumentParser(description = "Visualize reads in bam files (or other type of data) mapping to a specific region.")

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

parser$add_argument("-a", "--annotation",
                    help = "Comma-spearated list of extra anotation tracks to add (gff3/gtf format)")

parser$add_argument("-n", "--names",
                    help = "Comma separated list of Names of datatracks (optional. By default, name is extracted from filename using regex argument)",
                    default = "")

parser$add_argument("-g", "--gene_annotation",
                    help = "gtf file containing gene annotaion data (Use this option when using a species other than C. elegans)",
                    default = "")

args <- parser$parse_args()

# Get output image type
if (args$output == "") {
    image_type <- "png"
} else if (grepl("\\.png", args$output)){
    image_type <- "png"
} else if (grepl("\\.svg", args$output)){
    image_type <- "svg"
} else{
    print("ERROR: invalid extension for output filename (use .png or .svg)")
    q(save = "no", status = 1)
}

# Check input files
for (file in args$files){
    if (!file.exists(file)){
        args$files <- args$files[args$files != file]
        print(paste0("WARNING: Could not find the file ", file))
    }
}

# Set up data track names
if (args$names != ""){
    data_names <- unlist(strsplit(args$names, ","))
} else {
    data_names <- sapply(args$files, function(x) gsub(args$regex, "\\1", basename(x)))
}

# Set up output name
if (args$output == ""){
    outfile <- paste0(args$chrom, "_", as.character(args$start), "_", as.character(args$end), ".png")
} else {
    outfile <- args$output
}

# Annotation track
options(ucscChromosomeNames=FALSE)
ano_gtf <- args$gene_annotation
if (ano_gtf == ""){ 
    
    print("Attempting to get gene annotation from BioMart")
    
    # There's sometimes a problem with connecting to biomart, so need to keep trying if there is an error
    backup_grtrack <- FALSE
    try_counter <- 0
    while(!exists("mart") & try_counter < 10){
        tryCatch(mart <- useMart(biomart="ensembl", dataset="celegans_gene_ensembl"), error = function(x) {Sys.sleep(10)})
        try_counter <- try_counter + 1
    }
    
    if (!exists("mart")){
        print("Could not connect to Biomart, defaulting to gtf backup")
        ano_gtf <- "/nas/longleaf/home/sfrenk/proj/seq/WS251/genes.gtf"
    }
}

# Axis track
gtrack <- GenomeAxisTrack(fontsize = 10, labelPos = "above")

# Extra annotation tracks
ano_tracks <- list()

ano_files <- unlist(strsplit(args$annotation, ","))
if (length(ano_files) > 0){ 
    for (i in 1:length(ano_files)){
        gff <- import.gff(ano_files[i])
        gff$symbol <- gff$gene_id
        ano_tracks[i] <- GeneRegionTrack(gff, showId = TRUE, name = str_extract(args$annotation[i], "[^/]+$"), background.title = "transparent", fontsize = 10, shape = "arrow")
    }
}

## Plot tracks
#backup_grtrack <- TRUE
make_panel <- function(chrom, start, end){

    # Anotation track (gene models)
    if (ano_gtf == ""){
        grtrack <- BiomartGeneRegionTrack(start, end, mart, chrom, shape = "arrow")
    } else{
        g = as.data.frame(import(ano_gtf))
        g <- g %>% filter(seqnames == chrom)
        grtrack <- GeneRegionTrack(start = start, end = end, chromosome = chrom, rstarts = g$start, rends = g$end, transcript = g$transcript_id, strand = g$strand, name = "", symbol = g$gene_name, showId = TRUE, geneSymbol = TRUE, shape = "arrow")
    }
    
    # Assemble all tracks into one list
    tracks <- list(gtrack, grtrack)
    
    # Add any extra annotation tracks
    n_ano_tracks <- length(ano_tracks)
    if (n_ano_tracks > 0){
        for (i in 1:n_ano_tracks){
            tracks[2+i] <- ano_tracks[i]
        }
    }
    
    ntracks <- length(tracks)
    
    # Data tracks
    n_data_tracks <- length(args$files)
    for (i in 1:n_data_tracks){
        if (grepl(".bam", args$files[i])){
            if (args$ylim > 0){
                tracks[ntracks+i] <- AlignmentsTrack(args$files[i], name = data_names[i], ylim = c(0, args$ylim), type = args$type, frame = TRUE)
            } else {
                tracks[ntracks+i] <- AlignmentsTrack(args$files[i], name = data_names[i], type = args$type, frame = TRUE)
            }
        } else {
            # Bedgraph file
            tracks[ntracks+i] <- DataTrack(args$files[i], genome = "ce11", name = data_names[i], type = "l")
        }
    }

    # Plot
    if (image_type == "png"){
        png(outfile, width = 6, height = 6, units = "in", res = 300)
        plotTracks(tracks, from = start, to = end, chromosome = chrom, baseline = 0, sizes = c(1, 1, rep(0.5, n_ano_tracks), rep(2, n_data_tracks)))
        dev.off()
    } else if (image_type == "svg"){
        svg(outfile, width = 4, height = 4)
        plotTracks(tracks, from = start, to = end, chromosome = chrom, baseline = 0, sizes = c(1, 1, rep(0.5, n_ano_tracks), rep(2, n_data_tracks)))
        dev.off()
    }
}
make_panel(args$chrom, args$start, args$end)
