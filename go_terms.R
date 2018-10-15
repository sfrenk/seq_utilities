#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RDAVIDWebService))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = "Find GO terms for a list of genes")

parser$add_argument("file", help = "Input file (list of genes, one gene per line)")
parser$add_argument("-o", "--output", help = "output filename (default: go_terms.txt))", default = "go_terms.txt")
parser$add_argument("-c", "--convert", help = "Use this option to convert gene IDs from gene name to wormbase ID", default = FALSE, action = "store_true")

args <- parser$parse_args()

# Load data
genes <- read.table(args$file)$V1

# Convert gene IDs if necessary
if (args$convert){
    # Set up bioMart
    # There's sometimes a problem with connecting to biomart, so need to keep trying if there is an error
    try_counter <- 0
    while(!exists("mart") & try_counter < 10){
        tryCatch(mart <- useMart(biomart="ensembl", dataset="celegans_gene_ensembl"), error = function(x) {Sys.sleep(10)})
        try_counter <- try_counter + 1
    }
    
    if (!exists("mart")){
        print("ERROR: Could not connect to Biomart")
        quit(save = "no", status = 1, runLast = FALSE)
    }
    
    # Convert gene names to Wormbase IDs for DAVID
    results <- getBM(attributes = "wormbase_gene", 
                     filters = "external_gene_name", 
                     values = genes,
                     mart = mart)
    gene_list <- results$wormbase_gene
} else{
    gene_list <- genes
}

print("Finding GO terms..")

# Connect to DAVID web service
david <- DAVIDWebService$new(email="sfrenk@unc.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

# Define gene list for DAVID
go_result <- addList(david, gene_list, idType = "WORMBASE_GENE_ID", listName = "list", listType = "Gene")
# See proportion of genes mapped by david
print(paste0("Proportion of genes recognized by DAVID: ", go_result[[1]]))

# Specify annotation categories
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))

# Get functional annotation chart
#FuncAnnotChart <- getFunctionalAnnotationChart(david)
getFunctionalAnnotationChartFile(david, args$output)

# Get functional annotation cluster
termCluster <- getClusterReport(david, type="Term")
