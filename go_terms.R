#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RDAVIDWebService))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = "Find GO terms for up/downregulated genes in DESeq2 output")

parser$add_argument("file", help = "Input file (results object from DESeq2 as tab-delimitted file)")
parser$add_argument("-c", "--change_cutoff", help = "log2-fold change cutoff for DE genes (default: 2)", default = 2, type = "double")
parser$add_argument("-p", "--pval_cutoff", help = "adjusted p-value cutoff (default: 0.05)", default = 0.05, type = "double")
parser$add_argument("-o", "--output", help = "output basename (default: results)", default = "results")

args <- parser$parse_args()

# Load data
data <- read.table(args$file, sep = "\t", header = TRUE, row.names = 1)
data$gene <- rownames(data)

# Create lists of up and down-regulated genes
up_genes <- filter(data, padj < args$pval_cutoff, log2FoldChange > args$change_cutoff) %>% .$gene
down_genes <- filter(data, padj < args$pval_cutoff, log2FoldChange < -(args$change_cutoff)) %>% .$gene

data$de <- ifelse(data$gene %in% up_genes, "up", ifelse(data$gene %in% down_genes, "dn", ""))

# Set up bioMart
mart <- useMart(biomart="ensembl", dataset="celegans_gene_ensembl")

print("Finding GO terms..")

for (i in c("up", "dn")){
    
    # Connect to DAVID web service
    david <- DAVIDWebService$new(email="sfrenk@unc.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
    
    # Convert gene names to Wormbase IDs for DAVID
    results <- getBM(attributes = "wormbase_gene", 
                     filters = "external_gene_name", 
                     values = data %>% filter(de == i) %>% .$gene,
                     mart = mart)

    # Define gene list for DAVID
    result <- addList(david, results$wormbase_gene, idType = "WORMBASE_GENE_ID", listName = "list", listType = "Gene")
    # See proportion of genes mapped by david
    print(result)
    
    # Specify annotation categories
    setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
    
    # Get functional annotation chart
    #FuncAnnotChart <- getFunctionalAnnotationChart(david)
    getFunctionalAnnotationChartFile(david, paste0(args$output, "_", i, ".txt"))
    
    # Get functional annotation cluster
    termCluster <- getClusterReport(david, type="Term")

    print(head(summary(termCluster)))
}