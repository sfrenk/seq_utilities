library(dplyr)
library("RDAVIDWebService")
library(biomaRt)

data <- read.table("~/Documents/Sequencing2/nrde_1021_FINAL/nrde_late_vs_nrde_early/DESeq2_results.txt")
data$gene <- rownames(data)

all_genes <- data$gene

up_genes <- filter(data, padj < 0.05, log2FoldChange > 0) %>% select(gene) %>% .$gene
down_genes <- filter(data, padj < 0.05, log2FoldChange < 0) %>% select(gene) %>% .$gene

david<-DAVIDWebService$new(email="sfrenk@unc.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

mart <- useMart(biomart="ensembl", dataset="celegans_gene_ensembl")

results <- getBM(attributes = "wormbase_gene", 
                 filters = "external_gene_name", 
                 values = up_genes,
                 mart = mart)

head(results)

result <- addList(david, results$wormbase_gene, idType = "WORMBASE_GENE_ID", listName = "list", listType = "Gene")
termCLuster <- getClusterReport(david, type="Term")

