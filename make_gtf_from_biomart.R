#!/usr/bin/env R

# Make genes.gtf file from gene info downloaded from biomart

library(biomaRt)

mart <- useMart(biomart="ensembl", dataset="celegans_gene_ensembl")
gene_info <- getBM(attributes=c("ensembl_transcript_id","transcript_start","transcript_end","ensembl_exon_id","exon_chrom_start","exon_chrom_end","strand","chromosome_name","gene_biotype"),filters = c("ensembl_gene_id","biotype"), values=list(ensembl_id,"protein_coding"), mart=mart)