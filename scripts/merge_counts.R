#!/usr/bin/env Rscript


# Script to merge output count files from HTSeq

### Import packages ###
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(rtracklayer)
})

### Command arguments ###
args = commandArgs(trailingOnly = TRUE)

countfile_dir = args[1] #Count file  <- 'path/to/counts/'
sample_names = as.character(args[2]) #Sample names <- 'sample1,sample2,sample3'
gtf_file = args[3] # GTF file also use for mapping
outfile = args[4] #Output file <- 'path/to/counts_raw.tsv'

### Test Zone ###
# countfile_dir <- "~/Library/CloudStorage/OneDrive-Aarhusuniversitet/PhD project SLA/Projects/rna_seq_snake_pipeline/temp/count/"
# sample_names <- 'SRR14870728,SRR14870729'
# gtf_file <- "~/Library/CloudStorage/OneDrive-Aarhusuniversitet/PhD project SLA/Projects/rna_seq_snake_pipeline/temp/genome/Homo_sapiens.GRCh38.115.gtf"
# outfile <- "~/Library/CloudStorage/OneDrive-Aarhusuniversitet/PhD project SLA/Projects/rna_seq_snake_pipeline/temp/count/counts_raw.csv"

# Open first count file and make data for count data
sample_idx <- strsplit(sample_names, ',')[[1]]
sample_first <- sample_idx[1]
count_df <- read_tsv(paste0(countfile_dir, sample_first, '.tsv'), col_names = TRUE) %>% as.data.frame()
rownames(count_df) <- count_df[,1]
row_order <- rownames(count_df)
colnames(count_df) <- c('Gene ID', sample_first)

# Add data from other files
if (length(sample_idx) > 1) {
  for (id in sample_idx[2:length(sample_idx)]) {
    # Open count file and order to fit count dataframe
    count_file <- read_tsv(paste0(countfile_dir, id, '.tsv'), col_names = TRUE) %>% as.data.frame()
    rownames(count_file) <- count_file[,1]
    count_file_ordered <- count_file[row_order,]
    counts <- count_file_ordered[,2]
    
    # Merge with count df
    count_df[,id] <- counts
  }
}

# Add gene names and biotype to dataframe
gtf <- import(gtf_file)
gtf_genes <- subset(gtf, gtf$type == 'gene')

geneidx_df <- data.frame(
  'gene_id' = gtf_genes$gene_id,
  'gene_name' = gtf_genes$gene_name,
  'gene_biotype' = gtf_genes$gene_biotype
)
rownames(geneidx_df) <- geneidx_df$gene_id
rm(gtf, gtf_genes)

count_df[row_order,'Gene name'] <- geneidx_df[row_order,'gene_name']
count_df[row_order,'Gene biotype'] <- geneidx_df[row_order,'gene_biotype']


# Save count data as CSV file
write.csv(count_df[,c('Gene ID', 'Gene name', 'Gene biotype', sample_idx)], 
          file = outfile)
