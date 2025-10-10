#!/usr/bin/env Rscript


# Script to normalise counts usinge DeSeq2

### Import packages ###
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
})

### Command arguments ###
args = commandArgs(trailingOnly = TRUE)

countfile = args[1] #Count file  <- 'path/to/counts_raw.tsv'
out_dir = args[2] #Output file <- 'path/to/output'

### Test Zone ###
# countfile <- "~/Library/CloudStorage/OneDrive-Aarhusuniversitet/PhD project SLA/Projects/rna_seq_snake_pipeline/temp/count/counts_raw.csv"
# out_dir <- "~/Library/CloudStorage/OneDrive-Aarhusuniversitet/PhD project SLA/Projects/rna_seq_snake_pipeline/temp/count/"


# Load counts
counts <- read.csv(file = countfile,
                   header = TRUE)

# Remove non-genes counts
counts_fltr <- counts[1:(nrow(counts)-5),]

#### FOR TESTING ####
# counts_fltr$SRR14870729 <- sample(counts_fltr$SRR14870729)

row_order <- counts_fltr$Gene.ID
rownames(counts_fltr) <- counts_fltr$Gene.ID

count_values <- counts_fltr[,5:ncol(counts_fltr)]
rownames(count_values) <- counts_fltr$Gene.ID

# Load design factors
coldata <- data.frame(
  row.names = colnames(count_values),
  condition = factor(rep('A', length(colnames(count_values))))
)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = count_values, colData = coldata, design = ~1)

# Add size factors to dataset
dds <- estimateSizeFactors(dds)

# Get normalized counts
norm_counts <- counts(dds, normalized = TRUE) %>% as.data.frame()
norm_counts[,'Gene ID'] <- rownames(count_values)
norm_counts[row_order, 'Gene name'] <- counts_fltr[row_order, 'Gene.name']
norm_counts[row_order, 'Gene biotype'] <- counts_fltr[row_order, 'Gene.biotype']

# Get size factors
sf <- sizeFactors(dds)

# Save normalised counts as CSV file
write.csv(norm_counts[,c('Gene ID', 'Gene name', 'Gene biotype', colnames(count_values))], 
          file = paste0(out_dir, 'counts_normalised_deseq2.csv'))

# Save size factors
sf_df <- as.data.frame(sf)
colnames(sf_df) <- 'Size Factor'
write.csv(sf_df, 
          file = paste0(out_dir, 'deseq2_size_factors.csv'))
