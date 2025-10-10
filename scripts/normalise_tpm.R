#!/usr/bin/env Rscript

# Script to normalise counts usinge TPM method

suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicFeatures)
})

### Command arguments ###
args <- commandArgs(trailingOnly = TRUE)
countfile <- args[1]
gtffile <- args[2]
outfile <- args[3]

### Test Zone ###
# countfile <- "~/Library/CloudStorage/OneDrive-Aarhusuniversitet/PhD project SLA/Projects/rna_seq_snake_pipeline/temp/count/counts_raw.csv"
# gtffile <- "~/Library/CloudStorage/OneDrive-Aarhusuniversitet/PhD project SLA/Projects/rna_seq_snake_pipeline/temp/genome/Homo_sapiens.GRCh38.115.gtf"
# outfile <- "~/Library/CloudStorage/OneDrive-Aarhusuniversitet/PhD project SLA/Projects/rna_seq_snake_pipeline/temp/count/counts_normalised_tpm.csv"


# Load counts
counts <- read.csv(file = countfile,
                   header = TRUE)

# Remove non-genes counts
counts_fltr <- counts[1:(nrow(counts)-5),]
row_order <- counts_fltr$Gene.ID
rownames(counts_fltr) <- counts_fltr$Gene.ID

gene_ids <- counts_fltr[[2]]
count_values <- counts_fltr[,5:ncol(counts)]
rownames(count_values) <- gene_ids

# Extract gene lengths from GTF
txdb <- makeTxDbFromGFF(gtffile, format = "gtf")
gene_lengths <- exonsBy(txdb, by = "gene") %>%
  lapply(function(x) sum(width(reduce(x)))) %>%
  unlist()

# Keep only matching genes
gene_lengths <- gene_lengths[names(gene_lengths) %in% rownames(count_values)]
count_values <- count_values[names(gene_lengths), ]

# Calculate TPM
length_kb <- gene_lengths / 1000
rpk <- count_values[names(gene_lengths),]
for (sample in colnames(count_values)) {
  rpk[,sample] <- rpk[,sample]/length_kb
}
tpm <- rpk
for (sample in colnames(count_values)) {
  tpm[,sample] <- (tpm[,sample] * 1000000)/colSums(rpk)[sample]
}
tpm[,'Gene ID'] <- names(gene_lengths)

# Add gene names and biotypes
tpm[gene_ids,'Gene name'] <- counts_fltr[gene_ids, 'Gene.name']
tpm[gene_ids,'Gene biotype'] <- counts_fltr[gene_ids, 'Gene.biotype']


# ---- Write output ----
write.csv(tpm[, c('Gene ID', 'Gene name', 'Gene biotype', colnames(count_values))],
          file = paste0(outfile))

