#!/usr/bin/env Rscript
#===============================================================================
# scripts/02_dge/02_edgeR_DE.R
#
# PURPOSE:
#   Normalize raw gene counts, fit a GLM with batch and group,
#   and perform differential expression testing using edgeR.
#
# USAGE:
#   Rscript 02_edgeR_DE.R \
#     <counts_matrix.csv> \
#     <metadata.csv> \
#     <output_DE_results.csv>
#
# ARGS:
#   counts_matrix.csv  CSV of raw counts (genes × samples; colnames = sample IDs)
#   metadata.csv       CSV of sample metadata with columns: sampleID, group, batch
#   output_DE_results.csv
#
# OUTPUT:
#   A CSV with columns: gene, log2FC, PValue, FDR, etc., for each gene
#===============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  cat("Usage: Rscript 02_edgeR_DE.R counts_matrix.csv metadata.csv output_DE_results.csv\n")
  quit(status = 1)
}
counts_file <- args[1]
meta_file   <- args[2]
out_file    <- args[3]

# Load counts and metadata
cat("Loading counts from", counts_file, "...\n")
counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
cat("Loading metadata from", meta_file, "...\n")
meta   <- read.csv(meta_file, row.names = 1, check.names = FALSE)

# Ensure samples match between counts and metadata
if (!all(colnames(counts) %in% rownames(meta))) {
  stop("Some samples in counts_matrix.csv are missing from metadata.csv")
}
meta <- meta[colnames(counts), ]

# Create DGEList
dge <- DGEList(counts = counts, group = meta$group)

# Filter out lowly expressed genes
cat("Filtering lowly expressed genes...\n")
design0 <- model.matrix(~ batch + group, data = meta)
keep <- filterByExpr(dge, design = design0)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Calculate normalization factors
cat("Calculating TMM normalization factors...\n")
dge <- calcNormFactors(dge)

# Create design matrix for GLM (no intercept)
cat("Creating design matrix...\n")
design <- model.matrix(~ batch + group, data = meta)

# Estimate dispersions
cat("Estimating dispersions...\n")
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

# Fit GLM and perform likelihood ratio test
cat("Fitting GLM and testing for differential expression...\n")
fit <- glmFit(dge, design)
lrt <- glmLRT(fit)

# Extract results
cat("Extracting results and writing to", out_file, "...\n")
res <- topTags(lrt, n = nrow(dge$counts), sort.by = "none")$table
res$gene <- rownames(res)
# Reorder columns: gene first
res <- res[, c("gene", setdiff(colnames(res), "gene"))]

write.csv(res, file = out_file, row.names = FALSE)
cat("Done.\n")
