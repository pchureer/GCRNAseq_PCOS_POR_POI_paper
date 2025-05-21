#!/usr/bin/env Rscript
#===============================================================================
# scripts/04_heatmap/04_plot_heatmap.R
#
# PURPOSE:
#   Visualize log₂ fold-change values of the 17-gene core signature across
#   all six dual-control contrasts as a clustered heatmap.
#
# USAGE:
#   Rscript scripts/04_heatmap/04_plot_heatmap.R \
#     results/log2FC_17core.csv \
#     figures/core17_heatmap.png
#
# ARGS:
#   1) Path to CSV of log2FC values (rows=genes, cols=contrasts)
#   2) Output PNG filename
#===============================================================================

suppressPackageStartupMessages({
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat("Usage: 04_plot_heatmap.R <log2FC_matrix.csv> <output.png>\n")
  quit(status = 1)
}

log2fc_file <- args[1]
out_png     <- args[2]

# Read in the log2FC matrix
mat <- read.csv(log2fc_file, row.names = 1, check.names = FALSE)
mat <- as.matrix(mat)

# Optional: scale rows (genes) if desired
# mat <- t(scale(t(mat)))

# Define a color palette
palette <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# Create the heatmap and save as PNG
png(filename = out_png, width = 800, height = 600, res = 150)
pheatmap(
  mat,
  color           = palette,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  show_rownames   = TRUE,
  show_colnames   = TRUE,
  fontsize_row    = 8,
  fontsize_col    = 8,
  border_color    = NA,
  main            = "Core 17-Gene Signature log2FC Heatmap"
)
dev.off()
