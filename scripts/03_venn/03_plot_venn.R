#!/usr/bin/env Rscript
#===============================================================================
# scripts/03_venn/03_plot_venn.R
#
# PURPOSE:
#   Plot a 3‐set Venn diagram showing overlap of significantly differentially
#   expressed genes (FDR < 0.05) across PCOS, POR, and POI pooled‐control analyses.
#
# USAGE:
#   Rscript scripts/03_venn/03_plot_venn.R \
#     results/DE_results_PCOS_pooled.csv \
#     results/DE_results_POR_pooled.csv \
#     results/DE_results_POI_pooled.csv \
#     figures/venn_pooled.png
#
# ARGS:
#   1) CSV of DE results for PCOS vs pooled control
#   2) CSV of DE results for POR vs pooled control
#   3) CSV of DE results for POI vs pooled control
#   4) Output path for the Venn diagram PNG
#===============================================================================

suppressPackageStartupMessages({
  library(ggVennDiagram)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  cat("Usage: 03_plot_venn.R <PCOS.csv> <POR.csv> <POI.csv> <output.png>\n")
  quit(status = 1)
}

# Read in DE result tables
de_pcos <- read.csv(args[1], row.names = 1, check.names = FALSE)
de_por  <- read.csv(args[2], row.names = 1, check.names = FALSE)
de_poi  <- read.csv(args[3], row.names = 1, check.names = FALSE)

# Subset to significant genes (FDR < 0.05)
genes_pcos <- rownames(subset(de_pcos, FDR < 0.05))
genes_por  <- rownames(subset(de_por,  FDR < 0.05))
genes_poi  <- rownames(subset(de_poi,  FDR < 0.05))

# Create named list for Venn
gene_sets <- list(
  PCOS = genes_pcos,
  POR  = genes_por,
  POI  = genes_poi
)

# Generate Venn diagram
p <- ggVennDiagram(gene_sets, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Shared DE Genes (FDR < 0.05) Across PCOS, POR, and POI")

# Save to file
ggsave(filename = args[4], plot = p, width = 6, height = 6, dpi = 300)
