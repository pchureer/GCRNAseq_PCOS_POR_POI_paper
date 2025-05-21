#!/usr/bin/env Rscript
#===============================================================================
# scripts/06_enrichment/06_text_enrichment.R
#
# PURPOSE:
#   Generate a word-cloud of significantly enriched GO Biological Process terms.
#
# USAGE:
#   Rscript scripts/06_enrichment/06_text_enrichment.R \
#     <GO_enrichment.csv> \
#     <output_wordcloud.png>
#
# ARGS:
#   1) GO_enrichment.csv    CSV with columns: Term, Count, FDR
#   2) output_wordcloud.png Path to save the word-cloud PNG
#
# OUTPUT:
#   A word-cloud image where term size = gene count and color = –log10(FDR)
#===============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggwordcloud)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat("Usage: 06_text_enrichment.R <GO_enrichment.csv> <output_wordcloud.png>\n")
  quit(status = 1)
}

go_file  <- args[1]
out_png  <- args[2]

# Read GO enrichment results
go_df <- read.csv(go_file, stringsAsFactors = FALSE)

# Keep only significant terms
go_df <- subset(go_df, FDR < 0.05)

# Create labels and aesthetics
go_df$label <- paste0(go_df$Term, " (", go_df$Count, ")")
go_df$score <- -log10(go_df$FDR)

# Plot word-cloud
p <- ggplot(go_df, aes(label = label, size = Count, color = score)) +
  geom_text_wordcloud_area() +
  scale_size_area(max_size = 15) +
  scale_color_viridis_c(option = "plasma", name = "-log10(FDR)") +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("GO Biological Process Enrichment Word Cloud")

# Save to file
ggsave(filename = out_png, plot = p, width = 8, height = 6, dpi = 300)
