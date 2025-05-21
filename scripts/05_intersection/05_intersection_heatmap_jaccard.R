#!/usr/bin/env Rscript
#===============================================================================
# scripts/05_intersection/05_intersection_heatmap_jaccard.R
#
# PURPOSE:
#   Compute intersection counts and Jaccard similarity between six DE gene sets
#   (from pooled- and matched-control analyses for PCOS, POR, POI) and three
#   Open Targets reference lists, then plot two heatmaps.
#
# USAGE:
#   Rscript scripts/05_intersection/05_intersection_heatmap_jaccard.R \
#     results/DE_results_PCOS_pooled.csv \
#     results/DE_results_PCOS_matched.csv \
#     results/DE_results_POR_pooled.csv \
#     results/DE_results_POR_matched.csv \
#     results/DE_results_POI_pooled.csv \
#     results/DE_results_POI_matched.csv \
#     data/OpenTargets_PCOS.csv \
#     data/OpenTargets_OvarianDysfunction.csv \
#     data/OpenTargets_POI.csv \
#     figures/intersect_count_heatmap.png \
#     figures/jaccard_heatmap.png
#
# ARGS:
#   1-6) DE result CSVs (FDR < 0.05 gene sets) in this order:
#         PCOS_pooled, PCOS_matched,
#         POR_pooled,  POR_matched,
#         POI_pooled,  POI_matched
#   7-9) Open Targets reference gene lists (one gene per line, no header)
#         OpenTargets_PCOS, OpenTargets_OvarianDysfunction, OpenTargets_POI
#   10) Output PNG for intersection count heatmap
#   11) Output PNG for Jaccard similarity heatmap
#===============================================================================

suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11) {
  cat("Usage: 05_intersection_heatmap_jaccard.R DE1 DE2 DE3 DE4 DE5 DE6 OT1 OT2 OT3 intersect.png jaccard.png\n")
  quit(status = 1)
}

# Assign arguments
de_files     <- args[1:6]
ot_files     <- args[7:9]
out_count_png <- args[10]
out_jacc_png  <- args[11]

# Read DE gene sets (FDR < 0.05)
de_sets <- lapply(de_files, function(f) {
  df <- read.csv(f, row.names = 1, check.names = FALSE)
  rownames(subset(df, FDR < 0.05))
})
names(de_sets) <- c(
  "PCOS_pooled","PCOS_matched",
  "POR_pooled", "POR_matched",
  "POI_pooled", "POI_matched"
)

# Read Open Targets gene lists (one column, no header)
ot_sets <- lapply(ot_files, function(f) {
  unique(scan(f, what = character(), sep = "\n", quiet = TRUE))
})
names(ot_sets) <- c("OT_PCOS","OT_Ovary","OT_POI")

# Combine names and initialize matrices
all_names <- c(names(de_sets), names(ot_sets))
n <- length(all_names)
mat_count <- matrix(0, nrow = n, ncol = n,
                    dimnames = list(all_names, all_names))
mat_jacc  <- mat_count

# Compute intersection counts and Jaccard indices
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    A <- if (i <= 6) de_sets[[i]] else ot_sets[[i - 6]]
    B <- if (j <= 6) de_sets[[j]] else ot_sets[[j - 6]]
    inter <- length(intersect(A, B))
    uni   <- length(union(A, B))
    mat_count[i, j] <- inter
    mat_jacc[i, j]  <- if (uni > 0) inter / uni else NA
  }
}

# Melt for plotting
df_count <- melt(mat_count, varnames = c("Set1", "Set2"), value.name = "Count")
df_jacc  <- melt(mat_jacc,  varnames = c("Set1", "Set2"), value.name = "Jaccard")

# Plot intersection count heatmap
p1 <- ggplot(df_count, aes(x = Set1, y = Set2, fill = Count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Count), size = 3) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = NULL, y = NULL, fill = "Count",
       title = "DE vs. Reference: Intersection Counts") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(out_count_png, p1, width = 8, height = 6, dpi = 300)

# Plot Jaccard similarity heatmap
p2 <- ggplot(df_jacc, aes(x = Set1, y = Set2, fill = Jaccard)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Jaccard)), size = 3) +
  scale_fill_gradient(low = "white", high = "blue", na.value = "grey90") +
  labs(x = NULL, y = NULL, fill = "Jaccard",
       title = "DE vs. Reference: Jaccard Similarity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(out_jacc_png, p2, width = 8, height = 6, dpi = 300)
