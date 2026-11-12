#!/usr/bin/env Rscript

# pca.R - PCA on VST-normalized counts
#
# Input:  vst_counts.csv + sample_meta.tsv
# Output: PCA coordinates, plots colored by each metadata column

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
})

opts <- list(
  make_option("--vst",      type = "character", help = "VST counts CSV (genes x samples)"),
  make_option("--vst_raw",  type = "character", default = NULL,
              help = "Optional uncorrected VST CSV - when provided, paired before/after plots are generated"),
  make_option("--metadata", type = "character", help = "Sample metadata TSV"),
  make_option("--color_by", type = "character", default = NULL,
              help = "Comma-separated metadata columns to color by (default: all)"),
  make_option("--top_genes", type = "integer",  default = 500, help = "Top N variable genes for PCA"),
  make_option("--outdir",   type = "character",  default = ".", help = "Output directory")
)

args <- parse_args(OptionParser(option_list = opts))

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# Read metadata
meta <- read.delim(args$metadata, stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id

# Helper: run PCA on a VST matrix, return (pca_df, var_pct, loadings)
run_pca <- function(vst_path, meta, top_n) {
  vst <- read.csv(vst_path, row.names = 1, check.names = FALSE)
  shared <- intersect(colnames(vst), meta$sample_id)
  vst    <- vst[, shared]
  m      <- meta[shared, ]
  gene_vars <- apply(vst, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(top_n, nrow(vst))]
  vst_top   <- vst[top_genes, ]
  pr <- prcomp(t(vst_top), center = TRUE, scale. = TRUE)
  var_pct <- round(100 * summary(pr)$importance[2, ] * 100) / 100
  df <- as.data.frame(pr$x[, 1:min(5, ncol(pr$x))])
  df$sample_id <- rownames(df)
  df <- merge(df, m, by = "sample_id")
  loadings <- as.data.frame(pr$rotation[, 1:2])
  loadings$gene_id <- rownames(loadings)
  loadings <- loadings %>% arrange(desc(abs(PC1)))
  list(df = df, var_pct = var_pct, loadings = loadings)
}

# Always run PCA on the active VST
active <- run_pca(args$vst, meta, args$top_genes)

# Optionally also run on raw VST for before/after comparison
has_raw <- !is.null(args$vst_raw) && file.exists(args$vst_raw) &&
  normalizePath(args$vst_raw, mustWork = FALSE) != normalizePath(args$vst, mustWork = FALSE)
raw <- if (has_raw) run_pca(args$vst_raw, meta, args$top_genes) else NULL

write.csv(active$df, file.path(args$outdir, "pca_coordinates.csv"), row.names = FALSE)
if (has_raw) {
  write.csv(raw$df, file.path(args$outdir, "pca_coordinates_raw.csv"), row.names = FALSE)
}

# Determine columns to color by
if (!is.null(args$color_by)) {
  color_cols <- strsplit(args$color_by, ",")[[1]]
} else {
  exclude <- c("sample_id", "bam", "fastq_1", "fastq_2", "strandedness", "input_type")
  color_cols <- setdiff(colnames(meta), exclude)
}

make_plot <- function(df, var_pct, col, title) {
  ggplot(df, aes(x = PC1, y = PC2, color = .data[[col]])) +
    geom_point(size = 3) +
    geom_text(aes(label = sample_id), size = 2, vjust = -1, show.legend = FALSE) +
    labs(
      x     = paste0("PC1 (", var_pct[1], "%)"),
      y     = paste0("PC2 (", var_pct[2], "%)"),
      title = title,
      color = col
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
}

for (col in color_cols) {
  if (!col %in% colnames(active$df)) next

  if (has_raw) {
    p_raw    <- make_plot(raw$df,    raw$var_pct,    col, paste0("PCA (uncorrected VST) - ", col))
    p_active <- make_plot(active$df, active$var_pct, col, paste0("PCA (corrected VST) - ", col))
    ggsave(file.path(args$outdir, paste0("pca_", col, "_raw.pdf")),
           plot = p_raw, width = 8, height = 6)
    ggsave(file.path(args$outdir, paste0("pca_", col, "_compare.pdf")),
           plot = patchwork::wrap_plots(p_raw, p_active, ncol = 2),
           width = 14, height = 6)
  }
  ggsave(file.path(args$outdir, paste0("pca_", col, ".pdf")),
         plot = make_plot(active$df, active$var_pct, col, paste("PCA colored by", col)),
         width = 8, height = 6)
}

writeLines(
  paste0("PC", seq_along(active$var_pct), ": ", active$var_pct, "%"),
  file.path(args$outdir, "variance_explained.txt")
)
write.csv(active$loadings, file.path(args$outdir, "loadings_pc1_pc2.csv"), row.names = FALSE)

message("Done. Generated PCA plots", if (has_raw) " (with before/after pairs)" else "",
        " in ", args$outdir)
