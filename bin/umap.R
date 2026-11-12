#!/usr/bin/env Rscript

# umap.R - UMAP on VST-normalized counts
#
# Input:  vst_counts.csv + sample_meta.tsv
# Output: UMAP coordinates, plots colored by each metadata column

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(uwot)
  library(dplyr)
})

opts <- list(
  make_option("--vst",         type = "character", help = "VST counts CSV (genes x samples)"),
  make_option("--vst_raw",     type = "character", default = NULL,
              help = "Optional uncorrected VST CSV - when provided, paired before/after plots are generated"),
  make_option("--metadata",    type = "character", help = "Sample metadata TSV"),
  make_option("--color_by",    type = "character", default = NULL,
              help = "Comma-separated metadata columns to color by (default: all)"),
  make_option("--top_genes",   type = "integer",   default = 500, help = "Top N variable genes"),
  make_option("--n_neighbors", type = "integer",   default = 15,  help = "UMAP n_neighbors"),
  make_option("--min_dist",    type = "numeric",   default = 0.1, help = "UMAP min_dist"),
  make_option("--seed",        type = "integer",   default = 42,  help = "Random seed"),
  make_option("--outdir",      type = "character",  default = ".",  help = "Output directory")
)

args <- parse_args(OptionParser(option_list = opts))

set.seed(args$seed)
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

meta <- read.delim(args$metadata, stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id

run_umap <- function(vst_path, meta, top_n, n_neighbors_req, min_dist) {
  vst <- read.csv(vst_path, row.names = 1, check.names = FALSE)
  shared <- intersect(colnames(vst), meta$sample_id)
  vst    <- vst[, shared]
  m      <- meta[shared, ]
  gene_vars <- apply(vst, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(top_n, nrow(vst))]
  vst_top   <- vst[top_genes, ]
  n_neighbors <- min(n_neighbors_req, ncol(vst_top) - 1)
  u <- uwot::umap(t(vst_top), n_neighbors = n_neighbors, min_dist = min_dist,
                  metric = "euclidean", n_components = 2)
  df <- data.frame(sample_id = colnames(vst_top), UMAP1 = u[, 1], UMAP2 = u[, 2])
  merge(df, m, by = "sample_id")
}

active <- run_umap(args$vst, meta, args$top_genes, args$n_neighbors, args$min_dist)
has_raw <- !is.null(args$vst_raw) && file.exists(args$vst_raw) &&
  normalizePath(args$vst_raw, mustWork = FALSE) != normalizePath(args$vst, mustWork = FALSE)
raw <- if (has_raw) run_umap(args$vst_raw, meta, args$top_genes, args$n_neighbors, args$min_dist) else NULL

write.csv(active, file.path(args$outdir, "umap_coordinates.csv"), row.names = FALSE)
if (has_raw) {
  write.csv(raw, file.path(args$outdir, "umap_coordinates_raw.csv"), row.names = FALSE)
}

if (!is.null(args$color_by)) {
  color_cols <- strsplit(args$color_by, ",")[[1]]
} else {
  exclude <- c("sample_id", "bam", "fastq_1", "fastq_2", "strandedness", "input_type")
  color_cols <- setdiff(colnames(meta), exclude)
}

make_plot <- function(df, col, title) {
  ggplot(df, aes(x = UMAP1, y = UMAP2, color = .data[[col]])) +
    geom_point(size = 3) +
    geom_text(aes(label = sample_id), size = 2, vjust = -1, show.legend = FALSE) +
    labs(title = title, color = col) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
}

for (col in color_cols) {
  if (!col %in% colnames(active)) next

  if (has_raw) {
    p_raw    <- make_plot(raw,    col, paste0("UMAP (uncorrected VST) - ", col))
    p_active <- make_plot(active, col, paste0("UMAP (corrected VST) - ", col))
    ggsave(file.path(args$outdir, paste0("umap_", col, "_raw.pdf")),
           plot = p_raw, width = 8, height = 6)
    ggsave(file.path(args$outdir, paste0("umap_", col, "_compare.pdf")),
           plot = patchwork::wrap_plots(p_raw, p_active, ncol = 2),
           width = 14, height = 6)
  }
  ggsave(file.path(args$outdir, paste0("umap_", col, ".pdf")),
         plot = make_plot(active, col, paste("UMAP colored by", col)),
         width = 8, height = 6)
}

message("Done. Generated UMAP plots", if (has_raw) " (with before/after pairs)" else "",
        " in ", args$outdir)
