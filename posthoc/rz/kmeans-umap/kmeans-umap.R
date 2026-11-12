#!/usr/bin/env Rscript
# Re-render the k-means / HDBSCAN UMAP panels at paper type size.
#
# The pipeline (bin/kmeans_clustering.R) writes only the plot, not the
# embedding, so the coordinates have to be recomputed. That is reproducible
# because the pipeline seeds immediately before uwot::umap, and because
# cluster_assignments.csv preserves the exact row order of the matrix that was
# embedded (assignments_df$gene_id = rownames(vst_sub), and the annotation step
# restores the original order). Reproducing the pipeline's preprocessing
# byte-for-byte is therefore the whole job:
#
#   vst -> subset to the assignment gene order -> average by condition
#       -> z-score rows -> set.seed(seed) -> uwot::umap(n_neighbors, min_dist)
#
# Any drift in condition_from / condition_order changes the matrix and hence
# the embedding, so both are required arguments, taken from the run's YAML.

suppressPackageStartupMessages({
  library(optparse); library(ggplot2); library(uwot)
})
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "typography.R"))

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "palette.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--vst",      type = "character"),
  make_option("--meta",     type = "character"),
  make_option("--clusters", type = "character", help = "cluster_assignments.csv from the run"),
  make_option("--condition-from",  type = "character", default = "genotype_stage", dest = "condition_from"),
  make_option("--condition-order", type = "character", default = "", dest = "condition_order"),
  make_option("--seed",     type = "integer", default = 42, help = "Must match the pipeline's --seed"),
  make_option("--title",    type = "character", default = ""),
  make_option("--mode",     type = "character", default = "kmeans",
              help = paste("kmeans or hdbscan. hdbscan uses that script's point",
                           "style and greys cluster 0 (noise). NOTE: the pipeline",
                           "does NOT re-seed immediately before uwot::umap in the",
                           "hdbscan path, so an exact match is not guaranteed --",
                           "verify the result against the published panel.")),
  make_option("--font-scale", type = "double", default = 1, dest = "font_scale"),
  make_option("--point-size", type = "double", default = 1.2, dest = "point_size",
              help = "kmeans point size; the pipeline's 1.2 reads crowded on 5k genes"),
  make_option("--legend",   type = "character", default = "right",
              help = "right, or bottom (two rows of keys, for one-third width)"),
  make_option("--width",    type = "double", default = 8),
  make_option("--height",   type = "double", default = 6),
  make_option("--outdir",   type = "character")
)))
stopifnot(!is.null(opt$vst), !is.null(opt$meta), !is.null(opt$clusters), !is.null(opt$outdir))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# --- verbatim from bin/kmeans_clustering.R -----------------------------------
average_by_condition <- function(mat, metadata_path, condition_from, condition_order) {
  meta <- read.delim(metadata_path, stringsAsFactors = FALSE)
  cond_cols <- strsplit(condition_from, "_", fixed = TRUE)[[1]]
  meta$condition <- do.call(paste, c(meta[, cond_cols, drop = FALSE], sep = "_"))
  conditions <- if (nzchar(condition_order)) {
    strsplit(condition_order, ",", fixed = TRUE)[[1]]
  } else sort(unique(meta$condition))
  averaged <- lapply(conditions, function(cond) {
    samples <- intersect(meta$sample_id[meta$condition == cond], colnames(mat))
    if (length(samples) == 0) stop("no samples for condition: ", cond)
    rowMeans(mat[, samples, drop = FALSE])
  })
  avg <- do.call(cbind, averaged); colnames(avg) <- conditions; avg
}
zscore_rows <- function(mat) {
  scaled <- t(scale(t(mat)))
  scaled[complete.cases(scaled), , drop = FALSE]
}
# -----------------------------------------------------------------------------

vst <- read.csv(opt$vst, row.names = 1, check.names = FALSE)
km  <- read.csv(opt$clusters, stringsAsFactors = FALSE)

missing <- setdiff(km$gene_id, rownames(vst))
if (length(missing) > 0) stop(length(missing), " assignment genes absent from the VST matrix")

# row order must follow cluster_assignments.csv exactly -- UMAP is order-sensitive
vst_sub <- vst[km$gene_id, , drop = FALSE]
vst_sub <- average_by_condition(vst_sub, opt$meta, opt$condition_from, opt$condition_order)
vst_sub <- zscore_rows(as.matrix(vst_sub))

if (nrow(vst_sub) != nrow(km)) {
  stop("row count changed during preprocessing (", nrow(km), " -> ", nrow(vst_sub),
       "); the embedding would not match the published panel")
}

set.seed(opt$seed)
emb <- uwot::umap(vst_sub, n_neighbors = min(15, nrow(vst_sub) - 1),
                  min_dist = 0.1, ret_model = FALSE, verbose = FALSE)

plot_df <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2],
                      cluster = factor(km$cluster))
ttl <- if (nzchar(opt$title)) opt$title else
  paste0("K-means clusters (k = ", nlevels(plot_df$cluster), ")")

p <- if (opt$mode == "hdbscan") {
  n_clusters <- sum(levels(plot_df$cluster) != "0")
  n_noise    <- sum(plot_df$cluster == "0")
  if (!nzchar(opt$title))
    ttl <- paste0("HDBSCAN (", n_clusters, " clusters, ", n_noise, " noise)")
  ggplot(plot_df, aes(UMAP1, UMAP2, colour = cluster)) +
    geom_point(size = 0.5, alpha = 0.6) +
    # strong hues to the big clusters: colours are assigned by cluster size,
    # so the two large HDBSCAN clusters never land on yellow or black
    # HDBSCAN labels unclustered genes 0; key them as "noise", last, so 0 is not read as a cluster
    scale_color_manual(values = c("0" = "grey80",
      setNames(PAL$okabe[rank(-table(plot_df$cluster)[as.character(seq_len(n_clusters))], ties.method = "first")],
               as.character(seq_len(n_clusters)))),
      breaks = c(as.character(seq_len(n_clusters)), "0"),
      labels = c(as.character(seq_len(n_clusters)), "noise")) +
    labs(title = ttl, x = "UMAP1", y = "UMAP2") +
    theme_minimal(base_size = 12 * opt$font_scale)
} else {
  ggplot(plot_df, aes(UMAP1, UMAP2, colour = cluster)) +
    geom_point(size = opt$point_size, alpha = 0.8) +
    scale_color_manual(values = setNames(PAL$okabe[seq_len(nlevels(plot_df$cluster))], levels(plot_df$cluster))) +
    labs(title = ttl, x = "UMAP1", y = "UMAP2") +
    theme_minimal(base_size = 12 * opt$font_scale)
}
p <- p +
  # ticks on multiples of 5, a little air between panel and legend
  scale_x_continuous(breaks = scales::breaks_width(5)) +
  scale_y_continuous(breaks = scales::breaks_width(5)) +
  # a little air between the cluster keys (1, 2, 3 ...), not between plot and legend
  theme(legend.key.spacing.y = unit(5 * opt$font_scale, "pt"),
        legend.title = element_text(size = 9.6 * opt$font_scale),   # same 10pt as its keys
        legend.position = opt$legend)
if (opt$legend == "bottom")
  p <- p + guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2.5, alpha = 1)))

ggsave(file.path(opt$outdir, paste0(opt$mode, "_umap.pdf")), p, width = opt$width, height = opt$height)
ggsave(file.path(opt$outdir, paste0(opt$mode, "_umap.png")), p, width = opt$width, height = opt$height, dpi = 200)
message("wrote ", opt$outdir, "/", opt$mode, "_umap.{pdf,png}  (", nrow(plot_df), " genes, ",
        nlevels(plot_df$cluster), " clusters)")
