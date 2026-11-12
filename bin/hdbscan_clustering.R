#!/usr/bin/env Rscript

# hdbscan_clustering.R - HDBSCAN hard clustering on VST counts
#
# Input:  VST counts + filtered gene list
# Output: Cluster assignments, plots

suppressPackageStartupMessages({
  library(optparse)
  library(dbscan)
  library(ggplot2)
  library(dplyr)
  library(uwot)
  library(RSQLite)
})

.find_lib <- function() {
  for (p in c(Sys.which("gene_db_lib.R"), file.path("bin", "gene_db_lib.R")))
    if (nchar(p) > 0 && file.exists(p)) return(p)
  NULL
}
lib <- .find_lib()
if (!is.null(lib)) source(lib)

opts <- list(
  make_option("--vst",       type = "character", help = "VST counts CSV"),
  make_option("--genes",     type = "character", help = "Directory with DE gene lists or a single gene list file"),
  make_option("--metadata",  type = "character", default = NULL, help = "Sample metadata TSV"),
  make_option("--min_pts",         type = "integer",   default = 5, help = "HDBSCAN minPts parameter"),
  make_option("--condition_from",  type = "character", default = "genotype_timepoint", help = "Underscore-joined metadata columns for condition (e.g. genotype_stage)"),
  make_option("--condition_order", type = "character", default = NULL, help = "Comma-separated condition order for plots"),
  make_option("--seed",           type = "integer",   default = 42, help = "Random seed"),
  make_option("--gene_db",        type = "character",  default = NULL, help = "Path to gene_db.sqlite"),
  make_option("--outdir",         type = "character",  default = ".", help = "Output directory")
)

args <- parse_args(OptionParser(option_list = opts))

gdb <- if (exists("maybe_open_gene_db")) maybe_open_gene_db(args$gene_db) else NULL

set.seed(args$seed)
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

# Read data
vst <- read.csv(args$vst, row.names = 1, check.names = FALSE)

# Read gene list
if (file.info(args$genes)$isdir) {
  gene_files <- list.files(args$genes, pattern = "genes_de\\.txt$|all_de_genes\\.txt$",
                           full.names = TRUE, recursive = TRUE)
  genes <- unique(unlist(lapply(gene_files, readLines)))
} else {
  genes <- readLines(args$genes)
}
genes <- intersect(genes, rownames(vst))
message("Clustering ", length(genes), " genes")

# Subset and scale
vst_sub <- t(scale(t(vst[genes, , drop = FALSE])))  # z-score per gene

# Average by condition and order
if (!is.null(args$metadata) && file.exists(args$metadata)) {
  meta <- read.delim(args$metadata, stringsAsFactors = FALSE)
  cond_cols <- strsplit(args$condition_from, "_")[[1]]
  meta$condition <- do.call(paste, c(meta[, cond_cols, drop = FALSE], sep = "_"))

  # Use explicit order if provided, otherwise sort naturally
  if (!is.null(args$condition_order) && nchar(args$condition_order) > 0) {
    conditions <- strsplit(args$condition_order, ",")[[1]]
  } else {
    conditions <- sort(unique(meta$condition))
  }

  vst_avg <- do.call(cbind, lapply(conditions, function(cond) {
    samples <- intersect(meta$sample_id[meta$condition == cond], colnames(vst_sub))
    if (length(samples) == 0) return(NULL)
    rowMeans(vst_sub[, samples, drop = FALSE])
  }))
  colnames(vst_avg) <- conditions
  vst_sub <- vst_avg
}

# Run HDBSCAN
message("Running HDBSCAN (minPts = ", args$min_pts, ")...")
hdb <- hdbscan(vst_sub, minPts = args$min_pts)

# Results
results_df <- data.frame(
  gene_id     = genes,
  cluster     = hdb$cluster,
  membership  = hdb$membership_prob,
  outlier     = hdb$cluster == 0,
  stringsAsFactors = FALSE
)

n_clusters <- length(unique(hdb$cluster[hdb$cluster > 0]))
n_noise    <- sum(hdb$cluster == 0)
message("Found ", n_clusters, " clusters, ", n_noise, " noise points")

if (!is.null(gdb)) results_df <- annotate_df(gdb, results_df)
write.csv(results_df, file.path(args$outdir, "cluster_assignments.csv"), row.names = FALSE)

# Per-cluster gene lists
for (i in sort(unique(hdb$cluster))) {
  label <- if (i == 0) "noise" else paste0("cluster_", i)
  cluster_genes <- results_df$gene_id[results_df$cluster == i]
  writeLines(cluster_genes, file.path(args$outdir, paste0(label, "_genes.txt")))
}

# UMAP visualization of clusters
if (nrow(vst_sub) > 15) {
  n_neighbors <- min(15, nrow(vst_sub) - 1)
  umap_coords <- uwot::umap(vst_sub, n_neighbors = n_neighbors, min_dist = 0.1)

  plot_df <- data.frame(
    UMAP1   = umap_coords[, 1],
    UMAP2   = umap_coords[, 2],
    cluster = factor(hdb$cluster)
  )

  p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
    geom_point(size = 0.5, alpha = 0.6) +
    labs(title = paste0("HDBSCAN Clusters (", n_clusters, " clusters, ", n_noise, " noise)")) +
    theme_minimal(base_size = 12) +
    scale_color_manual(values = c("0" = "grey80", setNames(
      scales::hue_pal()(n_clusters), as.character(1:n_clusters)
    )))

  ggsave(file.path(args$outdir, "plots", "hdbscan_umap.pdf"),
         plot = p, width = 8, height = 6)
}

# Cluster expression profile plots (averaged by condition)
# Use the condition-averaged matrix if available, otherwise raw
plot_mat <- vst_sub
time_labels <- colnames(plot_mat)

for (i in sort(unique(hdb$cluster[hdb$cluster > 0]))) {
  cluster_genes <- results_df$gene_id[results_df$cluster == i]
  n_genes <- length(cluster_genes)

  # Standardize for plotting
  cluster_expr <- t(scale(t(plot_mat[cluster_genes, , drop = FALSE])))
  center <- colMeans(cluster_expr, na.rm = TRUE)

  cluster_long <- data.frame(
    gene = rep(cluster_genes, each = length(time_labels)),
    condition = rep(time_labels, times = n_genes),
    expr = as.vector(t(cluster_expr)),
    stringsAsFactors = FALSE
  )
  cluster_long$condition <- factor(cluster_long$condition, levels = time_labels)

  center_df <- data.frame(condition = factor(time_labels, levels = time_labels), expr = center)

  p <- ggplot() +
    geom_line(data = cluster_long, aes(x = condition, y = expr, group = gene),
              color = "grey60", linewidth = 0.3, alpha = 0.3) +
    geom_line(data = center_df, aes(x = condition, y = expr, group = 1),
              color = "red", linewidth = 1.2) +
    geom_point(data = center_df, aes(x = condition, y = expr), color = "red", size = 2) +
    labs(title = paste0("HDBSCAN Cluster ", i, " (", n_genes, " genes)"),
         x = "Condition", y = "Standardized Expression") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(args$outdir, "plots", paste0("cluster_", i, "_profile.pdf")),
         plot = p, width = 6, height = 4)
}

# Cluster size summary
summary_df <- results_df %>%
  group_by(cluster) %>%
  summarise(n_genes = n(), mean_membership = mean(membership), .groups = "drop")
write.csv(summary_df, file.path(args$outdir, "cluster_summary.csv"), row.names = FALSE)

if (!is.null(gdb)) close_gene_db(gdb)
message("Done. Output: ", args$outdir)
