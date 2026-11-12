#!/usr/bin/env Rscript

# hierarchical_clustering.R - Agglomerative hierarchical clustering on VST counts
#
# Input:  VST counts + filtered gene list
# Output: Cluster assignments, dendrogram, heatmap, per-cluster profile plots

suppressPackageStartupMessages({
  library(optparse)
  library(cluster)
  library(ggplot2)
  library(dplyr)
  library(uwot)
  library(RSQLite)
  library(ComplexHeatmap)
  library(circlize)
  library(dynamicTreeCut)
})

.find_lib <- function() {
  for (p in c(Sys.which("gene_db_lib.R"), file.path("bin", "gene_db_lib.R"))) {
    if (nchar(p) > 0 && file.exists(p)) return(p)
  }
  NULL
}

lib <- .find_lib()
if (!is.null(lib)) source(lib)

opts <- list(
  make_option("--vst", type = "character", help = "VST counts CSV"),
  make_option("--genes", type = "character", help = "Directory with DE gene lists or a single gene list file"),
  make_option("--metadata", type = "character", default = NULL, help = "Sample metadata TSV"),
  make_option("--num_clusters", type = "character", default = "auto", help = "Number of clusters or 'auto'"),
  make_option("--cluster_range", type = "character", default = "2,12", help = "Min,max for auto cluster selection"),
  make_option("--condition_from", type = "character", default = "genotype_timepoint", help = "Underscore-joined metadata columns for condition"),
  make_option("--condition_order", type = "character", default = NULL, help = "Comma-separated condition order for plots"),
  make_option("--linkage", type = "character", default = "ward.D2", help = "hclust linkage method"),
  make_option("--distance", type = "character", default = "euclidean", help = "Distance: euclidean | correlation"),
  make_option("--auto_method", type = "character", default = "dynamic_tree_cut",
              help = "Auto-k method: dynamic_tree_cut | silhouette"),
  make_option("--deep_split", type = "integer", default = 2,
              help = "dynamicTreeCut deepSplit (0-4); higher = more, smaller clusters"),
  make_option("--min_cluster_size", type = "integer", default = 30,
              help = "dynamicTreeCut minClusterSize"),
  make_option("--max_plot_genes", type = "integer", default = 2000, help = "Max gene rows in heatmap (downsamples by variance for the plot only)"),
  make_option("--seed", type = "integer", default = 42, help = "Random seed"),
  make_option("--gene_db", type = "character", default = NULL, help = "Path to gene_db.sqlite"),
  make_option("--outdir", type = "character", default = ".", help = "Output directory")
)

args <- parse_args(OptionParser(option_list = opts))

gdb <- if (exists("maybe_open_gene_db")) maybe_open_gene_db(args$gene_db) else NULL
if (!is.null(gdb) && exists("close_gene_db")) {
  on.exit(close_gene_db(gdb), add = TRUE)
}

set.seed(args$seed)
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

read_gene_list <- function(path) {
  if (isTRUE(file.info(path)$isdir)) {
    gene_files <- list.files(
      path,
      pattern = "genes_de\\.txt$|all_de_genes\\.txt$",
      full.names = TRUE,
      recursive = TRUE
    )
    unique(unlist(lapply(gene_files, readLines)))
  } else {
    readLines(path)
  }
}

average_by_condition <- function(mat, metadata_path, condition_from, condition_order) {
  if (is.null(metadata_path)) return(mat)
  if (!file.exists(metadata_path)) stop("Metadata file not found: ", metadata_path)

  meta <- read.delim(metadata_path, stringsAsFactors = FALSE)
  cond_cols <- strsplit(condition_from, "_", fixed = TRUE)[[1]]

  if (!"sample_id" %in% colnames(meta)) stop("Metadata must contain a sample_id column")
  missing_cols <- setdiff(cond_cols, colnames(meta))
  if (length(missing_cols) > 0) stop("Metadata missing condition columns: ", paste(missing_cols, collapse = ", "))

  meta$condition <- do.call(paste, c(meta[, cond_cols, drop = FALSE], sep = "_"))

  if (!is.null(condition_order) && nchar(condition_order) > 0) {
    requested <- strsplit(condition_order, ",", fixed = TRUE)[[1]]
    missing_conditions <- setdiff(requested, unique(meta$condition))
    if (length(missing_conditions) > 0) {
      stop("Requested conditions in --condition_order not found in metadata: ",
           paste(missing_conditions, collapse = ", "))
    }
    conditions <- requested
  } else {
    conditions <- sort(unique(meta$condition))
  }

  samples_by_condition <- lapply(conditions, function(cond) {
    intersect(meta$sample_id[meta$condition == cond], colnames(mat))
  })
  empty <- conditions[vapply(samples_by_condition, length, integer(1)) == 0]
  if (length(empty) > 0) {
    stop("Metadata conditions had no matching samples: ", paste(empty, collapse = ", "))
  }

  averaged <- lapply(samples_by_condition, function(samples) {
    rowMeans(mat[, samples, drop = FALSE])
  })
  avg_mat <- do.call(cbind, averaged)
  colnames(avg_mat) <- conditions
  avg_mat
}

zscore_rows <- function(mat) {
  scaled <- t(scale(t(mat)))
  keep <- complete.cases(scaled)
  dropped <- sum(!keep)
  if (dropped > 0) message("Dropping ", dropped, " genes with incomplete standardized profiles")
  scaled[keep, , drop = FALSE]
}

parse_cluster_range <- function(range_text) {
  parts <- as.integer(strsplit(range_text, ",", fixed = TRUE)[[1]])
  if (length(parts) != 2 || any(is.na(parts))) stop("--cluster_range must be min,max")
  if (parts[1] < 2 || parts[2] < parts[1]) stop("--cluster_range must satisfy min >= 2 and max >= min")
  parts
}

compute_distance <- function(mat, method) {
  if (method == "correlation") {
    cor_mat <- stats::cor(t(mat), use = "pairwise.complete.obs")
    stats::as.dist(1 - cor_mat)
  } else if (method == "euclidean") {
    stats::dist(mat, method = "euclidean")
  } else {
    stop("Unknown --distance: ", method)
  }
}

vst <- read.csv(args$vst, row.names = 1, check.names = FALSE)
genes <- intersect(read_gene_list(args$genes), rownames(vst))
if (length(genes) < 2) stop("Need at least 2 intersecting genes for hierarchical clustering")
message("Clustering ", length(genes), " genes")

vst_sub <- vst[genes, , drop = FALSE]
vst_sub <- average_by_condition(vst_sub, args$metadata, args$condition_from, args$condition_order)
if (ncol(vst_sub) < 2) {
  stop("Hierarchical clustering requires at least 2 sample or condition columns; got ", ncol(vst_sub))
}
vst_sub <- zscore_rows(as.matrix(vst_sub))
if (nrow(vst_sub) < 2) stop("Need at least 2 genes with complete profiles after preprocessing")

message("Computing ", args$distance, " distance and ", args$linkage, " linkage...")
d <- compute_distance(vst_sub, args$distance)
hc <- stats::hclust(d, method = args$linkage)

saveRDS(hc, file.path(args$outdir, "hclust.rds"))

pick_k_silhouette <- function(d, hc, k_range) {
  rows <- lapply(k_range, function(k) {
    cl <- stats::cutree(hc, k = k)
    if (length(unique(cl)) < 2) return(NULL)
    sil <- cluster::silhouette(cl, d)
    data.frame(k = k, mean_silhouette = mean(sil[, "sil_width"]))
  })
  do.call(rbind, rows)
}

cut_dynamic <- function(hc, d, deep_split, min_cluster_size) {
  cl <- dynamicTreeCut::cutreeDynamic(
    dendro = hc,
    distM = as.matrix(d),
    method = "hybrid",
    deepSplit = deep_split,
    minClusterSize = min_cluster_size,
    pamRespectsDendro = FALSE,
    verbose = 0
  )
  cl <- as.integer(cl)
  # Drop unassigned (label 0) by reassigning to a new singleton-merge cluster:
  # in practice we keep them as cluster 0 and warn.
  n_unassigned <- sum(cl == 0)
  if (n_unassigned > 0) {
    message("dynamicTreeCut left ", n_unassigned, " genes unassigned (label 0); kept as cluster 0")
  }
  cl
}

if (args$num_clusters == "auto") {
  if (args$auto_method == "dynamic_tree_cut") {
    message("Auto-cutting tree with dynamicTreeCut (deepSplit=", args$deep_split,
            ", minClusterSize=", args$min_cluster_size, ")")
    clusters <- cut_dynamic(hc, d, args$deep_split, args$min_cluster_size)
    num_clusters <- length(unique(clusters))
    if (num_clusters < 2) {
      stop("dynamicTreeCut produced fewer than 2 clusters; try lowering --min_cluster_size or raising --deep_split")
    }
    # Selection summary: just record the chosen k for parity with other methods.
    write.csv(
      data.frame(method = "dynamic_tree_cut",
                 deep_split = args$deep_split,
                 min_cluster_size = args$min_cluster_size,
                 num_clusters = num_clusters),
      file.path(args$outdir, "k_selection_summary.csv"),
      row.names = FALSE
    )
  } else if (args$auto_method == "silhouette") {
    rng <- parse_cluster_range(args$cluster_range)
    candidate_k <- seq(rng[1], rng[2])
    candidate_k <- candidate_k[candidate_k < nrow(vst_sub)]
    if (length(candidate_k) == 0) stop("No valid k in --cluster_range for ", nrow(vst_sub), " genes")

    selection_df <- pick_k_silhouette(d, hc, candidate_k)
    num_clusters <- selection_df$k[which.max(selection_df$mean_silhouette)]
    selection_df$selected <- selection_df$k == num_clusters
    write.csv(selection_df, file.path(args$outdir, "k_selection_summary.csv"), row.names = FALSE)

    p_sel <- ggplot(selection_df, aes(x = k, y = mean_silhouette)) +
      geom_line(color = "#2C3E50") +
      geom_point(aes(color = selected), size = 2.5) +
      geom_vline(xintercept = num_clusters, color = "red", linetype = "dashed") +
      scale_color_manual(values = c("FALSE" = "#2C3E50", "TRUE" = "red")) +
      scale_x_continuous(breaks = selection_df$k) +
      labs(title = "Hierarchical cluster selection (mean silhouette)",
           x = "Number of clusters", y = "Mean silhouette width") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none")
    ggsave(file.path(args$outdir, "plots", "silhouette_selection.pdf"),
           plot = p_sel, width = 7, height = 4.5)

    clusters <- stats::cutree(hc, k = num_clusters)
  } else {
    stop("Unknown --auto_method: ", args$auto_method)
  }
} else {
  num_clusters <- as.integer(args$num_clusters)
  if (is.na(num_clusters) || num_clusters < 2) stop("--num_clusters must be integer >= 2 or 'auto'")
  if (num_clusters >= nrow(vst_sub)) stop("--num_clusters must be smaller than number of genes")
  message("Cutting tree at fixed k = ", num_clusters)
  clusters <- stats::cutree(hc, k = num_clusters)
}

cluster_levels <- sort(unique(clusters))
num_clusters <- length(cluster_levels)
message("Final number of clusters: ", num_clusters,
        " (labels: ", paste(cluster_levels, collapse = ","), ")")
sil <- cluster::silhouette(clusters, d)

assignments_df <- data.frame(
  gene_id = rownames(vst_sub),
  cluster = clusters,
  silhouette_width = sil[, "sil_width"],
  stringsAsFactors = FALSE
)

if (!is.null(gdb) && exists("annotate_df")) {
  assignments_df$row_index <- seq_len(nrow(assignments_df))
  assignments_df <- annotate_df(gdb, assignments_df)
  assignments_df <- assignments_df[order(assignments_df$row_index), , drop = FALSE]
  assignments_df$row_index <- NULL
}

write.csv(assignments_df, file.path(args$outdir, "cluster_assignments.csv"), row.names = FALSE)

cluster_summary <- assignments_df %>%
  group_by(cluster) %>%
  summarise(n_genes = n(),
            mean_silhouette = mean(silhouette_width, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(cluster)
write.csv(cluster_summary, file.path(args$outdir, "cluster_summary.csv"), row.names = FALSE)

for (i in cluster_levels) {
  cluster_genes <- assignments_df$gene_id[assignments_df$cluster == i]
  writeLines(cluster_genes, file.path(args$outdir, paste0("cluster_", i, "_genes.txt")))
}

# --- Dendrogram ---
pdf(file.path(args$outdir, "plots", "dendrogram.pdf"), width = 10, height = 5)
op <- par(mar = c(2, 4, 2, 1))
plot(hc, labels = FALSE, hang = -1,
     main = paste0("Hierarchical clustering (", args$linkage, " / ", args$distance, ")"),
     xlab = "", sub = "")
# rect.hclust only works for a uniform horizontal cut. dynamicTreeCut makes
# variable-depth cuts, so the boxes wouldn't match the actual clustering.
if (args$num_clusters != "auto" || args$auto_method == "silhouette") {
  stats::rect.hclust(hc, k = num_clusters, border = "red")
}
par(op)
dev.off()

# --- Heatmap with row dendrogram ---
plot_mat <- vst_sub
plot_clusters <- clusters
if (nrow(plot_mat) > args$max_plot_genes) {
  message("Downsampling heatmap to top ", args$max_plot_genes, " by variance (full table preserved)")
  vars <- apply(plot_mat, 1, stats::var)
  keep_idx <- order(vars, decreasing = TRUE)[seq_len(args$max_plot_genes)]
  plot_mat <- plot_mat[keep_idx, , drop = FALSE]
  plot_clusters <- plot_clusters[keep_idx]
}

cluster_colors <- setNames(
  grDevices::hcl.colors(num_clusters, palette = "Dark 3"),
  as.character(cluster_levels)
)
row_anno <- ComplexHeatmap::rowAnnotation(
  cluster = as.character(plot_clusters),
  col = list(cluster = cluster_colors),
  show_annotation_name = FALSE
)

heatmap_col <- circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

heatmap_distance_fn <- if (args$distance == "correlation") {
  function(m) stats::as.dist(1 - stats::cor(t(m), use = "pairwise.complete.obs"))
} else {
  function(m) stats::dist(m, method = "euclidean")
}

pdf(file.path(args$outdir, "plots", "heatmap.pdf"), width = 8, height = 10)
ht <- ComplexHeatmap::Heatmap(
  plot_mat,
  name = "z-score",
  col = heatmap_col,
  cluster_rows = TRUE,
  clustering_distance_rows = heatmap_distance_fn,
  clustering_method_rows = args$linkage,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_row_dend = TRUE,
  row_split = plot_clusters,
  left_annotation = row_anno,
  use_raster = TRUE,
  raster_quality = 2,
  column_title = paste0("Hierarchical clusters (k = ", num_clusters, ")")
)
ComplexHeatmap::draw(ht, merge_legend = TRUE)
dev.off()

# --- Per-cluster centroid profile plots ---
condition_labels <- colnames(vst_sub)
for (i in cluster_levels) {
  cluster_genes <- assignments_df$gene_id[assignments_df$cluster == i]
  cluster_expr <- vst_sub[cluster_genes, , drop = FALSE]
  center <- colMeans(cluster_expr, na.rm = TRUE)

  cluster_long <- data.frame(
    gene = rep(cluster_genes, each = length(condition_labels)),
    condition = rep(condition_labels, times = length(cluster_genes)),
    expr = as.vector(t(cluster_expr)),
    stringsAsFactors = FALSE
  )
  cluster_long$condition <- factor(cluster_long$condition, levels = condition_labels)
  center_df <- data.frame(
    condition = factor(condition_labels, levels = condition_labels),
    expr = center
  )

  p_profile <- ggplot() +
    geom_line(data = cluster_long,
              aes(x = condition, y = expr, group = gene),
              color = "grey60", linewidth = 0.3, alpha = 0.3) +
    geom_line(data = center_df, aes(x = condition, y = expr, group = 1),
              color = "red", linewidth = 1.1) +
    geom_point(data = center_df, aes(x = condition, y = expr),
               color = "red", size = 2) +
    labs(title = paste0("Hierarchical Cluster ", i, " (", length(cluster_genes), " genes)"),
         x = "Condition", y = "Standardized Expression") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(args$outdir, "plots", paste0("cluster_", i, "_profile.pdf")),
         plot = p_profile, width = 6, height = 4)
}

message("Done. Output: ", args$outdir)
