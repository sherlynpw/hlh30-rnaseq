#!/usr/bin/env Rscript

# kmeans_clustering.R - K-means clustering on VST counts
#
# Input:  VST counts + filtered gene list
# Output: Cluster assignments, elbow-based k selection, plots

suppressPackageStartupMessages({
  library(optparse)
  library(cluster)
  library(ggplot2)
  library(dplyr)
  library(uwot)
  library(RSQLite)
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
  make_option("--nstart", type = "integer", default = 25, help = "Number of random starts for kmeans"),
  make_option("--iter_max", type = "integer", default = 100, help = "Maximum number of iterations for kmeans"),
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
  if (is.null(metadata_path)) {
    return(mat)
  }
  if (!file.exists(metadata_path)) {
    stop("Metadata file not found: ", metadata_path)
  }

  meta <- read.delim(metadata_path, stringsAsFactors = FALSE)
  cond_cols <- strsplit(condition_from, "_", fixed = TRUE)[[1]]

  if (!"sample_id" %in% colnames(meta)) {
    stop("Metadata must contain a sample_id column")
  }
  missing_cols <- setdiff(cond_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Metadata missing condition columns: ", paste(missing_cols, collapse = ", "))
  }

  meta$condition <- do.call(paste, c(meta[, cond_cols, drop = FALSE], sep = "_"))

  if (!is.null(condition_order) && nchar(condition_order) > 0) {
    requested <- strsplit(condition_order, ",", fixed = TRUE)[[1]]
    missing_conditions <- setdiff(requested, unique(meta$condition))
    if (length(missing_conditions) > 0) {
      stop(
        "Requested conditions in --condition_order were not found in metadata: ",
        paste(missing_conditions, collapse = ", ")
      )
    }
    conditions <- requested
  } else {
    conditions <- sort(unique(meta$condition))
  }

  samples_by_condition <- lapply(conditions, function(cond) {
    intersect(meta$sample_id[meta$condition == cond], colnames(mat))
  })
  missing_sample_conditions <- conditions[vapply(samples_by_condition, length, integer(1)) == 0]
  if (length(missing_sample_conditions) > 0) {
    stop(
      "Metadata conditions had no matching samples in the expression matrix: ",
      paste(missing_sample_conditions, collapse = ", ")
    )
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
  if (dropped > 0) {
    message("Dropping ", dropped, " genes with incomplete standardized profiles")
  }
  scaled[keep, , drop = FALSE]
}

parse_cluster_range <- function(range_text) {
  parts <- as.integer(strsplit(range_text, ",", fixed = TRUE)[[1]])
  if (length(parts) != 2 || any(is.na(parts))) {
    stop("--cluster_range must be formatted as min,max")
  }
  if (parts[1] < 2 || parts[2] < parts[1]) {
    stop("--cluster_range must satisfy min >= 2 and max >= min")
  }
  parts
}

run_kmeans <- function(mat, k, seed, nstart, iter_max) {
  set.seed(seed)
  stats::kmeans(mat, centers = k, nstart = nstart, iter.max = iter_max)
}

compute_embedding <- function(mat, seed) {
  n_points <- nrow(mat)
  if (n_points >= 3) {
    n_neighbors <- min(15, n_points - 1)
    set.seed(seed)
    embed <- tryCatch(
      uwot::umap(mat, n_neighbors = n_neighbors, min_dist = 0.1, ret_model = FALSE, verbose = FALSE),
      error = function(e) NULL
    )
    if (!is.null(embed) && ncol(embed) >= 2) {
      return(embed[, 1:2, drop = FALSE])
    }
  }

  if (n_points >= 2) {
    pcs <- stats::prcomp(mat, center = FALSE, scale. = FALSE)$x
    if (ncol(pcs) == 1) {
      return(cbind(Dim1 = pcs[, 1], Dim2 = rep(0, n_points)))
    }
    return(pcs[, 1:2, drop = FALSE])
  }

  cbind(Dim1 = 0, Dim2 = 0)
}

vst <- read.csv(args$vst, row.names = 1, check.names = FALSE)
genes <- intersect(read_gene_list(args$genes), rownames(vst))

if (length(genes) < 2) {
  stop("Need at least 2 intersecting genes for k-means clustering")
}

message("Clustering ", length(genes), " genes")

vst_sub <- vst[genes, , drop = FALSE]
vst_sub <- average_by_condition(vst_sub, args$metadata, args$condition_from, args$condition_order)
if (ncol(vst_sub) < 2) {
  stop(
    "K-means clustering requires at least 2 samples or condition profiles after preprocessing; ",
    "got ", ncol(vst_sub), ". Add more samples/conditions or adjust metadata grouping."
  )
}
vst_sub <- zscore_rows(as.matrix(vst_sub))

if (nrow(vst_sub) < 2) {
  stop("Need at least 2 genes with complete profiles after preprocessing")
}
if (ncol(vst_sub) < 1) {
  stop("No expression columns available after preprocessing")
}

distance_matrix <- stats::dist(vst_sub)

pick_elbow_k <- function(selection_df) {
  if (nrow(selection_df) == 1) {
    return(selection_df$k[[1]])
  }
  if (nrow(selection_df) == 2) {
    return(selection_df$k[[2]])
  }

  x <- selection_df$k
  y <- selection_df$tot_withinss

  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))

  start_point <- c(x_norm[[1]], y_norm[[1]])
  end_point <- c(x_norm[[length(x_norm)]], y_norm[[length(y_norm)]])
  line_vec <- end_point - start_point
  line_len <- sqrt(sum(line_vec^2))

  if (line_len == 0) {
    return(selection_df$k[[1]])
  }

  distances <- vapply(seq_along(x_norm), function(i) {
    point <- c(x_norm[[i]], y_norm[[i]])
    rel_point <- point - start_point
    abs(line_vec[[1]] * rel_point[[2]] - line_vec[[2]] * rel_point[[1]]) / line_len
  }, numeric(1))

  selection_df$k[[which.max(distances)]]
}

if (args$num_clusters == "auto") {
  range_parts <- parse_cluster_range(args$cluster_range)
  candidate_k <- seq(range_parts[1], range_parts[2])
  candidate_k <- candidate_k[candidate_k < nrow(vst_sub)]

  if (length(candidate_k) == 0) {
    stop("No valid k values in --cluster_range for ", nrow(vst_sub), " genes")
  }

  selection_rows <- lapply(candidate_k, function(k) {
    km <- run_kmeans(vst_sub, k, args$seed, args$nstart, args$iter_max)
    data.frame(
      k = k,
      tot_withinss = km$tot.withinss,
      stringsAsFactors = FALSE
    )
  })
  selection_df <- bind_rows(selection_rows) %>% arrange(k)

  num_clusters <- pick_elbow_k(selection_df)
  selection_df <- selection_df %>% mutate(selected = k == num_clusters)
  write.csv(selection_df, file.path(args$outdir, "k_selection_summary.csv"), row.names = FALSE)

  p_sel <- ggplot(selection_df, aes(x = k, y = tot_withinss)) +
    geom_line(color = "#2C3E50") +
    geom_point(aes(color = selected), size = 2.5) +
    geom_vline(xintercept = num_clusters, color = "red", linetype = "dashed") +
    scale_color_manual(values = c("FALSE" = "#2C3E50", "TRUE" = "red")) +
    scale_x_continuous(breaks = selection_df$k) +
    labs(
      title = "K-means cluster selection by elbow",
      x = "Number of clusters",
      y = "Total within-cluster sum of squares"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

  ggsave(
    file.path(args$outdir, "plots", "elbow_selection.pdf"),
    plot = p_sel,
    width = 7,
    height = 4.5
  )
} else {
  num_clusters <- as.integer(args$num_clusters)
  if (is.na(num_clusters) || num_clusters < 2) {
    stop("--num_clusters must be an integer >= 2 or 'auto'")
  }
  if (num_clusters >= nrow(vst_sub)) {
    stop("--num_clusters must be smaller than the number of clustered genes")
  }

  km_fixed <- run_kmeans(vst_sub, num_clusters, args$seed, args$nstart, args$iter_max)
  selection_df <- data.frame(
    k = num_clusters,
    tot_withinss = km_fixed$tot.withinss,
    selected = TRUE,
    stringsAsFactors = FALSE
  )
  write.csv(selection_df, file.path(args$outdir, "k_selection_summary.csv"), row.names = FALSE)

  p_sel <- ggplot(selection_df, aes(x = k, y = tot_withinss)) +
    geom_point(color = "red", size = 2.5) +
    geom_vline(xintercept = num_clusters, color = "red", linetype = "dashed") +
    scale_x_continuous(breaks = selection_df$k) +
    labs(
      title = "K-means cluster selection by elbow",
      x = "Number of clusters",
      y = "Total within-cluster sum of squares"
    ) +
    theme_minimal(base_size = 12)

  ggsave(
    file.path(args$outdir, "plots", "elbow_selection.pdf"),
    plot = p_sel,
    width = 7,
    height = 4.5
  )
}

message("Running k-means with k = ", num_clusters, "...")
km <- run_kmeans(vst_sub, num_clusters, args$seed, args$nstart, args$iter_max)
silhouette_info <- cluster::silhouette(km$cluster, distance_matrix)

assignments_df <- data.frame(
  gene_id = rownames(vst_sub),
  cluster = km$cluster,
  silhouette_width = silhouette_info[, "sil_width"],
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
  summarise(
    n_genes = n(),
    mean_silhouette = mean(silhouette_width, na.rm = TRUE),
    withinss = km$withinss[as.integer(first(cluster))],
    .groups = "drop"
  ) %>%
  arrange(cluster)

write.csv(cluster_summary, file.path(args$outdir, "cluster_summary.csv"), row.names = FALSE)

for (i in seq_len(num_clusters)) {
  cluster_genes <- assignments_df$gene_id[assignments_df$cluster == i]
  writeLines(cluster_genes, file.path(args$outdir, paste0("cluster_", i, "_genes.txt")))
}

embedding <- compute_embedding(vst_sub, args$seed)
plot_df <- data.frame(
  UMAP1 = embedding[, 1],
  UMAP2 = embedding[, 2],
  cluster = factor(assignments_df$cluster, levels = seq_len(num_clusters))
)

p_umap <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 1.2, alpha = 0.8) +
  labs(
    title = paste0("K-means clusters (k = ", num_clusters, ")"),
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(args$outdir, "plots", "kmeans_umap.pdf"),
  plot = p_umap,
  width = 8,
  height = 6
)

plot_mat <- vst_sub
condition_labels <- colnames(plot_mat)

for (i in seq_len(num_clusters)) {
  cluster_genes <- assignments_df$gene_id[assignments_df$cluster == i]
  cluster_expr <- plot_mat[cluster_genes, , drop = FALSE]
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
    geom_line(
      data = cluster_long,
      aes(x = condition, y = expr, group = gene),
      color = "grey60",
      linewidth = 0.3,
      alpha = 0.3
    ) +
    geom_line(
      data = center_df,
      aes(x = condition, y = expr, group = 1),
      color = "red",
      linewidth = 1.1
    ) +
    geom_point(data = center_df, aes(x = condition, y = expr), color = "red", size = 2) +
    labs(
      title = paste0("K-means Cluster ", i, " (", length(cluster_genes), " genes)"),
      x = "Condition",
      y = "Standardized Expression"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    file.path(args$outdir, "plots", paste0("cluster_", i, "_profile.pdf")),
    plot = p_profile,
    width = 6,
    height = 4
  )
}

message("Done. Output: ", args$outdir)
