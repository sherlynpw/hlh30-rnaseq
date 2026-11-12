#!/usr/bin/env Rscript

# fuzzy_clustering.R - Mfuzz soft fuzzy c-means clustering
#
# Input:  VST counts + filtered gene list
# Output: Cluster assignments, core genes, time-series plots

suppressPackageStartupMessages({
  library(optparse)
  library(Mfuzz)
  library(Biobase)
  library(dplyr)
  library(ggplot2)
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
  make_option("--num_clusters", type = "character", default = "auto", help = "Number of clusters or 'auto'"),
  make_option("--cluster_range", type = "character", default = "2,12", help = "Min,max for auto cluster selection"),
  make_option("--sort_by",   type = "character", default = "pvalue", help = "Sort criterion"),
  make_option("--min_membership",  type = "numeric",  default = 0.5, help = "Membership threshold for core genes"),
  make_option("--condition_from",  type = "character", default = "genotype_timepoint", help = "Underscore-joined metadata columns for condition (e.g. genotype_stage)"),
  make_option("--condition_order", type = "character", default = NULL, help = "Comma-separated condition order for plots"),
  make_option("--seed",           type = "integer",   default = 42, help = "Random seed"),
  make_option("--gene_db",        type = "character",  default = NULL, help = "Path to gene_db.sqlite"),
  make_option("--outdir",         type = "character",  default = ".", help = "Output directory")
)

args <- parse_args(OptionParser(option_list = opts))

gdb <- NULL
if (!is.null(args$gene_db) && file.exists(args$gene_db)) {
  gdb <- open_gene_db(args$gene_db)
}

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

# Subset to DE genes
vst_sub <- vst[genes, , drop = FALSE]

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
  vst_sub <- as.data.frame(vst_avg)
}

# Create ExpressionSet for Mfuzz
mat <- as.matrix(vst_sub)
eset <- new("ExpressionSet", exprs = mat)

# Standardize
eset <- standardise(eset)

# Determine number of clusters
if (args$num_clusters == "auto") {
  range_parts <- as.integer(strsplit(args$cluster_range, ",")[[1]])
  k_min <- range_parts[1]
  k_max <- range_parts[2]

  message("Auto-selecting clusters (range: ", k_min, "-", k_max, ")...")

  # Use Dmin method
  dmin_results <- Dmin(eset, m = 1.25, crange = seq(k_min, k_max, 1), repeats = 3)

  # Find elbow: largest drop in Dmin
  dmin_diffs <- -diff(dmin_results)
  k_optimal <- which.max(dmin_diffs) + k_min

  message("Selected k = ", k_optimal)

  # Save selection plot
  pdf(file.path(args$outdir, "plots", "dmin_plot.pdf"), width = 6, height = 4)
  plot(k_min:k_max, dmin_results, type = "b", xlab = "Number of clusters",
       ylab = "Dmin", main = "Cluster Number Selection")
  abline(v = k_optimal, col = "red", lty = 2)
  dev.off()

  num_clusters <- k_optimal
} else {
  num_clusters <- as.integer(args$num_clusters)
}

# Run fuzzy c-means
message("Running Mfuzz with k = ", num_clusters, "...")
cl <- mfuzz(eset, c = num_clusters, m = 1.25)

# Extract results
membership <- cl$membership
cluster_assignments <- apply(membership, 1, which.max)

# Build results table
results_df <- data.frame(
  gene_id    = rownames(membership),
  cluster    = cluster_assignments,
  max_membership = apply(membership, 1, max),
  stringsAsFactors = FALSE
)
results_df <- cbind(results_df, as.data.frame(membership))

# Annotate
if (!is.null(gdb)) results_df <- annotate_df(gdb, results_df)

write.csv(results_df, file.path(args$outdir, "cluster_assignments.csv"), row.names = FALSE)

# Core genes (membership >= threshold)
core_genes <- results_df %>% filter(max_membership >= args$min_membership)
write.csv(core_genes, file.path(args$outdir, "core_genes.csv"), row.names = FALSE)

message("Core genes (membership >= ", args$min_membership, "): ", nrow(core_genes))

# Generate cluster profile plots
centers <- cl$centers
time_labels <- colnames(centers)

for (i in 1:num_clusters) {
  cluster_genes <- results_df$gene_id[results_df$cluster == i]
  n_genes <- length(cluster_genes)
  n_core  <- sum(results_df$cluster == i & results_df$max_membership >= args$min_membership)

  # Plot cluster center + member gene traces
  cluster_data <- as.data.frame(t(exprs(eset)[cluster_genes, , drop = FALSE]))
  cluster_data$time <- time_labels
  cluster_long <- tidyr::pivot_longer(cluster_data, -time, names_to = "gene", values_to = "expr")
  cluster_long$membership <- membership[cluster_long$gene, i]
  cluster_long$time <- factor(cluster_long$time, levels = time_labels)

  center_df <- data.frame(
    time = factor(time_labels, levels = time_labels),
    expr = centers[i, ]
  )

  p <- ggplot() +
    geom_line(data = cluster_long, aes(x = time, y = expr, group = gene, alpha = membership),
              color = "grey60", linewidth = 0.3) +
    geom_line(data = center_df, aes(x = time, y = expr), color = "red", linewidth = 1.2) +
    geom_point(data = center_df, aes(x = time, y = expr), color = "red", size = 2) +
    scale_alpha_continuous(range = c(0.05, 0.5)) +
    labs(
      title = paste0("Cluster ", i, " (", n_genes, " genes, ", n_core, " core)"),
      x = "Condition", y = "Standardized Expression"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

  ggsave(file.path(args$outdir, "plots", paste0("cluster_", i, ".pdf")),
         plot = p, width = 6, height = 4)
}

# Per-cluster gene lists
for (i in 1:num_clusters) {
  cluster_genes <- results_df$gene_id[results_df$cluster == i]
  core <- results_df$gene_id[results_df$cluster == i & results_df$max_membership >= args$min_membership]
  writeLines(cluster_genes, file.path(args$outdir, paste0("cluster_", i, "_genes.txt")))
  writeLines(core, file.path(args$outdir, paste0("cluster_", i, "_core_genes.txt")))
}

message("Done. ", num_clusters, " clusters in ", args$outdir)
