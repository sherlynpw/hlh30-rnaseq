#!/usr/bin/env Rscript

# sample_correlation.R - Pairwise sample-sample correlation heatmap.
#
# More robust than PCA for "are these two groups similar" questions: PCA only
# shows dominant axes of variance, which can inflate apparent distance when
# a strong non-biological axis (lab, protocol) dominates PC1/PC2 even if the
# underlying expression vectors are highly concordant. Pairwise correlation
# directly measures concordance gene-by-gene.
#
# Inputs:
#   --vst       VST-normalized counts (matrix, rows = genes, cols = samples)
#   --metadata  Sample meta TSV (sample_id + factor columns)
#   --top_n     Use top-N most-variable genes (default 2000) - whole matrix
#               makes the correlation dominated by housekeeping flat structure
#   --method    Pearson or Spearman (default spearman; rank-based, robust to
#               outliers and library-size-normalization wobbles)
#   --outdir
#
# Outputs:
#   sample_correlation.csv            full pairwise matrix
#   sample_correlation.pdf            clustered heatmap with values overlaid
#   sample_correlation_summary.tsv    per-(group_a, group_b) mean correlation
#                                     (uses condition by default for grouping)

suppressPackageStartupMessages({
  library(optparse)
  library(pheatmap)
  library(RColorBrewer)
})

opts <- list(
  make_option("--vst",      type = "character", help = "VST counts CSV"),
  make_option("--metadata", type = "character", help = "Sample metadata TSV/CSV"),
  make_option("--top_n",    type = "integer",   default = 0,
              help = "Use top-N genes by --selection (0 = all genes; default 0). When >0 with --selection=variable, selection-induced anti-correlation can occur if your samples form distinct groups (the top-variable genes are by definition the ones differentiating the groups)."),
  make_option("--selection", type = "character", default = "expressed",
              help = "How to pick the top-N genes when top_n>0: 'expressed' (highest mean expression, no group bias) | 'variable' (highest variance, biased toward inter-group differences if groups exist). Default: expressed."),
  make_option("--method",   type = "character", default = "spearman",
              help = "Correlation method: pearson | spearman (default spearman)"),
  make_option("--group_by", type = "character", default = "condition",
              help = "Metadata column to group samples for summary table (default condition)"),
  make_option("--outdir",   type = "character", default = ".",
              help = "Output directory")
)
args <- parse_args(OptionParser(option_list = opts))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# Read VST counts (rows = genes, cols = samples)
vst <- as.matrix(read.csv(args$vst, row.names = 1, check.names = FALSE))

# Read metadata (try TSV first, fall back to CSV)
meta <- tryCatch(read.delim(args$metadata, stringsAsFactors = FALSE),
                 error = function(e) read.csv(args$metadata, stringsAsFactors = FALSE))
rownames(meta) <- meta$sample_id

# Align meta rows to VST columns
common <- intersect(colnames(vst), rownames(meta))
vst  <- vst[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]

# Gene selection
if (args$top_n > 0) {
  if (args$selection == "variable") {
    score <- apply(vst, 1, var)
    sel_name <- "variable"
  } else if (args$selection == "expressed") {
    score <- rowMeans(vst)
    sel_name <- "expressed"
  } else {
    stop("Unknown --selection: ", args$selection, " (use 'variable' or 'expressed')")
  }
  top_idx <- order(score, decreasing = TRUE)[seq_len(min(args$top_n, length(score)))]
  vst_top <- vst[top_idx, , drop = FALSE]
  message(sprintf("Using top %d / %d most %s genes for correlation",
                  length(top_idx), nrow(vst), sel_name))
} else {
  vst_top <- vst
  message(sprintf("Using all %d genes for correlation", nrow(vst)))
}

# Pairwise correlation
cm <- cor(vst_top, method = args$method)
write.csv(cm, file.path(args$outdir, "sample_correlation.csv"))

# Heatmap. Annotation: any sensible factor column we can find.
ann_cols <- intersect(c("condition", "timepoint", "source", "genotype"),
                      colnames(meta))
ann_df <- if (length(ann_cols) > 0) meta[, ann_cols, drop = FALSE] else NULL

pdf(file.path(args$outdir, "sample_correlation.pdf"),
    width = max(7, 0.4 * ncol(cm) + 3),
    height = max(6, 0.35 * ncol(cm) + 2))
pheatmap(cm,
         color = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 7,
         annotation_row = ann_df,
         annotation_col = ann_df,
         main = sprintf("Sample-sample %s correlation (%d genes)",
                        args$method, nrow(vst_top)))
invisible(dev.off())

# Per-group-pair summary: mean off-diagonal correlation between group_a and
# group_b samples. Useful to check that two groups look similar.
if (args$group_by %in% colnames(meta)) {
  groups <- meta[[args$group_by]]
  ug <- sort(unique(groups))
  rows <- list()
  for (i in seq_along(ug)) {
    for (j in i:length(ug)) {
      ga <- ug[i]; gb <- ug[j]
      si <- which(groups == ga)
      sj <- which(groups == gb)
      vals <- cm[si, sj, drop = FALSE]
      if (i == j) {
        # within-group: take upper triangle only
        if (length(si) >= 2) {
          vals <- vals[upper.tri(vals)]
        } else {
          vals <- NA_real_
        }
      } else {
        vals <- as.vector(vals)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        group_a = ga, group_b = gb,
        n_pairs = length(vals),
        mean    = if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE),
        median  = if (all(is.na(vals))) NA_real_ else median(vals, na.rm = TRUE),
        min     = if (all(is.na(vals))) NA_real_ else min(vals, na.rm = TRUE),
        max     = if (all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  summary_df <- do.call(rbind, rows)
  write.table(summary_df,
              file = file.path(args$outdir, "sample_correlation_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  message("Per-group summary (mean correlation):")
  for (r in seq_len(nrow(summary_df))) {
    message(sprintf("  %-20s vs %-20s  n=%d  mean=%.3f",
                    summary_df$group_a[r], summary_df$group_b[r],
                    summary_df$n_pairs[r], summary_df$mean[r]))
  }
} else {
  message("Skipping group summary: column '", args$group_by, "' not in metadata")
}

message("Done. Output: ", args$outdir)
