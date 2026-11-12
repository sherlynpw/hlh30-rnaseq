#!/usr/bin/env Rscript

# cv_density.R - Coefficient of variation density plots
#
# Input:  VST counts CSV
# Output: CV density plot (one curve per condition)

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
})

opts <- list(
  make_option("--counts",   type = "character", help = "VST counts CSV (genes x samples)"),
  make_option("--metadata",       type = "character", default = NULL, help = "Sample metadata TSV (optional, for per-condition CV)"),
  make_option("--condition_from", type = "character", default = "genotype_timepoint", help = "Underscore-joined metadata columns for condition"),
  make_option("--outdir",   type = "character", default = ".", help = "Output directory")
)

args <- parse_args(OptionParser(option_list = opts))

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# Read VST counts
vst <- read.csv(args$counts, row.names = 1, check.names = FALSE)

compute_cv <- function(mat) {
  # CV = sd / mean * 100 for each gene (row)
  means <- rowMeans(mat)
  sds   <- apply(mat, 1, sd)
  # Avoid division by zero
  cv <- ifelse(means > 0, (sds / means) * 100, NA)
  cv[!is.na(cv)]
}

if (!is.null(args$metadata) && file.exists(args$metadata)) {
  # Per-condition CV
  meta <- read.delim(args$metadata, stringsAsFactors = FALSE)
  cond_cols <- strsplit(args$condition_from, "_")[[1]]
  meta$condition <- do.call(paste, c(meta[, cond_cols, drop = FALSE], sep = "_"))

  cv_list <- list()
  for (cond in unique(meta$condition)) {
    samples <- meta$sample_id[meta$condition == cond]
    samples <- intersect(samples, colnames(vst))
    if (length(samples) < 2) next
    cv_vals <- compute_cv(vst[, samples, drop = FALSE])
    cv_list[[cond]] <- data.frame(cv = cv_vals, condition = cond)
  }
  cv_df <- bind_rows(cv_list)

  p <- ggplot(cv_df, aes(x = cv, color = condition)) +
    geom_density(linewidth = 0.8) +
    labs(
      x     = "Coefficient of Variation (%)",
      y     = "Density",
      title = "Gene Expression CV Density by Condition",
      color = "Condition"
    ) +
    theme_minimal(base_size = 12) +
    xlim(0, quantile(cv_df$cv, 0.99))
} else {
  # Single CV across all samples
  cv_vals <- compute_cv(vst)
  cv_df   <- data.frame(cv = cv_vals)

  p <- ggplot(cv_df, aes(x = cv)) +
    geom_density(linewidth = 0.8, fill = "steelblue", alpha = 0.3) +
    labs(
      x     = "Coefficient of Variation (%)",
      y     = "Density",
      title = "Gene Expression CV Density"
    ) +
    theme_minimal(base_size = 12) +
    xlim(0, quantile(cv_df$cv, 0.99))
}

ggsave(
  file.path(args$outdir, "cv_density_plot.pdf"),
  plot = p, width = 8, height = 5
)

write.csv(cv_df, file.path(args$outdir, "cv_data.csv"), row.names = FALSE)

message("Done. Output: ", args$outdir)
