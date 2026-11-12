#!/usr/bin/env Rscript

# ruvseq.R - Estimate factors of unwanted variation with RUVSeq
#
# Default method: RUVs (uses biological replicates of the same condition
# as anchors for unwanted variation). Outputs a TSV of W factors per sample
# that can be merged into DESeq2 design and passed as covariates to
# limma::removeBatchEffect for clustering/visualization.

suppressPackageStartupMessages({
  library(optparse)
  library(RUVSeq)
  library(EDASeq)
})

opts <- list(
  make_option("--counts",         type = "character", help = "Gene count matrix TSV (genes x samples)"),
  make_option("--metadata",       type = "character", help = "Sample metadata TSV"),
  make_option("--method",         type = "character", default = "RUVs", help = "RUVs | RUVg | RUVr"),
  make_option("--k",              type = "integer",   default = 1, help = "Number of factors of unwanted variation"),
  make_option("--condition_from", type = "character", default = "genotype_timepoint",
              help = "Column name(s) joined by underscore to define condition for anchor groups"),
  make_option("--anchor_groups",  type = "character", default = NULL,
              help = "Optional path to TSV with columns sample_id,group; samples sharing a group are anchors. Overrides condition-based grouping."),
  make_option("--min_count",      type = "integer",   default = 10, help = "Minimum total count per gene"),
  make_option("--outdir",         type = "character", default = ".", help = "Output directory")
)

args <- parse_args(OptionParser(option_list = opts))

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# Read inputs
counts <- read.delim(args$counts, row.names = 1, check.names = FALSE)
meta   <- read.delim(args$metadata, stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id

shared <- intersect(colnames(counts), meta$sample_id)
stopifnot("No matching samples between counts and metadata" = length(shared) > 0)
counts <- as.matrix(counts[, shared, drop = FALSE])
meta   <- meta[shared, , drop = FALSE]

# Filter low-count genes
keep <- rowSums(counts) >= args$min_count
counts <- counts[keep, ]
storage.mode(counts) <- "integer"
message("Genes after filter: ", nrow(counts))

# Build condition column for diagnostics and (default) anchor grouping
cond_cols <- strsplit(args$condition_from, "_")[[1]]
meta$condition <- apply(meta[, cond_cols, drop = FALSE], 1, paste, collapse = "_")

# Build replicate matrix for RUVs: rows = anchor groups, cols = sample indices
build_scIdx <- function(group_vec, sample_ids) {
  groups <- split(seq_along(sample_ids), group_vec)
  groups <- groups[lengths(groups) >= 2]
  if (length(groups) == 0) {
    stop("No anchor group has >=2 samples - RUVs needs replicate groups")
  }
  max_n <- max(lengths(groups))
  scIdx <- matrix(-1, nrow = length(groups), ncol = max_n)
  for (i in seq_along(groups)) {
    idx <- groups[[i]]
    scIdx[i, seq_along(idx)] <- idx
  }
  scIdx
}

# Resolve anchor grouping
if (!is.null(args$anchor_groups) && file.exists(args$anchor_groups)) {
  ag <- read.delim(args$anchor_groups, stringsAsFactors = FALSE)
  stopifnot("anchor_groups TSV must have sample_id,group columns" =
              all(c("sample_id", "group") %in% colnames(ag)))
  group_vec <- ag$group[match(colnames(counts), ag$sample_id)]
  group_vec[is.na(group_vec)] <- paste0("__solo_", seq_len(sum(is.na(group_vec))))
  message("Using explicit anchor_groups (", length(unique(group_vec)), " groups)")
} else {
  group_vec <- meta$condition
  message("Using condition as anchor group (", length(unique(group_vec)), " groups)")
}

set <- newSeqExpressionSet(counts, phenoData = data.frame(
  condition = meta$condition,
  row.names = colnames(counts)
))
set <- betweenLaneNormalization(set, which = "upper")

# Run chosen RUV variant
if (args$method == "RUVs") {
  scIdx <- build_scIdx(group_vec, colnames(counts))
  message("RUVs: ", nrow(scIdx), " anchor groups, k=", args$k)
  set_ruv <- RUVs(set, cIdx = rownames(counts), k = args$k, scIdx = scIdx)
} else if (args$method == "RUVg") {
  stop("RUVg not yet implemented - pass control_genes via a future flag")
} else if (args$method == "RUVr") {
  stop("RUVr not yet implemented")
} else {
  stop("Unknown --method: ", args$method)
}

W <- pData(set_ruv)[, grep("^W_", colnames(pData(set_ruv))), drop = FALSE]
W$sample_id <- rownames(W)
W <- W[, c("sample_id", setdiff(colnames(W), "sample_id"))]

# Diagnostics: correlation of each W with condition (one-way ANOVA R^2)
cond_factor <- factor(meta$condition[match(W$sample_id, meta$sample_id)])
diag_rows <- lapply(grep("^W_", colnames(W), value = TRUE), function(w) {
  r2 <- summary(lm(W[[w]] ~ cond_factor))$r.squared
  data.frame(factor = w, r2_vs_condition = r2)
})
diag <- do.call(rbind, diag_rows)
message("W vs condition R^2:")
print(diag)

write.table(W,    file = file.path(args$outdir, "w_factors.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(diag, file = file.path(args$outdir, "ruv_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Diagnostic plots
pdf(file.path(args$outdir, "ruv_diagnostics.pdf"), width = 9, height = 7)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

cols_cond <- as.integer(cond_factor)

plotRLE(set,     outline = FALSE, col = cols_cond, main = "RLE - before RUV")
plotRLE(set_ruv, outline = FALSE, col = cols_cond, main = "RLE - after RUV")
plotPCA(set,     col = cols_cond, cex = 1.0, main = "PCA - before RUV")
plotPCA(set_ruv, col = cols_cond, cex = 1.0, main = "PCA - after RUV")

par(op)
dev.off()

message("Done. W factors: ", file.path(args$outdir, "w_factors.tsv"))
