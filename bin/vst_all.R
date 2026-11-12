#!/usr/bin/env Rscript

# vst_all.R - all-sample VST for PCA / exploration in grouped (split) runs.
#
# Design ~1: no DE testing, just one variance-stabilized matrix across every
# sample on a common scale. The cross-lab batch effect is intentionally
# retained - this matrix is for PCA inspection, not statistics.

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
})

opts <- list(
  make_option("--counts",    type = "character", help = "Gene count matrix TSV"),
  make_option("--metadata",  type = "character", help = "Sample metadata TSV"),
  make_option("--min_count", type = "integer",   default = 10),
  make_option("--outdir",    type = "character", default = ".")
)
args <- parse_args(OptionParser(option_list = opts))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

counts <- read.delim(args$counts, row.names = 1, check.names = FALSE)
meta   <- read.delim(args$metadata, stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id

shared <- intersect(colnames(counts), meta$sample_id)
stopifnot("No matching samples between counts and metadata" = length(shared) > 0)
counts <- counts[, shared, drop = FALSE]
meta   <- meta[shared, , drop = FALSE]
counts <- counts[rowSums(counts) >= args$min_count, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(round(as.matrix(counts)), colData = meta, design = ~ 1)
dds <- estimateSizeFactors(dds)
v <- tryCatch(vst(dds, blind = TRUE),
              error = function(e) varianceStabilizingTransformation(dds, blind = TRUE))
m <- assay(v)

# No batch correction here, so corrected == raw
write.csv(m, file.path(args$outdir, "vst_counts.csv"))
write.csv(m, file.path(args$outdir, "vst_counts_raw.csv"))
message("vst_all: ", nrow(m), " genes x ", ncol(m), " samples")
