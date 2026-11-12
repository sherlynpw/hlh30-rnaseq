#!/usr/bin/env Rscript
# Split compare_de's "differential" genes by the sign of the two age changes.
# compare_de names them "higher_in_<contrast>" by the larger signed log2FC, so
# a gene that falls in both genotypes but falls more in the mutant counts as
# "higher in WT". Here a gene is "more up" in the genotype with the larger
# log2FC when both rise, "more down" in the genotype with the more negative
# log2FC when both fall, and "opposite" when the signs differ.
#
# Fig 9F/9G use more_up_in_<WT> plus opposite_up_in_<WT> (the genes that rise
# more in WT, or rise in WT and fall in the mutant), written as
# *_genes_differential_up_more_in_<WT>_or_opposite.txt.
#
#   Rscript split-differential.R --compared out/<run>/<pair>_compared_genes.csv

suppressPackageStartupMessages({ library(optparse) })
opt <- parse_args(OptionParser(option_list = list(
  make_option("--compared", type = "character", help = "*_compared_genes.csv from compare_de.R"),
  make_option("--outdir",   type = "character", default = NULL, help = "default: the folder of --compared")
)))
outdir <- if (is.null(opt$outdir)) dirname(opt$compared) else opt$outdir

d <- read.csv(opt$compared, stringsAsFactors = FALSE, check.names = FALSE)
d <- d[d$category == "differential", ]
contrasts <- sub("^log2FoldChange\\.", "", grep("^log2FoldChange\\.", names(d), value = TRUE))
wt <- grep("^(WT|N2)_", contrasts, value = TRUE)
stopifnot(length(contrasts) == 2, length(wt) == 1)
mut <- setdiff(contrasts, wt)
w <- d[[paste0("log2FoldChange.", wt)]]; m <- d[[paste0("log2FoldChange.", mut)]]

pair <- sub("_compared_genes\\.csv$", "", basename(opt$compared))
sets <- list()
sets[[paste0("more_up_in_", wt)]]     <- d$gene_id[w > 0 & m > 0 & w > m]
sets[[paste0("more_up_in_", mut)]]    <- d$gene_id[w > 0 & m > 0 & m > w]
sets[[paste0("more_down_in_", wt)]]   <- d$gene_id[w < 0 & m < 0 & w < m]
sets[[paste0("more_down_in_", mut)]]  <- d$gene_id[w < 0 & m < 0 & m < w]
sets[[paste0("opposite_up_in_", wt)]] <- d$gene_id[w > 0 & m < 0]
sets[[paste0("opposite_up_in_", mut)]] <- d$gene_id[w < 0 & m > 0]
stopifnot(sum(lengths(sets)) == nrow(d))
sets[[paste0("up_more_in_", wt, "_or_opposite")]] <-
  c(sets[[paste0("more_up_in_", wt)]], sets[[paste0("opposite_up_in_", wt)]])

for (nm in names(sets)) {
  if (length(sets[[nm]]) == 0) { message("    0  ", nm, " (not written)"); next }
  f <- file.path(outdir, sprintf("%s_genes_differential_%s.txt", pair, nm))
  writeLines(sets[[nm]], f)
  message(sprintf("%5d  %s", length(sets[[nm]]), basename(f)))
}
