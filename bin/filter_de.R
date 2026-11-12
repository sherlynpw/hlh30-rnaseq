#!/usr/bin/env Rscript

# filter_de.R - Filter DESeq2 results by padj and log2FC thresholds
#
# Reads all CSVs from input directory, filters each, writes flat gene list files.
# Enrichment scripts downstream find all *genes*.txt files recursively.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})

opts <- list(
  make_option("--input_dir", type = "character", default = ".", help = "Directory with DESeq2 result CSVs"),
  make_option("--padj",      type = "numeric",   default = 0.05, help = "Adjusted p-value threshold"),
  make_option("--lfc",       type = "numeric",   default = 1.0,  help = "Absolute log2FC threshold"),
  make_option("--outdir",    type = "character",  default = "filtered", help = "Output directory")
)

args <- parse_args(OptionParser(option_list = opts))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

result_files <- list.files(args$input_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot("No result files found" = length(result_files) > 0)

all_de_genes <- c()

for (f in result_files) {
  comparison <- tools::file_path_sans_ext(basename(f))
  res <- read.csv(f, stringsAsFactors = FALSE)

  de <- res %>% filter(!is.na(padj), padj < args$padj, abs(log2FoldChange) >= args$lfc)
  de_up   <- de %>% filter(log2FoldChange > 0) %>% arrange(desc(log2FoldChange))
  de_down <- de %>% filter(log2FoldChange < 0) %>% arrange(log2FoldChange)

  # Flat filenames - each comparison + direction is unique
  writeLines(de$gene_id,      file.path(args$outdir, paste0(comparison, "_genes_de.txt")))
  writeLines(de_up$gene_id,   file.path(args$outdir, paste0(comparison, "_genes_up.txt")))
  writeLines(de_down$gene_id, file.path(args$outdir, paste0(comparison, "_genes_down.txt")))
  write.csv(de, file.path(args$outdir, paste0(comparison, "_filtered_results.csv")), row.names = FALSE)

  all_de_genes <- union(all_de_genes, de$gene_id)
  message(comparison, ": ", nrow(de), " DE (", nrow(de_up), " up, ", nrow(de_down), " down)")
}

writeLines(all_de_genes, file.path(args$outdir, "all_de_genes.txt"))

writeLines(c(
  paste("Comparisons:", length(result_files)),
  paste("padj threshold:", args$padj),
  paste("lfc threshold:", args$lfc)
), file.path(args$outdir, "filter_summary.txt"))

message("Total unique DE genes: ", length(all_de_genes))
message("Done. Output: ", args$outdir)
