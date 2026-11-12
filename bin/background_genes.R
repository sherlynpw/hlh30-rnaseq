#!/usr/bin/env Rscript

# background_genes.R - union of gene_id across DESeq2 result CSVs = the set of
# genes actually tested in the analysis. Used as the enrichment universe so
# GO/KEGG test against expressed/tested genes, not the whole annotated genome
# (which inflates broad, vague terms).

suppressPackageStartupMessages(library(optparse))

opts <- list(
  make_option("--input_dir", type = "character", default = "."),
  make_option("--outdir",    type = "character", default = ".")
)
args <- parse_args(OptionParser(option_list = opts))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(args$input_dir, pattern = "\\.csv$", full.names = TRUE)
genes <- character(0)
for (f in files) {
  df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(df) && "gene_id" %in% names(df))
    genes <- union(genes, df$gene_id)
}
genes <- genes[!is.na(genes) & nzchar(genes)]
if (length(genes) == 0)
  stop("No gene_id values found in CSVs under ", args$input_dir)

writeLines(genes, file.path(args$outdir, "background.txt"))
message("background.txt: ", length(genes), " tested genes")
