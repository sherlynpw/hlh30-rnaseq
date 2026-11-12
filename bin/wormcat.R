#!/usr/bin/env Rscript

# wormcat.R - WormCat enrichment
# Accepts either DESeq2 result CSVs (filters internally) or gene list .txt files

suppressPackageStartupMessages({
  library(optparse)
  library(wormcat)
  library(dplyr)
  library(RSQLite)
  library(parallel)
})

.find_lib <- function() {
  for (p in c(Sys.which("gene_db_lib.R"), file.path("bin", "gene_db_lib.R")))
    if (nchar(p) > 0 && file.exists(p)) return(p)
  NULL
}
lib <- .find_lib()
if (!is.null(lib)) source(lib)

opts <- list(
  make_option("--input",     type = "character", help = "Directory with DESeq2 CSVs or gene list .txt files"),
  make_option("--filter_by", type = "character", default = NULL, help = "Filter column: padj, pvalue, or NULL"),
  make_option("--threshold", type = "numeric",   default = 0.05, help = "Filter threshold"),
  make_option("--lfc",       type = "numeric",   default = 0, help = "Min |log2FC| (0 = no filter)"),
  make_option("--cpus",      type = "integer",   default = 1, help = "Parallel workers (mclapply)"),
  make_option("--outdir",    type = "character",  default = ".", help = "Output directory")
)

args <- parse_args(OptionParser(option_list = opts))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

gene_lists <- find_gene_lists(args$input, args$filter_by, args$threshold, args$lfc)

if (length(gene_lists) == 0) {
  write.csv(data.frame(message = "No gene lists found"),
            file.path(args$outdir, "wormcat_no_results.csv"), row.names = FALSE)
  quit(save = "no", status = 0)
}

message("Processing ", length(gene_lists), " gene lists with ", args$cpus, " worker(s)")

process_one <- function(label) {
  genes <- gene_lists[[label]]
  out_sub <- file.path(args$outdir, label)

  if (length(genes) < 3) {
    dir.create(out_sub, showWarnings = FALSE)
    write.csv(data.frame(message = paste("Too few genes:", length(genes))),
              file.path(out_sub, "wormcat_skipped.csv"), row.names = FALSE)
    return(invisible(NULL))
  }

  message("WormCat: ", label, " (", length(genes), " genes)")

  tmp_input <- tempfile(fileext = ".csv")
  write.csv(data.frame(Wormbase.ID = genes), tmp_input, row.names = FALSE)

  tryCatch({
    worm_cat_fun(
      file_to_process = tmp_input,
      title           = label,
      output_dir      = out_sub,
      input_type      = "Wormbase.ID",
      zip_files       = FALSE
    )
    out_files <- list.files(out_sub, recursive = TRUE)
    message("  [", label, "] Done - ", length(out_files), " files")
  }, error = function(e) {
    message("  [", label, "] ERROR: ", e$message)
    dir.create(out_sub, showWarnings = FALSE)
    writeLines(e$message, file.path(out_sub, "ERROR.txt"))
  })
  invisible(NULL)
}

n_workers <- max(1L, as.integer(args$cpus))
if (n_workers > 1L) {
  mclapply(names(gene_lists), process_one, mc.cores = n_workers, mc.preschedule = FALSE)
} else {
  lapply(names(gene_lists), process_one)
}

if (length(list.files(args$outdir, pattern = "\\.csv$", recursive = TRUE)) == 0) {
  write.csv(data.frame(message = "No WormCat results"),
            file.path(args$outdir, "wormcat_no_results.csv"), row.names = FALSE)
}

message("Done. Output: ", args$outdir)
