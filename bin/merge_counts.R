#!/usr/bin/env Rscript

# merge_counts.R - Merge per-sample featureCounts into a single count matrix
#
# Input:  Multiple featureCounts .txt files + samplesheet YAML
# Output: gene_counts.tsv (genes x samples) + sample_meta.tsv

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

opts <- list(
  make_option(c("--counts_dir"), type = "character", default = ".", help = "Directory containing featureCounts files"),
  make_option(c("--pattern"),    type = "character", default = "*.featureCounts.txt", help = "Glob pattern for count files"),
  make_option(c("--metadata"),   type = "character", help = "Samplesheet YAML file"),
  make_option(c("--outdir"),     type = "character", default = ".", help = "Output directory")
)

args <- parse_args(OptionParser(option_list = opts))

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# Find count files
count_files <- list.files(args$counts_dir, pattern = glob2rx(args$pattern), full.names = TRUE)
if (length(count_files) == 0) {
  # Try current directory
  count_files <- list.files(".", pattern = glob2rx(args$pattern), full.names = TRUE)
}
stopifnot("No count files found" = length(count_files) > 0)

# Read featureCounts files and merge
# featureCounts format: Geneid, Chr, Start, End, Strand, Length, Count
read_featurecounts <- function(f) {
  df <- read.delim(f, comment.char = "#", stringsAsFactors = FALSE)
  sample_name <- gsub("\\.featureCounts\\.txt$", "", basename(f))
  # Last column is the count
  counts <- data.frame(
    gene_id = df$Geneid,
    count   = df[[ncol(df)]],
    stringsAsFactors = FALSE
  )
  colnames(counts)[2] <- sample_name
  counts
}

message("Reading ", length(count_files), " count files...")
count_list <- lapply(count_files, read_featurecounts)

# Merge all into one matrix
count_matrix <- count_list[[1]]
if (length(count_list) > 1) {
  for (i in 2:length(count_list)) {
    count_matrix <- merge(count_matrix, count_list[[i]], by = "gene_id", all = TRUE)
  }
}

# Replace NAs with 0
count_matrix[is.na(count_matrix)] <- 0

# Set gene_id as rownames
rownames(count_matrix) <- count_matrix$gene_id
count_matrix$gene_id <- NULL

# Ensure integer counts
count_matrix[] <- lapply(count_matrix, as.integer)

message("Count matrix: ", nrow(count_matrix), " genes x ", ncol(count_matrix), " samples")

# Write count matrix
write.table(
  count_matrix,
  file      = file.path(args$outdir, "gene_counts.tsv"),
  sep       = "\t",
  quote     = FALSE,
  col.names = NA  # includes rownames header
)

# Parse samplesheet and write flat metadata TSV
if (!is.null(args$metadata) && file.exists(args$metadata)) {
  samplesheet <- yaml::read_yaml(args$metadata)
  # Dynamically extract all metadata columns (skip file path fields)
  skip_keys <- c("bam", "fastq_1", "fastq_2")
  meta_rows <- lapply(samplesheet$samples, function(s) {
    fields <- s[!names(s) %in% skip_keys]
    # Apply defaults
    if (is.null(fields$source))       fields$source       <- "unknown"
    if (is.null(fields$strandedness)) fields$strandedness  <- "auto"
    as.data.frame(fields, stringsAsFactors = FALSE)
  })
  meta_df <- do.call(rbind, meta_rows)

  # Only keep samples that are in the count matrix
  meta_df <- meta_df[meta_df$sample_id %in% colnames(count_matrix), ]

  write.table(
    meta_df,
    file      = file.path(args$outdir, "sample_meta.tsv"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  message("Sample metadata: ", nrow(meta_df), " samples")
}

message("Done. Output: ", args$outdir)
