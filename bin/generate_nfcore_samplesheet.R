#!/usr/bin/env Rscript

# generate_nfcore_samplesheet.R - Convert our YAML samplesheet to nf-core CSV format
#
# Called by run.sh to generate nf-core/rnaseq compatible samplesheets
# Splits FASTQ samples into separate samplesheets for 3' and standard RNA-seq

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

opts <- list(
  make_option("--input",  type = "character", help = "Pipeline samplesheet YAML"),
  make_option("--outdir", type = "character", default = ".", help = "Output directory for CSV samplesheets")
)

args <- parse_args(OptionParser(option_list = opts))

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

samplesheet <- yaml::read_yaml(args$input)

# Split by input type
samples_3prime   <- list()
samples_standard <- list()

for (s in samplesheet$samples) {
  if (s$input_type == "fastq_3prime") {
    samples_3prime <- c(samples_3prime, list(s))
  } else if (s$input_type == "fastq_standard") {
    samples_standard <- c(samples_standard, list(s))
  }
  # plasmidsaurus_bam samples are handled directly by main.nf, not nf-core
}

# Write nf-core CSV samplesheet
write_nfcore_csv <- function(samples, output_file) {
  if (length(samples) == 0) return(invisible(NULL))

  lines <- c("sample,fastq_1,fastq_2,strandedness")
  for (s in samples) {
    fq1    <- s$fastq_1
    fq2    <- ifelse(is.null(s$fastq_2), "", s$fastq_2)
    strand <- ifelse(is.null(s$strandedness), "auto", s$strandedness)
    lines  <- c(lines, paste(s$sample_id, fq1, fq2, strand, sep = ","))
  }

  writeLines(lines, output_file)
  message("Wrote ", length(samples), " samples to ", output_file)
}

write_nfcore_csv(samples_3prime,   file.path(args$outdir, "samplesheet_3prime.csv"))
write_nfcore_csv(samples_standard, file.path(args$outdir, "samplesheet_standard.csv"))

# Report
n_plasm <- sum(sapply(samplesheet$samples, function(s) s$input_type == "plasmidsaurus_bam"))
message("Summary:")
message("  Plasmidsaurus BAMs (skip alignment): ", n_plasm)
message("  3' RNA-seq FASTQ:  ", length(samples_3prime))
message("  Standard FASTQ:    ", length(samples_standard))
