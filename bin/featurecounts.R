#!/usr/bin/env Rscript

# featurecounts.R - Run Rsubread featureCounts on a BAM file

suppressPackageStartupMessages(library(Rsubread))

args <- commandArgs(trailingOnly = TRUE)
# args: bam gtf output_prefix strand_specific is_paired_end nthreads

bam           <- args[1]
gtf           <- args[2]
output_prefix <- args[3]
strand        <- as.integer(args[4])
paired        <- as.logical(args[5])
nthreads      <- as.integer(args[6])

fc <- featureCounts(
  files               = bam,
  annot.ext           = gtf,
  isGTFAnnotationFile = TRUE,
  GTF.featureType     = "exon",
  GTF.attrType        = "gene_id",
  strandSpecific      = strand,
  isPairedEnd         = paired,
  nthreads            = nthreads
)

# Write counts
counts_df <- data.frame(
  Geneid = rownames(fc$counts),
  fc$annotation,
  Count  = fc$counts[, 1]
)
write.table(counts_df, paste0(output_prefix, ".featureCounts.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Write summary
write.table(fc$stat, paste0(output_prefix, ".featureCounts.txt.summary"),
            sep = "\t", row.names = FALSE, quote = FALSE)

message("Done: ", output_prefix)
