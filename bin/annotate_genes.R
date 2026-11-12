#!/usr/bin/env Rscript

# annotate_genes.R - Gene annotation from gene_db.sqlite (no BioMart calls)

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(RSQLite)
})

.find_lib <- function() {
  for (p in c(Sys.which("gene_db_lib.R"), file.path("bin", "gene_db_lib.R")))
    if (nchar(p) > 0 && file.exists(p)) return(p)
  stop("Cannot find gene_db_lib.R")
}
source(.find_lib())

opts <- list(
  make_option("--genes",   type = "character", help = "Gene list file or directory"),
  make_option("--gene_db", type = "character", default = NULL, help = "Path to gene_db.sqlite"),
  make_option("--dataset", type = "character", default = "celegans_gene_ensembl", help = "BioMart dataset (fallback)"),
  make_option("--outdir",  type = "character", default = ".", help = "Output directory")
)

args <- parse_args(OptionParser(option_list = opts))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# Find gene files
genes_base <- args$genes
if (file.info(args$genes)$isdir) {
  gene_files <- list.files(args$genes, pattern = "genes.*\\.txt$", full.names = TRUE, recursive = TRUE)
} else {
  gene_files <- args$genes
  genes_base <- dirname(args$genes)
}

# Collect all unique genes
all_genes <- unique(unlist(lapply(gene_files, function(f) {
  g <- readLines(f); g[nchar(g) > 0]
})))

message("Annotating ", length(all_genes), " unique genes...")

# Get annotations from SQLite or BioMart fallback
if (!is.null(args$gene_db) && file.exists(args$gene_db)) {
  gdb <- open_gene_db(args$gene_db)
  annotations <- get_annotations(gdb, all_genes)
  close_gene_db(gdb)
  message("  Using gene_db.sqlite")
} else {
  message("  No gene_db, falling back to BioMart...")
  annotations <- tryCatch({
    suppressPackageStartupMessages(library(biomaRt))
    mart <- tryCatch(
      useEnsembl(biomart = "ensembl", dataset = args$dataset),
      error = function(e) useEnsembl(biomart = "ensembl", dataset = args$dataset, mirror = "uswest")
    )
    ann <- getBM(
      attributes = c("wormbase_gene", "external_gene_name", "gene_biotype", "description"),
      filters = "wormbase_gene", values = all_genes, mart = mart
    )
    ann$description <- gsub("\\s*\\[Source:.*\\]$", "", ann$description)
    ann[!duplicated(ann$wormbase_gene), ]
  }, error = function(e) {
    message("  BioMart failed: ", e$message)
    message("  Writing gene IDs without annotations")
    data.frame(wormbase_gene = all_genes, external_gene_name = all_genes,
               gene_biotype = "unknown", description = all_genes, stringsAsFactors = FALSE)
  })
}

# Write combined annotations
write.csv(annotations, file.path(args$outdir, "annotated_genes.csv"), row.names = FALSE)

# Write per-file annotated versions (recursive-aware: preserve path in output name)
for (f in gene_files) {
  rel_path <- sub(paste0("^", normalizePath(genes_base), "/?"), "", normalizePath(f))
  label <- tools::file_path_sans_ext(gsub("/", "_", rel_path))

  genes <- readLines(f)
  genes <- genes[nchar(genes) > 0]

  annotated <- data.frame(wormbase_gene = genes, stringsAsFactors = FALSE) %>%
    left_join(annotations, by = "wormbase_gene")

  write.csv(annotated, file.path(args$outdir, paste0(label, "_annotated.csv")), row.names = FALSE)
}

message("Done. Annotated ", nrow(annotations), "/", length(all_genes), " genes in ", args$outdir)
