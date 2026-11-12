#!/usr/bin/env Rscript

# build_gene_db.R - Build a SQLite database from go-kegg cache files
#
# Imports existing cache (47K files) into a single portable SQLite database.
#
# Tables:
#   - gene_mapping:   WBGene <-> Entrez ID (from org.Ce.eg.db)
#   - annotations:    WBGene -> gene name, descriptions, type, operon (from cache)
#   - wb_to_entrez:   WBGene -> Entrez (from cache)
#   - entrez_to_wb:   Entrez -> WBGene (from cache)
#
# Usage:
#   Rscript bin/build_gene_db.R --cache /path/to/go-kegg/.cache --output gene_db.sqlite

suppressPackageStartupMessages({
  library(optparse)
  library(RSQLite)
  library(jsonlite)
})

opts <- list(
  make_option("--cache",  type = "character", help = "Path to go-kegg .cache directory"),
  make_option("--output", type = "character", default = "gene_db.sqlite", help = "Output SQLite file")
)

args <- parse_args(OptionParser(option_list = opts))

if (is.null(args$cache) || !dir.exists(args$cache)) {
  stop("Please provide --cache /path/to/go-kegg/.cache")
}

if (file.exists(args$output)) file.remove(args$output)
db <- dbConnect(SQLite(), args$output)

# --- Table 1: Annotations (from cache/annotations/*.json) ---
anno_dir <- file.path(args$cache, "annotations")
message("Importing annotations from ", anno_dir, "...")
json_files <- list.files(anno_dir, pattern = "\\.json$", full.names = TRUE)

cached_annos <- lapply(json_files, function(f) {
  tryCatch({
    j <- fromJSON(f)
    data.frame(
      wormbase_gene          = gsub("\\.json$", "", basename(f)),
      gene_name              = j$gene_name %||% NA_character_,
      gene_class_description = j$gene_class_description %||% NA_character_,
      concise_description    = j$concise_description %||% NA_character_,
      legacy_description     = j$legacy_description %||% NA_character_,
      sequence_name          = j$sequence_name %||% NA_character_,
      gene_type              = j$gene_type %||% NA_character_,
      operon                 = j$operon %||% NA_character_,
      other_names            = j$other_names %||% NA_character_,
      status                 = j$status %||% NA_character_,
      stringsAsFactors       = FALSE
    )
  }, error = function(e) NULL)
})
anno_df <- do.call(rbind, Filter(Negate(is.null), cached_annos))
dbWriteTable(db, "annotations", anno_df, overwrite = TRUE)
message("  ", nrow(anno_df), " annotations")

# --- Table 2: WormBase -> Entrez (from cache/entrez/*.txt) ---
entrez_dir <- file.path(args$cache, "entrez")
message("Importing WB->Entrez from ", entrez_dir, "...")
entrez_files <- list.files(entrez_dir, pattern = "\\.txt$", full.names = TRUE)

wb2e <- lapply(entrez_files, function(f) {
  tryCatch({
    entrez_id <- trimws(readLines(f, n = 1, warn = FALSE))
    data.frame(
      wormbase_gene = gsub("\\.txt$", "", basename(f)),
      entrez_id     = entrez_id,
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
})
wb2e_df <- do.call(rbind, Filter(Negate(is.null), wb2e))
dbWriteTable(db, "wb_to_entrez", wb2e_df, overwrite = TRUE)
message("  ", nrow(wb2e_df), " mappings")

# --- Table 3: Entrez -> WormBase (from cache/entrez-to-wormbase/*.txt) ---
rev_dir <- file.path(args$cache, "entrez-to-wormbase")
message("Importing Entrez->WB from ", rev_dir, "...")
rev_files <- list.files(rev_dir, pattern = "\\.txt$", full.names = TRUE)

e2wb <- lapply(rev_files, function(f) {
  tryCatch({
    wb_id <- trimws(readLines(f, n = 1, warn = FALSE))
    data.frame(
      entrez_id     = gsub("\\.txt$", "", basename(f)),
      wormbase_gene = wb_id,
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
})
e2wb_df <- do.call(rbind, Filter(Negate(is.null), e2wb))
dbWriteTable(db, "entrez_to_wb", e2wb_df, overwrite = TRUE)
message("  ", nrow(e2wb_df), " mappings")

# --- Indices ---
message("Creating indices...")
dbExecute(db, "CREATE INDEX idx_anno_wb ON annotations(wormbase_gene)")
dbExecute(db, "CREATE INDEX idx_wb2e_wb ON wb_to_entrez(wormbase_gene)")
dbExecute(db, "CREATE INDEX idx_wb2e_e ON wb_to_entrez(entrez_id)")
dbExecute(db, "CREATE INDEX idx_e2wb_e ON entrez_to_wb(entrez_id)")
dbExecute(db, "CREATE INDEX idx_e2wb_wb ON entrez_to_wb(wormbase_gene)")

# --- Summary ---
for (t in dbListTables(db)) {
  n <- dbGetQuery(db, paste("SELECT COUNT(*) as n FROM", t))$n
  message("  ", t, ": ", n, " rows")
}

db_size <- file.info(args$output)$size / 1024 / 1024
message("")
message("Done: ", args$output, " (", round(db_size, 1), " MB)")

dbDisconnect(db)
