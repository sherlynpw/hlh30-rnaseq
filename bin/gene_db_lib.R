# gene_db_lib.R - Helper functions to query the gene SQLite database
#
# Tables in gene_db.sqlite (built by build_gene_db.R from go-kegg cache):
#   wb_to_entrez:   wormbase_gene, entrez_id
#   entrez_to_wb:   entrez_id, wormbase_gene
#   annotations:    wormbase_gene, gene_name, gene_class_description,
#                   concise_description, legacy_description, sequence_name,
#                   gene_type, operon, other_names, status

suppressPackageStartupMessages(library(RSQLite))

open_gene_db <- function(path) {
  if (!file.exists(path)) stop("Gene database not found: ", path)
  dbConnect(SQLite(), path, flags = SQLITE_RO)
}

close_gene_db <- function(db) {
  dbDisconnect(db)
}

# WormBase IDs -> Entrez IDs
wb_to_entrez <- function(db, wb_ids) {
  placeholders <- paste(rep("?", length(wb_ids)), collapse = ",")
  query <- paste0("SELECT wormbase_gene AS WORMBASE, entrez_id AS ENTREZID FROM wb_to_entrez WHERE wormbase_gene IN (", placeholders, ")")
  dbGetQuery(db, query, params = as.list(wb_ids))
}

# Entrez IDs -> WormBase IDs
entrez_to_wb <- function(db, entrez_ids) {
  placeholders <- paste(rep("?", length(entrez_ids)), collapse = ",")
  query <- paste0("SELECT entrez_id AS ENTREZID, wormbase_gene AS WORMBASE FROM entrez_to_wb WHERE entrez_id IN (", placeholders, ")")
  dbGetQuery(db, query, params = as.list(as.character(entrez_ids)))
}

# Get annotations for WormBase IDs
get_annotations <- function(db, wb_ids) {
  placeholders <- paste(rep("?", length(wb_ids)), collapse = ",")
  query <- paste0("SELECT * FROM annotations WHERE wormbase_gene IN (", placeholders, ")")
  dbGetQuery(db, query, params = as.list(wb_ids))
}

# --- Annotation utilities ---

# Open gene_db if path is provided, return NULL otherwise
maybe_open_gene_db <- function(path) {
  if (!is.null(path) && nchar(path) > 0 && file.exists(path)) {
    message("Using gene_db for annotation: ", path)
    open_gene_db(path)
  } else {
    NULL
  }
}

# Annotate a data.frame: add gene_name and description columns
# df must have a column with WormBase IDs (default: "gene_id")
# Safe to call with db=NULL (returns df unchanged)
annotate_df <- function(db, df, id_col = "gene_id") {
  if (is.null(db) || nrow(df) == 0 || !id_col %in% colnames(df)) return(df)
  wb_ids <- unique(df[[id_col]])
  annos <- get_annotations(db, wb_ids)
  if (nrow(annos) == 0) return(df)
  annos_slim <- annos[, c("wormbase_gene", "gene_name", "concise_description"), drop = FALSE]
  annos_slim <- annos_slim[!duplicated(annos_slim$wormbase_gene), ]
  merged <- merge(df, annos_slim, by.x = id_col, by.y = "wormbase_gene", all.x = TRUE, sort = FALSE)
  merged$gene_name[is.na(merged$gene_name)] <- merged[[id_col]][is.na(merged$gene_name)]
  merged
}

# Annotate by Entrez ID: look up WormBase ID first, then get annotations
annotate_entrez_df <- function(db, df, id_col = "ENTREZID") {
  if (is.null(db) || nrow(df) == 0 || !id_col %in% colnames(df)) return(df)
  entrez_ids <- unique(df[[id_col]])
  wb_map <- entrez_to_wb(db, entrez_ids)
  if (nrow(wb_map) == 0) return(df)
  annos <- get_annotations(db, wb_map$WORMBASE)
  if (nrow(annos) == 0) return(df)
  annos_slim <- annos[, c("wormbase_gene", "gene_name", "concise_description"), drop = FALSE]
  annos_slim <- annos_slim[!duplicated(annos_slim$wormbase_gene), ]
  wb_map <- merge(wb_map, annos_slim, by.x = "WORMBASE", by.y = "wormbase_gene", all.x = TRUE)
  merged <- merge(df, wb_map[, c("ENTREZID", "WORMBASE", "gene_name", "concise_description")],
                  by.x = id_col, by.y = "ENTREZID", all.x = TRUE, sort = FALSE)
  merged
}

# Write an annotated gene list: gene_id -> gene_id, gene_name, description
annotate_gene_list <- function(db, gene_ids, output_file) {
  df <- data.frame(gene_id = gene_ids, stringsAsFactors = FALSE)
  df <- annotate_df(db, df)
  write.csv(df, output_file, row.names = FALSE)
}

# WormBase IDs -> KEGG IDs (CELE_ + sequence_name)
wb_to_kegg <- function(db, wb_ids) {
  placeholders <- paste(rep("?", length(wb_ids)), collapse = ",")
  query <- paste0("SELECT wormbase_gene, sequence_name FROM annotations WHERE wormbase_gene IN (", placeholders, ") AND sequence_name IS NOT NULL AND sequence_name != ''")
  result <- dbGetQuery(db, query, params = as.list(wb_ids))
  if (nrow(result) > 0) {
    result$kegg_id <- paste0("CELE_", result$sequence_name)
  }
  result
}

# --- Filtering helpers ---

# Extract gene lists from DESeq2 results CSVs with configurable threshold
# Returns a named list: label -> list(de=genes, up=genes, down=genes)
extract_gene_lists_from_results <- function(results_dir, filter_by = "padj", threshold = 0.05, lfc = 0) {
  result_files <- list.files(results_dir, pattern = "\\.csv$", full.names = TRUE)
  all_lists <- list()

  for (f in result_files) {
    comparison <- tools::file_path_sans_ext(basename(f))
    res <- read.csv(f, stringsAsFactors = FALSE)

    # Filter by specified column
    if (filter_by == "pvalue") {
      de <- res[!is.na(res$pvalue) & res$pvalue < threshold, ]
    } else {
      de <- res[!is.na(res$padj) & res$padj < threshold, ]
    }

    # Additional LFC filter if specified
    if (lfc > 0) {
      de <- de[abs(de$log2FoldChange) >= lfc, ]
    }

    if (nrow(de) == 0) next

    # Sort by stat
    if ("stat" %in% colnames(de)) de <- de[order(-de$stat), ]

    de_up   <- de[de$log2FoldChange > 0, ]
    de_down <- de[de$log2FoldChange < 0, ]

    all_lists[[paste0(comparison, "_genes_de")]] <- de$gene_id
    if (nrow(de_up) > 0)   all_lists[[paste0(comparison, "_genes_up")]] <- de_up$gene_id
    if (nrow(de_down) > 0) all_lists[[paste0(comparison, "_genes_down")]] <- de_down$gene_id
  }

  all_lists
}

# Find gene lists from either: DESeq2 results CSVs (filtered) or pre-made .txt files
# Returns named list: label -> gene_ids
find_gene_lists <- function(input_dir, filter_by = NULL, threshold = 0.05, lfc = 0) {
  # Check if directory has CSVs (DESeq2 results) or TXT files (gene lists)
  csv_files <- list.files(input_dir, pattern = "\\.csv$", recursive = TRUE)
  txt_files <- list.files(input_dir, pattern = "genes.*\\.txt$", recursive = TRUE, full.names = TRUE)

  if (!is.null(filter_by) && length(csv_files) > 0) {
    # Mode 1: Filter DESeq2 results CSVs
    message("Filtering DESeq2 results by ", filter_by, " < ", threshold,
            if (lfc > 0) paste0(", |lfc| >= ", lfc) else "")
    return(extract_gene_lists_from_results(input_dir, filter_by, threshold, lfc))
  }

  if (length(txt_files) > 0) {
    # Mode 2: Read pre-made gene list .txt files
    message("Reading ", length(txt_files), " gene list files")
    all_lists <- list()
    for (f in txt_files) {
      rel_path <- sub(paste0("^", normalizePath(input_dir), "/?"), "", normalizePath(f))
      label <- tools::file_path_sans_ext(gsub("/", "_", rel_path))
      genes <- readLines(f)
      genes <- genes[nchar(genes) > 0]
      if (length(genes) > 0) all_lists[[label]] <- genes
    }
    return(all_lists)
  }

  message("No CSV or TXT gene list files found in ", input_dir)
  list()
}

# Get all WormBase -> Entrez mappings (for building gene sets)
get_all_mappings <- function(db) {
  dbGetQuery(db, "SELECT wormbase_gene AS WORMBASE, entrez_id AS ENTREZID FROM wb_to_entrez WHERE entrez_id IS NOT NULL AND entrez_id != ''")
}
