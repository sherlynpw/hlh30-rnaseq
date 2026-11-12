#!/usr/bin/env Rscript

# gsea.R - GSEA per comparison using fgsea
#
# Collections come from two sources:
#   --builtin             comma list of built-ins: go, kegg
#   --collections_manifest path to a JSON file describing external collections:
#       [
#         {"name": "hallmark", "file": "hallmark.gmt", "format": "gmt", "map_ortholog": true},
#         {"name": "wormcat_cat1", "file": "wormcat_cat1.gmt", "map_ortholog": false}
#       ]
#   The "file" field is resolved relative to --collections_dir.

suppressPackageStartupMessages({
  library(optparse)
  library(fgsea)
  library(clusterProfiler)
  library(org.Ce.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(dplyr)
  library(RSQLite)
  library(BiocParallel)
  library(jsonlite)
})

.find_lib <- function() {
  for (p in c(Sys.which("gene_db_lib.R"), file.path("bin", "gene_db_lib.R")))
    if (nchar(p) > 0 && file.exists(p)) return(p)
  NULL
}
lib <- .find_lib()
if (!is.null(lib)) source(lib)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || !nzchar(x)) y else x

opts <- list(
  make_option("--de_results_dir", type = "character", default = ".", help = "Directory with DESeq2 result CSVs"),
  make_option("--gene_db",    type = "character", default = NULL, help = "Path to gene_db.sqlite"),
  make_option("--ranking",    type = "character", default = "stat", help = "Ranking column"),
  make_option("--org_db",     type = "character", default = "org.Ce.eg.db", help = "OrgDb package"),
  make_option("--kegg_org",   type = "character", default = "cel", help = "KEGG organism code"),
  make_option("--builtin",    type = "character", default = "",
              help = "Comma-separated built-in collections: go, kegg"),
  make_option("--collections_manifest", type = "character", default = NULL,
              help = "Path to JSON manifest of external GMT collections"),
  make_option("--collections_dir", type = "character", default = ".",
              help = "Directory containing collection files referenced by manifest"),
  make_option("--ortholog_map", type = "character", default = NULL,
              help = "Optional TSV/CSV with source->target ortholog mapping"),
  make_option("--ortholog_from", type = "character", default = "wormbase_gene",
              help = "Column in ortholog_map matching DE result gene IDs"),
  make_option("--ortholog_to", type = "character", default = "human_symbol",
              help = "Column in ortholog_map matching external gene-set IDs"),
  make_option("--min_size",   type = "integer",   default = 15, help = "Min gene set size"),
  make_option("--max_size",   type = "integer",   default = 500, help = "Max gene set size"),
  make_option("--top_n",      type = "integer",   default = 20, help = "Top N pathways for plots"),
  make_option("--outdir",     type = "character",  default = ".", help = "Output directory")
)

args <- parse_args(OptionParser(option_list = opts))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

org_db <- get(args$org_db)

split_csv <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character())
  vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

is_true <- function(x) {
  if (is.logical(x)) return(isTRUE(x))
  tolower(trimws(as.character(x %||% "false"))) %in% c("true", "t", "1", "yes", "y")
}

read_gmt_sets <- function(path) {
  if (is.null(path) || !file.exists(path)) stop("GMT file not found: ", path)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  sets <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    genes <- unique(fields[-c(1, 2)])
    genes[nzchar(genes)]
  })
  names(sets) <- vapply(strsplit(lines, "\t", fixed = TRUE), `[`, character(1), 1)
  sets[lengths(sets) > 0]
}

read_gene_list_sets <- function(path, set_name) {
  if (is.null(path) || !file.exists(path)) stop("Gene list file not found: ", path)
  genes <- unique(trimws(readLines(path, warn = FALSE)))
  genes <- genes[nzchar(genes)]
  if (length(genes) == 0) return(list())
  out <- list(genes)
  names(out) <- set_name
  out
}

read_mapping_table <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read.csv(path, stringsAsFactors = FALSE)
  } else {
    read.delim(path, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
  }
}

map_external_sets_to_source_ids <- function(pathways, mapping_df, source_col, target_col) {
  if (is.null(mapping_df)) return(pathways)
  if (!all(c(source_col, target_col) %in% colnames(mapping_df))) {
    stop("Ortholog map must contain columns: ", source_col, ", ", target_col)
  }

  mapping_df <- mapping_df[, c(source_col, target_col)]
  mapping_df <- mapping_df[!is.na(mapping_df[[source_col]]) & !is.na(mapping_df[[target_col]]), , drop = FALSE]
  mapping_df[[source_col]] <- trimws(as.character(mapping_df[[source_col]]))
  mapping_df[[target_col]] <- trimws(as.character(mapping_df[[target_col]]))
  mapping_df <- mapping_df[nzchar(mapping_df[[source_col]]) & nzchar(mapping_df[[target_col]]), , drop = FALSE]

  target_to_source <- split(mapping_df[[source_col]], mapping_df[[target_col]])

  mapped <- lapply(pathways, function(genes) {
    unique(unlist(target_to_source[genes], use.names = FALSE))
  })
  mapped[lengths(mapped) > 0]
}

build_builtin_sets <- function(kind, org_db, wb_entrez_map, kegg_org) {
  if (kind == "go") {
    message("Building GO gene sets...")
    return(tryCatch({
      go2gene <- AnnotationDbi::select(org_db, keys = keys(org_db, "WORMBASE"),
                                       keytype = "WORMBASE", columns = "GO")
      go2gene <- go2gene[!is.na(go2gene$GO), c("GO", "WORMBASE")]
      split(go2gene$WORMBASE, go2gene$GO)
    }, error = function(e) {
      message("Warning: ", e$message)
      list()
    }))
  }

  if (kind == "kegg") {
    message("Building KEGG gene sets...")
    return(tryCatch({
      kegg_data <- download_KEGG(kegg_org)
      if (length(wb_entrez_map) == 0) stop("No WB-Entrez mappings")
      entrez_to_wb_map <- setNames(names(wb_entrez_map), wb_entrez_map)
      kegg_merged <- kegg_data$KEGGPATHID2EXTID
      kegg_merged$wb <- entrez_to_wb_map[kegg_merged$to]
      kegg_merged <- kegg_merged[!is.na(kegg_merged$wb), ]
      split(kegg_merged$wb, kegg_merged$from)
    }, error = function(e) {
      message("Warning: ", e$message)
      list()
    }))
  }

  stop("Unknown builtin collection: ", kind)
}

build_external_sets <- function(name, path, format, mapping_df, source_col, target_col) {
  if (is.null(path) || !file.exists(path)) {
    stop("Missing file for collection '", name, "': ", path)
  }
  format <- tolower(format %||% "gmt")
  if (!format %in% c("gmt", "list")) {
    stop("Unsupported gene-set format for '", name, "': ", format)
  }

  raw_sets <- if (format == "list") {
    read_gene_list_sets(path, name)
  } else {
    read_gmt_sets(path)
  }

  if (!is.null(mapping_df)) {
    raw_sets <- map_external_sets_to_source_ids(raw_sets, mapping_df, source_col, target_col)
  }

  raw_sets[lengths(raw_sets) > 0]
}

# Build WormBase -> Entrez mapping
if (!is.null(args$gene_db) && file.exists(args$gene_db)) {
  gdb <- open_gene_db(args$gene_db)
  all_mappings <- get_all_mappings(gdb)
  close_gene_db(gdb)
  wb_entrez_map <- setNames(all_mappings$ENTREZID, all_mappings$WORMBASE)
} else {
  mapping <- AnnotationDbi::select(org_db, keys = keys(org_db, "WORMBASE"),
                                    keytype = "WORMBASE", columns = "ENTREZID")
  mapping <- mapping[!is.na(mapping$ENTREZID), ]
  wb_entrez_map <- setNames(mapping$ENTREZID, mapping$WORMBASE)
}

mapping_df <- read_mapping_table(args$ortholog_map)
gene_set_collections <- list()

# Built-ins
for (kind in split_csv(args$builtin)) {
  if (!kind %in% c("go", "kegg")) {
    message("Warning: unknown builtin '", kind, "', skipping")
    next
  }
  gene_set_collections[[kind]] <- build_builtin_sets(kind, org_db, wb_entrez_map, args$kegg_org)
}

# External collections from manifest
if (!is.null(args$collections_manifest) && nzchar(args$collections_manifest) &&
    file.exists(args$collections_manifest)) {
  manifest <- jsonlite::fromJSON(args$collections_manifest, simplifyDataFrame = FALSE)
  for (entry in manifest) {
    name <- entry$name
    if (is.null(name) || !nzchar(name)) {
      message("Warning: manifest entry missing 'name', skipping")
      next
    }
    if (name %in% names(gene_set_collections)) {
      message("Warning: collection '", name, "' already defined (builtin?), overwriting")
    }
    file_path <- file.path(args$collections_dir, entry$file)
    fmt <- entry$format %||% "gmt"
    use_map <- if (is_true(entry$map_ortholog)) mapping_df else NULL
    gene_set_collections[[name]] <- tryCatch({
      build_external_sets(name, file_path, fmt, use_map, args$ortholog_from, args$ortholog_to)
    }, error = function(e) {
      message("Warning [", name, "]: ", e$message)
      list()
    })
  }
}

for (nm in names(gene_set_collections)) {
  message("Collection ", nm, ": ", length(gene_set_collections[[nm]]), " pathways")
}

# Find result files
result_files <- list.files(args$de_results_dir, pattern = "\\.csv$", full.names = TRUE)

for (f in result_files) {
  comparison <- tools::file_path_sans_ext(basename(f))
  message("GSEA: ", comparison)

  res <- read.csv(f, stringsAsFactors = FALSE)

  out_sub <- file.path(args$outdir, comparison)
  dir.create(out_sub, recursive = TRUE, showWarnings = FALSE)

  if (args$ranking == "pvalue_lfc") {
    res$rank_val <- -log10(pmax(res$pvalue, 1e-300)) * sign(res$log2FoldChange)
  } else {
    res$rank_val <- res[[args$ranking]]
  }

  if (is.null(res$rank_val) || all(is.na(res$rank_val))) {
    message("  Ranking column '", args$ranking, "' not found or all NA, skipping")
    write.csv(data.frame(message = paste("Ranking column", args$ranking, "not found or all NA")),
              file.path(out_sub, "gsea_skipped.csv"), row.names = FALSE)
    next
  }
  res <- res[!is.na(res$rank_val), ]
  ranked <- setNames(res$rank_val, res$gene_id)
  ranked <- sort(ranked, decreasing = TRUE)

  for (collection in names(gene_set_collections)) {
    collection_sets <- gene_set_collections[[collection]]
    out_collection <- file.path(out_sub, collection)
    dir.create(out_collection, recursive = TRUE, showWarnings = FALSE)

    if (length(collection_sets) == 0) {
      write.csv(data.frame(message = paste("No pathways available for", collection)),
                file.path(out_collection, "gsea_skipped.csv"), row.names = FALSE)
      next
    }

    tryCatch({
      fgsea_res <- fgseaMultilevel(
        pathways = collection_sets,
        stats = ranked,
        minSize = args$min_size,
        maxSize = args$max_size,
        nproc = 1,
        BPPARAM = SerialParam(progressbar = FALSE)
      )

      if (nrow(fgsea_res) > 0) {
        fgsea_res$comparison <- comparison
        fgsea_res$collection <- collection
        fgsea_res$leadingEdge <- sapply(fgsea_res$leadingEdge, paste, collapse = ";")
        fgsea_df <- as.data.frame(fgsea_res)

        write.csv(fgsea_df, file.path(out_collection, "gsea_results.csv"), row.names = FALSE)

        sig <- fgsea_df %>% filter(padj < 0.25) %>% arrange(padj) %>% head(args$top_n)
        if (nrow(sig) > 0) {
          p <- ggplot(sig, aes(x = NES, y = reorder(pathway, NES), fill = padj)) +
            geom_col() + scale_fill_gradient(low = "red", high = "blue") +
            labs(x = "Normalized Enrichment Score", y = NULL,
                 title = paste("GSEA:", comparison, "(", collection, ")"), fill = "padj") +
            theme_minimal(base_size = 9) +
            geom_vline(xintercept = 0, linetype = "dashed")
          ggsave(file.path(out_collection, "gsea_barplot.pdf"),
                 plot = p, width = 10, height = max(3, nrow(sig) * 0.25))
        }

        n_sig <- sum(fgsea_df$padj < 0.25, na.rm = TRUE)
        message("  [", collection, "] ", nrow(fgsea_df), " pathways, ", n_sig,
                " significant (padj<0.25)")
      } else {
        write.csv(data.frame(message = "No results"),
                  file.path(out_collection, "gsea_results.csv"), row.names = FALSE)
      }
    }, error = function(e) {
      message("  [", collection, "] Error: ", e$message)
      write.csv(data.frame(message = e$message),
                file.path(out_collection, "gsea_failed.csv"), row.names = FALSE)
    })
  }
}

message("Done. Output: ", args$outdir)
