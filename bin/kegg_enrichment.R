#!/usr/bin/env Rscript

# kegg_enrichment.R - KEGG pathway enrichment
# Accepts either DESeq2 result CSVs (filters internally) or gene list .txt files

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(org.Ce.eg.db)
  library(ggplot2)
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
  make_option("--gene_db",   type = "character",  default = NULL, help = "Path to gene_db.sqlite"),
  make_option("--universe",  type = "character",  default = NULL, help = "Background gene list (WormBase IDs) = enrichment universe"),
  make_option("--kegg_org",  type = "character",  default = "cel", help = "KEGG organism code"),
  make_option("--pvalue",    type = "numeric",   default = 0.05, help = "enricher p-value cutoff"),
  make_option("--qvalue",    type = "numeric",   default = 0.2, help = "enricher q-value cutoff"),
  make_option("--top_n",     type = "integer",   default = 20, help = "Top N for plots"),
  make_option("--cpus",      type = "integer",   default = 1, help = "Parallel workers (mclapply)"),
  make_option("--outdir",    type = "character",  default = ".", help = "Output directory")
)

args <- parse_args(OptionParser(option_list = opts))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# ID conversion (WormBase -> KEGG CELE_ format)
if (!is.null(args$gene_db) && file.exists(args$gene_db)) {
  gdb <- open_gene_db(args$gene_db)
  convert_ids <- function(wb_ids) {
    m <- wb_to_kegg(gdb, wb_ids)
    m[!is.na(m$kegg_id), ]
  }
} else {
  convert_ids <- function(wb_ids) {
    m <- AnnotationDbi::select(org.Ce.eg.db, keys = wb_ids, keytype = "WORMBASE", columns = "ENTREZID")
    m[!is.na(m$ENTREZID), ]
  }
}

# Download KEGG data once
message("Downloading KEGG pathway data for ", args$kegg_org, "...")
kegg_data <- tryCatch(download_KEGG(args$kegg_org), error = function(e) {
  message("KEGG download failed: ", e$message); NULL
})

# Background universe (genes tested in this analysis)
universe_kegg <- NULL; universe_entrez <- NULL
if (!is.null(args$universe) && nzchar(args$universe) && file.exists(args$universe)) {
  ug <- readLines(args$universe)
  ug <- ug[!is.na(ug) & nzchar(ug)]
  if (length(ug) > 0) {
    um <- convert_ids(ug)
    if ("kegg_id"  %in% colnames(um)) universe_kegg   <- unique(um$kegg_id)
    if ("ENTREZID" %in% colnames(um)) universe_entrez <- unique(um$ENTREZID)
    message("Universe: ", nrow(um), " background genes")
  }
}

# Get gene lists
gene_lists <- find_gene_lists(args$input, args$filter_by, args$threshold, args$lfc)

if (length(gene_lists) == 0) {
  write.csv(data.frame(message = "No gene lists found"),
            file.path(args$outdir, "kegg_no_results.csv"), row.names = FALSE)
  if (exists("gdb")) close_gene_db(gdb)
  quit(save = "no", status = 0)
}

message("Processing ", length(gene_lists), " gene lists with ", args$cpus, " worker(s)")

process_one <- function(label) {
  genes <- gene_lists[[label]]

  out_sub <- file.path(args$outdir, label)
  dir.create(out_sub, recursive = TRUE, showWarnings = FALSE)

  if (length(genes) < 3) {
    write.csv(data.frame(message = paste("Too few genes:", length(genes))),
              file.path(out_sub, "kegg_skipped.csv"), row.names = FALSE)
    return(invisible(NULL))
  }

  message("KEGG: ", label, " (", length(genes), " genes)")

  mapping <- convert_ids(genes)
  if (nrow(mapping) == 0) {
    write.csv(data.frame(message = "No ID mappings"),
              file.path(out_sub, "kegg_skipped.csv"), row.names = FALSE)
    return(invisible(NULL))
  }

  tryCatch({
    if (!is.null(kegg_data) && "kegg_id" %in% colnames(mapping)) {
      kk <- enricher(gene = unique(mapping$kegg_id), universe = universe_kegg,
                     TERM2GENE = kegg_data$KEGGPATHID2EXTID,
                     TERM2NAME = kegg_data$KEGGPATHID2NAME,
                     pAdjustMethod = "BH", pvalueCutoff = args$pvalue, qvalueCutoff = args$qvalue)
    } else {
      entrez_ids <- if ("ENTREZID" %in% colnames(mapping)) unique(mapping$ENTREZID) else NULL
      if (is.null(entrez_ids)) { message("  No Entrez IDs, skipping"); return(invisible(NULL)) }
      kk <- enrichKEGG(gene = entrez_ids, organism = args$kegg_org, universe = universe_entrez,
                       pAdjustMethod = "BH", pvalueCutoff = args$pvalue, qvalueCutoff = args$qvalue)
    }

    if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
      res_df <- as.data.frame(kk) %>%
        mutate(gene_list = label,
               fold_enrichment = as.numeric(sub("/.*", "", GeneRatio)) /
                                 as.numeric(sub(".*/", "", GeneRatio)) *
                                 as.numeric(sub(".*/", "", BgRatio)) /
                                 as.numeric(sub("/.*", "", BgRatio)))
      write.csv(res_df, file.path(out_sub, "kegg_results.csv"), row.names = FALSE)

      top <- res_df %>% arrange(p.adjust) %>% head(args$top_n)
      if (nrow(top) > 0) {
        p <- ggplot(top, aes(x = fold_enrichment, y = reorder(Description, fold_enrichment),
                             size = Count, color = p.adjust)) +
          geom_point() + scale_color_gradient(low = "red", high = "blue") +
          labs(x = "Fold Enrichment", y = NULL, title = paste("KEGG -", label)) +
          theme_minimal(base_size = 9)
        ggsave(file.path(out_sub, "kegg_dotplot.pdf"), plot = p,
               width = 10, height = max(3, nrow(top) * 0.25))
      }
    } else {
      write.csv(data.frame(message = "No significant pathways"),
                file.path(out_sub, "kegg_results.csv"), row.names = FALSE)
    }
  }, error = function(e) {
    message("  [", label, "] Error: ", e$message)
    write.csv(data.frame(message = paste("Error:", e$message)),
              file.path(out_sub, "kegg_error.csv"), row.names = FALSE)
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
  write.csv(data.frame(message = "No KEGG results"), file.path(args$outdir, "kegg_no_results.csv"), row.names = FALSE)
}

if (exists("gdb")) close_gene_db(gdb)
message("Done. Output: ", args$outdir)
