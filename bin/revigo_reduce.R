#!/usr/bin/env Rscript

# revigo_reduce.R - REVIGO-style reduction of GO enrichment results via rrvgo.
#
# Reads all go_*.csv files produced by go_enrichment.R (recursively under
# --input), groups by (ontology, gene_list), and for each group:
#   1. Computes pairwise semantic similarity between GO terms (GOSemSim)
#   2. Reduces the set to representatives at a similarity threshold
#   3. Emits a reduced CSV + 5 plots
#
# Why per (ontology, gene_list): semantic similarity is only meaningful WITHIN
# a single ontology (BP terms aren't comparable to CC terms), and you don't
# want to merge GO results from different gene lists into a single reduction.
#
# Output layout mirrors the GO module: <gene_list>/revigo_<ont>_*.pdf|.csv
#
# Inputs:
#   --input    Directory with go_*.csv (recursive search). The CSVs are the
#              clusterProfiler output, with columns: ID, Description, pvalue,
#              p.adjust, GeneRatio, BgRatio, Count, ontology, gene_list, ...
#   --org_db   org.Ce.eg.db (or equivalent OrgDb package name)
#   --threshold rrvgo similarity threshold: 0.5 (broad), 0.7 (default),
#               0.9 (fine-grained). Lower = fewer, broader clusters.
#   --ontologies BP,MF,CC (default: BP,MF,CC; reduction only runs on these)
#   --outdir
#
# Outputs (per gene_list, per ontology):
#   <gene_list>/revigo_<ont>_reduced.csv   - all input GO terms with rrvgo
#                                            cluster assignment + original
#                                            enrichGO columns (Description,
#                                            GeneRatio, BgRatio, Count,
#                                            p.adjust, FoldEnrichment, geneID).
#                                            Rows sorted: cluster, then
#                                            representative first, then p.adjust.
#                                            is_representative=TRUE for the
#                                            cluster head (go == parent).
#   <gene_list>/revigo_<ont>_treemap.pdf   - cluster treemap (area = -log10(padj))
#   <gene_list>/revigo_<ont>_scatter.pdf   - 2D semantic-space scatter
#   <gene_list>/revigo_<ont>_wordcloud.pdf - representative terms as word cloud
#   <gene_list>/revigo_<ont>_heatmap.pdf   - full similarity matrix annotated by cluster
#   <gene_list>/revigo_<ont>_heatmap_representatives.pdf
#                                          - compact heatmap of just the cluster
#                                            representatives (one row/col per
#                                            rrvgo cluster head); much more
#                                            readable than the full matrix for
#                                            big GO sets
#   <gene_list>/revigo_<ont>_dotplot.pdf   - GO-style dotplot of representative
#                                            terms only (x=FoldEnrichment,
#                                            color=p.adjust, size=Count)

suppressPackageStartupMessages({
  library(optparse)
  library(rrvgo)
  library(dplyr)
  library(ggplot2)
})

opts <- list(
  make_option("--input",      type = "character", default = ".",
              help = "Directory with go_*.csv files (searched recursively)"),
  make_option("--org_db",     type = "character", default = "org.Ce.eg.db",
              help = "OrgDb package name (default org.Ce.eg.db)"),
  make_option("--threshold",  type = "numeric",   default = 0.7,
              help = "rrvgo similarity threshold 0.4-0.9 (default 0.7)"),
  make_option("--ontologies", type = "character", default = "BP,MF,CC",
              help = "Comma-separated GO ontologies to reduce (default BP,MF,CC)"),
  make_option("--min_terms",  type = "integer",   default = 3,
              help = "Skip groups with fewer than N GO terms (default 3)"),
  make_option("--cpus",       type = "integer",   default = 1,
              help = "Worker processes for parallel group reduction (default 1; mclapply)"),
  make_option("--outdir",     type = "character", default = ".",
              help = "Output directory")
)
args <- parse_args(OptionParser(option_list = opts))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# Load OrgDb dynamically
if (!requireNamespace(args$org_db, quietly = TRUE)) {
  stop("OrgDb package not installed: ", args$org_db)
}
suppressPackageStartupMessages(library(args$org_db, character.only = TRUE))

ontologies <- trimws(strsplit(args$ontologies, ",", fixed = TRUE)[[1]])

# Parse "n/N" string -> numeric ratio. Returns NA on parse failure.
parse_ratio <- function(s) {
  m <- regmatches(s, regexec("^(\\d+)/(\\d+)$", s))
  sapply(m, function(mi) {
    if (length(mi) != 3) return(NA_real_)
    num <- as.numeric(mi[2]); den <- as.numeric(mi[3])
    if (is.na(den) || den == 0) return(NA_real_)
    num / den
  })
}

# Find all go_*.csv files. Format expected: clusterProfiler enrichGO columns +
# our added 'ontology' and 'gene_list' columns (set by go_enrichment.R).
files <- list.files(args$input, pattern = "^go_.*\\.csv$",
                    recursive = TRUE, full.names = TRUE)
files <- files[!grepl("go_no_results|go_skipped", basename(files))]

if (length(files) == 0) {
  message("No GO result CSVs found under ", args$input, " - nothing to reduce.")
  writeLines("No GO results to reduce", file.path(args$outdir, "revigo_no_input.txt"))
  quit(save = "no", status = 0)
}

message("Found ", length(files), " GO result file(s); reducing per (ontology, gene_list)")

# Aggregate. Keep all columns we want downstream:
#   - ID/Description/GeneRatio/BgRatio/Count/p.adjust  -> dotplot + ranking
#   - geneID  -> per-term member gene list (slash-separated Entrez IDs from
#                clusterProfiler), copied through into revigo_*_reduced.csv
#                so users can trace from a representative term to its genes
keep_cols <- c("ID", "Description", "GeneRatio", "BgRatio", "Count",
               "pvalue", "p.adjust", "geneID", "ontology", "gene_list")
all_rows <- list()
for (f in files) {
  df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) next
  needed <- c("ID", "p.adjust", "ontology", "gene_list")
  if (!all(needed %in% colnames(df))) {
    message("  Skipping ", basename(f), " - missing cols: ",
            paste(setdiff(needed, colnames(df)), collapse = ", "))
    next
  }
  all_rows[[length(all_rows) + 1L]] <-
    df[, intersect(keep_cols, colnames(df)), drop = FALSE]
}

if (length(all_rows) == 0) {
  message("No usable GO results - skipping.")
  writeLines("No usable GO results", file.path(args$outdir, "revigo_no_input.txt"))
  quit(save = "no", status = 0)
}

combined <- dplyr::bind_rows(all_rows)
combined <- combined[combined$ontology %in% ontologies, ]
combined <- combined[!is.na(combined$p.adjust), ]

if (nrow(combined) == 0) {
  message("No rows match requested ontologies ", args$ontologies)
  writeLines("No matching rows", file.path(args$outdir, "revigo_no_input.txt"))
  quit(save = "no", status = 0)
}

# FoldEnrichment from GeneRatio / BgRatio (matches what enrichplot uses)
if (all(c("GeneRatio", "BgRatio") %in% colnames(combined))) {
  combined$FoldEnrichment <- parse_ratio(combined$GeneRatio) /
                             parse_ratio(combined$BgRatio)
}

# Split into (ontology, gene_list) groups
combined$group_key <- paste(combined$ontology, combined$gene_list, sep = "__")
groups <- split(combined, combined$group_key)

make_dotplot <- function(reduced, grp_df, ont, gene_list, top_n = 20) {
  # Filter to representative terms (parents of their own cluster).
  reps_ids <- reduced$go[reduced$go == reduced$parent]
  rep_df <- grp_df[grp_df$ID %in% reps_ids, , drop = FALSE]
  if (nrow(rep_df) == 0) return(NULL)

  # Need at least Description + p.adjust to plot
  if (!"Description" %in% colnames(rep_df)) {
    rep_df$Description <- rep_df$ID
  }

  rep_df <- rep_df %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n)

  has_fe <- "FoldEnrichment" %in% colnames(rep_df) &&
            any(is.finite(rep_df$FoldEnrichment))
  has_count <- "Count" %in% colnames(rep_df) &&
               any(!is.na(rep_df$Count))

  rep_df$Description <- factor(rep_df$Description,
                               levels = rev(rep_df$Description))

  x_aes <- if (has_fe) "FoldEnrichment" else "p.adjust"
  x_lab <- if (has_fe) "Fold Enrichment" else "p.adjust"

  p <- ggplot(rep_df, aes(x = .data[[x_aes]], y = Description,
                          color = p.adjust)) +
    (if (has_count) geom_point(aes(size = Count)) else geom_point(size = 3)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(title = paste0("REVIGO dotplot - GO ", ont, " - ", gene_list,
                        " (n=", nrow(rep_df), " representatives)"),
         x = x_lab, y = NULL,
         color = "p.adjust", size = "Count") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_line(color = "grey92"))

  p
}

run_one <- function(grp_name, grp_df) {
  ont       <- unique(grp_df$ontology)[1]
  gene_list <- unique(grp_df$gene_list)[1]

  # Per-gene-list subdir (mirrors the GO module's structure)
  gene_dir <- file.path(args$outdir,
                        gsub("[^A-Za-z0-9_.-]+", "_", gene_list))
  dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

  # File-name prefix per ontology
  base <- file.path(gene_dir, paste0("revigo_", ont))

  if (nrow(grp_df) < args$min_terms) {
    message("  [", ont, " ", gene_list, "] only ", nrow(grp_df),
            " terms (<", args$min_terms, "), skipping")
    return(invisible(NULL))
  }

  message("[", ont, " ", gene_list, "] reducing ", nrow(grp_df), " GO terms")

  go_ids <- grp_df$ID
  scores <- setNames(-log10(pmax(grp_df$p.adjust, .Machine$double.eps)), go_ids)

  sim <- tryCatch(
    calculateSimMatrix(go_ids, orgdb = args$org_db, ont = ont, method = "Rel"),
    error = function(e) {
      message("    calculateSimMatrix failed: ", conditionMessage(e)); NULL
    }
  )
  if (is.null(sim) || any(dim(sim) < 2)) {
    message("    similarity matrix too small, skipping")
    return(invisible(NULL))
  }

  reduced <- tryCatch(
    reduceSimMatrix(sim, scores = scores, threshold = args$threshold,
                    orgdb = args$org_db),
    error = function(e) {
      message("    reduceSimMatrix failed: ", conditionMessage(e)); NULL
    }
  )
  if (is.null(reduced) || nrow(reduced) == 0) {
    message("    no reduced terms produced")
    return(invisible(NULL))
  }

  reduced$ontology  <- ont
  reduced$gene_list <- gene_list

  # Enrich the reduced table with original enrichGO columns. rrvgo returns:
  #   go, parent, parentTerm, score, size, term, cluster
  # We add (from the original GO CSV, joined on ID == go):
  #   Description, GeneRatio, BgRatio, Count, pvalue, p.adjust, geneID,
  #   FoldEnrichment, is_representative
  # geneID is slash-separated Entrez IDs (clusterProfiler convention).
  go_extra_cols <- intersect(
    c("Description", "GeneRatio", "BgRatio", "Count", "pvalue", "p.adjust",
      "geneID", "FoldEnrichment"),
    colnames(grp_df))
  reduced_full <- merge(
    reduced,
    unique(grp_df[, c("ID", go_extra_cols), drop = FALSE]),
    by.x = "go", by.y = "ID", all.x = TRUE, sort = FALSE
  )
  reduced_full$is_representative <- reduced_full$go == reduced_full$parent

  # Sort: representatives first within each cluster, then by p.adjust
  reduced_full <- reduced_full[order(reduced_full$cluster,
                                     -reduced_full$is_representative,
                                     reduced_full$p.adjust), ]

  write.csv(reduced_full, paste0(base, "_reduced.csv"), row.names = FALSE)

  # Treemap
  tryCatch({
    pdf(paste0(base, "_treemap.pdf"), width = 10, height = 7)
    treemapPlot(reduced, title = paste("REVIGO treemap:", ont, gene_list))
    invisible(dev.off())
  }, error = function(e) message("    treemap failed: ", conditionMessage(e)))

  # Scatter
  tryCatch({
    p <- scatterPlot(sim, reduced, labelSize = 3)
    ggsave(paste0(base, "_scatter.pdf"), p, width = 10, height = 8)
  }, error = function(e) message("    scatter failed: ", conditionMessage(e)))

  # Word cloud
  tryCatch({
    pdf(paste0(base, "_wordcloud.pdf"), width = 8, height = 6)
    wordcloudPlot(reduced, min.freq = 1, colors = "black")
    invisible(dev.off())
  }, error = function(e) message("    wordcloud failed: ", conditionMessage(e)))

  # Full heatmap. Generous per-term sizing + font floor at 8 so labels stay
  # legible even for big GO sets. Caps at 60in (PDF readers struggle past that;
  # if you have >250 terms the per-cluster summary below is the better view).
  tryCatch({
    n_terms <- nrow(sim)
    side    <- min(60, max(8, 0.30 * n_terms + 3))
    fontsz  <- max(8, 12 - 0.02 * n_terms)
    pdf(paste0(base, "_heatmap.pdf"), width = side, height = side)
    heatmapPlot(sim, reduced,
                annotateParent = TRUE,
                annotationLabel = "parentTerm",
                fontsize = fontsz)
    invisible(dev.off())
  }, error = function(e) message("    heatmap failed: ", conditionMessage(e)))

  # Representatives-only heatmap: one row/col per rrvgo cluster head. Always
  # compact (typically 5-20 terms) and readable. The full heatmap above shows
  # within-cluster structure; this one shows between-cluster similarity.
  tryCatch({
    rep_ids <- unique(reduced$parent)
    common  <- intersect(rep_ids, rownames(sim))
    if (length(common) >= 2) {
      sim_reps <- sim[common, common, drop = FALSE]
      reduced_reps <- reduced[reduced$go %in% common, , drop = FALSE]
      n <- length(common)
      side2 <- max(7, 0.45 * n + 3)
      pdf(paste0(base, "_heatmap_representatives.pdf"),
          width = side2, height = side2)
      heatmapPlot(sim_reps, reduced_reps,
                  annotateParent = TRUE,
                  annotationLabel = "parentTerm",
                  fontsize = 11)
      invisible(dev.off())
    }
  }, error = function(e) message("    heatmap_reps failed: ", conditionMessage(e)))

  # GO-style dotplot of representative terms
  tryCatch({
    dp <- make_dotplot(reduced, grp_df, ont, gene_list)
    if (!is.null(dp)) {
      ggsave(paste0(base, "_dotplot.pdf"), dp,
             width = 10, height = max(4, 0.3 * length(unique(reduced$parent)) + 3))
    }
  }, error = function(e) message("    dotplot failed: ", conditionMessage(e)))

  invisible(NULL)
}

# Parallelize across groups via mclapply (forks). Each group's calculateSimMatrix
# + reduceSimMatrix + plotting runs on its own core; useful when there are many
# (ontology, gene_list) combinations (e.g. 5 contrasts x 3 ontologies = 15 groups
# for DESEQ-source REVIGO).
n_workers <- max(1L, min(as.integer(args$cpus), length(groups)))
if (n_workers > 1L && .Platform$OS.type != "windows") {
  suppressPackageStartupMessages(library(parallel))
  message("Running ", length(groups), " group(s) on ", n_workers, " workers (mclapply)")
  invisible(parallel::mclapply(
    names(groups),
    function(k) run_one(k, groups[[k]]),
    mc.cores = n_workers, mc.preschedule = FALSE
  ))
} else {
  for (k in names(groups)) {
    run_one(k, groups[[k]])
  }
}

message("Done. Output: ", args$outdir)
