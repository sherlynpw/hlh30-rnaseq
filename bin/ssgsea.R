#!/usr/bin/env Rscript

# ssgsea.R - Per-sample gene-set scoring via GSVA's ssGSEA method
#
# Inputs:
#   --vst_counts          DESeq2 VST matrix (RUV/batch-corrected if available)
#   --sample_meta         Sample metadata TSV (rows = samples, must have sample_id)
#   --builtin             Same builtins as gsea.R (go, kegg)
#   --collections_manifest JSON manifest of external GMTs (same format as gsea.R)
#
# Outputs per collection:
#   <collection>/scores.csv               pathways x samples
#   <collection>/pathway_pca.pdf          PCA in pathway space, colored per meta column
#   <collection>/top_var_heatmap.pdf      heatmap of top-N most variable pathways
#   <collection>/limma/<coef>.csv         toptables when --formula is provided

suppressPackageStartupMessages({
  library(optparse)
  library(GSVA)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(jsonlite)
  library(RSQLite)
  library(AnnotationDbi)
  library(org.Ce.eg.db)
  library(clusterProfiler)
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
  make_option("--vst_counts",   type = "character", help = "VST counts CSV (genes x samples)"),
  make_option("--sample_meta",  type = "character", help = "Sample metadata TSV"),
  make_option("--gene_db",      type = "character", default = NULL, help = "gene_db.sqlite"),
  make_option("--org_db",       type = "character", default = "org.Ce.eg.db"),
  make_option("--kegg_org",     type = "character", default = "cel"),
  make_option("--builtin",      type = "character", default = "",
              help = "Comma-separated built-in collections: go, kegg"),
  make_option("--collections_manifest", type = "character", default = NULL),
  make_option("--collections_dir", type = "character", default = "."),
  make_option("--ortholog_map", type = "character", default = NULL),
  make_option("--ortholog_from", type = "character", default = "wormbase_gene"),
  make_option("--ortholog_to",   type = "character", default = "human_symbol"),
  make_option("--min_size",     type = "integer", default = 15),
  make_option("--max_size",     type = "integer", default = 500),
  make_option("--method",       type = "character", default = "ssgsea",
              help = "GSVA method: ssgsea, gsva, plage, zscore"),
  make_option("--formula",      type = "character", default = NULL,
              help = "limma design formula e.g. '~ genotype * timepoint + replicate' (NULL = skip stats)"),
  make_option("--coefficients", type = "character", default = "",
              help = "Comma-separated coefficient names to extract (empty = all non-intercept)"),
  make_option("--top_n",        type = "integer", default = 30,
              help = "Top-N most variable pathways for heatmap"),
  make_option("--color_by",     type = "character", default = "",
              help = "Comma-separated meta columns to color PCA by"),
  make_option("--reference_levels", type = "character", default = "",
              help = "Comma-separated col=level pairs setting factor reference (e.g. 'genotype=N2,timepoint=L4')"),
  make_option("--top_pathways", type = "integer", default = 12,
              help = "Default top-N pathways for boxplots/trajectories (panel grids)"),
  make_option("--dotplot_top_n", type = "integer", default = -1,
              help = "Top-N pathways for dotplot (-1 = use --top_pathways)"),
  make_option("--volcano_label_n", type = "integer", default = -1,
              help = "Top-N pathways to label on volcano (-1 = use --top_pathways)"),
  make_option("--trajectory_x",     type = "character", default = "",
              help = "Meta column for x-axis of trajectory plots (e.g. 'timepoint')"),
  make_option("--trajectory_group", type = "character", default = "",
              help = "Meta column for line color in trajectory plots (e.g. 'genotype')"),
  make_option("--box_x", type = "character", default = "",
              help = "Comma-separated meta columns; combined (interaction) as x-axis for boxplots"),
  make_option("--contrasts_json", type = "character", default = NULL,
              help = "Path to JSON file: { name: 'coef + coef:term', ... }; tested via limma::makeContrasts"),
  make_option("--de_results_dir", type = "character", default = NULL,
              help = "Directory containing DESeq2 result CSVs (default + compound contrasts). When a CSV named after the ssGSEA term exists here, gene-driver bar plots are rendered for top pathways."),
  make_option("--de_mapping_json", type = "character", default = NULL,
              help = "JSON file mapping ssGSEA term name -> DE CSV basename (without .csv). Used when names don't match directly, e.g. 'genotypemut' -> 'mut_L4_vs_N2_L4'."),
  make_option("--driver_top_n", type = "integer", default = 15,
              help = "Top genes (by |log2FC|) to show per pathway in gene-driver plots"),
  make_option("--column_order_cols", type = "character", default = "",
              help = "Comma-separated meta columns combined with '_' to build the sample group key"),
  make_option("--column_order", type = "character", default = "",
              help = "Comma-separated group values in desired order; disables column clustering when set"),
  make_option("--outdir",       type = "character", default = ".")
)

args <- parse_args(OptionParser(option_list = opts))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
org_db <- get(args$org_db)

split_csv <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character())
  v <- trimws(strsplit(x, ",", fixed = TRUE)[[1]]); v[nzchar(v)]
}
is_true <- function(x) {
  if (is.logical(x)) return(isTRUE(x))
  tolower(trimws(as.character(x %||% "false"))) %in% c("true", "t", "1", "yes", "y")
}

# ---- Collection loaders (mirrors gsea.R) ----

read_gmt_sets <- function(path) {
  lines <- readLines(path, warn = FALSE); lines <- lines[nzchar(lines)]
  sets <- lapply(lines, function(l) {
    f <- strsplit(l, "\t", fixed = TRUE)[[1]]
    g <- unique(f[-c(1, 2)]); g[nzchar(g)]
  })
  names(sets) <- vapply(strsplit(lines, "\t", fixed = TRUE), `[`, character(1), 1)
  sets[lengths(sets) > 0]
}

read_gene_list_sets <- function(path, set_name) {
  g <- unique(trimws(readLines(path, warn = FALSE))); g <- g[nzchar(g)]
  if (!length(g)) return(list())
  setNames(list(g), set_name)
}

read_mapping_table <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  if (grepl("\\.csv$", path, ignore.case = TRUE)) read.csv(path, stringsAsFactors = FALSE)
  else read.delim(path, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
}

map_to_source_ids <- function(pathways, mapping_df, src, tgt) {
  if (is.null(mapping_df)) return(pathways)
  if (!all(c(src, tgt) %in% colnames(mapping_df)))
    stop("Ortholog map missing columns: ", src, ", ", tgt)
  m <- mapping_df[, c(src, tgt)]
  m <- m[!is.na(m[[src]]) & !is.na(m[[tgt]]), , drop = FALSE]
  m[[src]] <- trimws(as.character(m[[src]])); m[[tgt]] <- trimws(as.character(m[[tgt]]))
  m <- m[nzchar(m[[src]]) & nzchar(m[[tgt]]), , drop = FALSE]
  t2s <- split(m[[src]], m[[tgt]])
  out <- lapply(pathways, function(g) unique(unlist(t2s[g], use.names = FALSE)))
  out[lengths(out) > 0]
}

build_builtin <- function(kind, wb_entrez_map, kegg_org) {
  if (kind == "go") return(tryCatch({
    g2g <- AnnotationDbi::select(org_db, keys = keys(org_db, "WORMBASE"),
                                 keytype = "WORMBASE", columns = "GO")
    g2g <- g2g[!is.na(g2g$GO), c("GO", "WORMBASE")]
    split(g2g$WORMBASE, g2g$GO)
  }, error = function(e) { message("GO build failed: ", e$message); list() }))
  if (kind == "kegg") return(tryCatch({
    kd <- download_KEGG(kegg_org)
    if (!length(wb_entrez_map)) stop("No WB-Entrez map")
    e2wb <- setNames(names(wb_entrez_map), wb_entrez_map)
    km <- kd$KEGGPATHID2EXTID; km$wb <- e2wb[km$to]
    km <- km[!is.na(km$wb), ]; split(km$wb, km$from)
  }, error = function(e) { message("KEGG build failed: ", e$message); list() }))
  stop("Unknown builtin: ", kind)
}

# ---- Build collections ----

if (!is.null(args$gene_db) && file.exists(args$gene_db)) {
  gdb <- open_gene_db(args$gene_db); m <- get_all_mappings(gdb); close_gene_db(gdb)
  wb_entrez_map <- setNames(m$ENTREZID, m$WORMBASE)
} else {
  m <- AnnotationDbi::select(org_db, keys = keys(org_db, "WORMBASE"),
                              keytype = "WORMBASE", columns = "ENTREZID")
  m <- m[!is.na(m$ENTREZID), ]
  wb_entrez_map <- setNames(m$ENTREZID, m$WORMBASE)
}

mapping_df <- read_mapping_table(args$ortholog_map)
collections <- list()
for (k in split_csv(args$builtin)) {
  if (!k %in% c("go", "kegg")) { message("Unknown builtin: ", k); next }
  collections[[k]] <- build_builtin(k, wb_entrez_map, args$kegg_org)
}

if (!is.null(args$collections_manifest) && nzchar(args$collections_manifest) &&
    file.exists(args$collections_manifest)) {
  mf <- jsonlite::fromJSON(args$collections_manifest, simplifyDataFrame = FALSE)
  for (e in mf) {
    nm <- e$name; if (is.null(nm) || !nzchar(nm)) next
    fp <- file.path(args$collections_dir, e$file)
    fmt <- tolower(e$format %||% "gmt")
    use_map <- if (is_true(e$map_ortholog)) mapping_df else NULL
    collections[[nm]] <- tryCatch({
      raw <- if (fmt == "list") read_gene_list_sets(fp, nm) else read_gmt_sets(fp)
      map_to_source_ids(raw, use_map, args$ortholog_from, args$ortholog_to)
    }, error = function(err) { message("[", nm, "] ", err$message); list() })
  }
}

for (n in names(collections))
  message("Collection ", n, ": ", length(collections[[n]]), " pathways")

# Optional DE mapping (ssGSEA term name -> DE CSV basename)
de_mapping <- if (!is.null(args$de_mapping_json) && nzchar(args$de_mapping_json) &&
                  file.exists(args$de_mapping_json)) {
  jsonlite::fromJSON(args$de_mapping_json, simplifyDataFrame = FALSE)
} else list()
de_dir <- if (!is.null(args$de_results_dir) && nzchar(args$de_results_dir) &&
              dir.exists(args$de_results_dir)) args$de_results_dir else NULL
if (!is.null(de_dir)) {
  message("Gene-driver plots: DE dir = ", de_dir,
          " (", length(list.files(de_dir, pattern = "\\.csv$")), " CSVs)")
}

# Inline reorder_within / scale_x_reordered (avoid tidytext dependency).
# Lets us order bars per-facet inside a facet_wrap without dragging in extra packages.
.reorder_within <- function(x, by, within, fun = mean, sep = "___") {
  reorder(paste(x, within, sep = sep), by, FUN = fun)
}
.scale_x_reordered <- function(..., sep = "___") {
  scale_x_discrete(labels = function(x) gsub(paste0(sep, ".+$"), "", x), ...)
}

# ---- Load expression + meta ----

mat <- as.matrix(read.csv(args$vst_counts, row.names = 1, check.names = FALSE))
mat <- mat[rowSums(!is.na(mat)) == ncol(mat), , drop = FALSE]
meta <- read.delim(args$sample_meta, stringsAsFactors = FALSE, check.names = FALSE)
if (!"sample_id" %in% colnames(meta)) stop("sample_meta missing sample_id column")
meta <- meta[match(colnames(mat), meta$sample_id), , drop = FALSE]
if (any(is.na(meta$sample_id))) stop("sample_meta does not cover all vst columns")
rownames(meta) <- meta$sample_id

# Apply user-specified factor reference levels so limma coefficient names are predictable.
for (pair in split_csv(args$reference_levels)) {
  kv <- strsplit(pair, "=", fixed = TRUE)[[1]]
  if (length(kv) != 2) { message("Bad reference_levels entry: ", pair); next }
  col <- trimws(kv[1]); ref <- trimws(kv[2])
  if (!col %in% colnames(meta)) { message("reference_levels: col not in meta: ", col); next }
  vals <- as.character(meta[[col]])
  if (!ref %in% vals) { message("reference_levels: '", ref, "' not in meta$", col); next }
  meta[[col]] <- relevel(factor(vals), ref = ref)
  message("Factor ", col, " releveled with reference=", ref)
}

# Reorder samples by user-specified condition order (e.g. clustering.condition_order).
ordered_cols <- FALSE
order_cols <- split_csv(args$column_order_cols)
order_vals <- split_csv(args$column_order)
if (length(order_vals) && length(order_cols)) {
  missing_meta <- setdiff(order_cols, colnames(meta))
  if (length(missing_meta)) {
    message("column_order_cols not in meta, skipping reorder: ", paste(missing_meta, collapse = ","))
  } else {
    sample_key <- do.call(paste, c(lapply(order_cols, function(c) as.character(meta[[c]])), sep = "_"))
    key_rank <- match(sample_key, order_vals)
    if (any(is.na(key_rank))) {
      unmatched <- unique(sample_key[is.na(key_rank)])
      message("column_order: samples with keys not in order list (will trail): ",
              paste(unmatched, collapse = ", "))
    }
    # Stable order: matched samples by rank, unmatched at the end in original order
    new_order <- order(key_rank, seq_along(key_rank), na.last = TRUE)
    meta <- meta[new_order, , drop = FALSE]
    mat  <- mat[, meta$sample_id, drop = FALSE]
    rownames(meta) <- meta$sample_id
    ordered_cols <- TRUE
    message("Columns reordered by ", paste(order_cols, collapse = "_"),
            " using ", length(order_vals), "-level order")
  }
}

message("Expression matrix: ", nrow(mat), " genes x ", ncol(mat), " samples")

# ---- ssGSEA per collection ----

run_one <- function(name, sets) {
  out_dir <- file.path(args$outdir, name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!length(sets)) {
    message("[", name, "] no pathways, skipping"); return(invisible())
  }
  message("[", name, "] scoring ", length(sets), " pathways via ", args$method)

  # GSVA 1.50+ uses param objects; older versions accept the legacy signature.
  scores <- tryCatch({
    par <- ssgseaParam(exprData = mat, geneSets = sets,
                       minSize = args$min_size, maxSize = args$max_size)
    gsva(par, verbose = FALSE)
  }, error = function(e) {
    gsva(mat, sets, method = args$method,
         min.sz = args$min_size, max.sz = args$max_size, verbose = FALSE)
  })

  write.csv(scores, file.path(out_dir, "scores.csv"))

  # Top-N variable heatmap
  if (nrow(scores) >= 2) {
    rv <- apply(scores, 1, var, na.rm = TRUE)
    keep <- head(order(rv, decreasing = TRUE), min(args$top_n, nrow(scores)))
    if (length(keep) >= 2) {
      ann_cols <- setdiff(colnames(meta), c("sample_id"))
      ann_df <- if (length(ann_cols)) meta[, ann_cols, drop = FALSE] else NULL
      pheatmap(scores[keep, , drop = FALSE],
               annotation_col = ann_df,
               cluster_cols = !ordered_cols,
               main = paste0(name, " ssGSEA (top ", length(keep), " variable)"),
               filename = file.path(out_dir, "top_var_heatmap.pdf"),
               width = max(7, 0.25 * ncol(scores) + 4),
               height = max(6, 0.2 * length(keep) + 3),
               fontsize_row = 6, silent = TRUE)
    }
  }

  # PCA in pathway space
  if (ncol(scores) >= 3 && nrow(scores) >= 2) {
    pc <- prcomp(t(scores), scale. = FALSE)
    pdf(file.path(out_dir, "pathway_pca.pdf"), width = 7, height = 6)
    color_cols <- split_csv(args$color_by)
    if (!length(color_cols)) color_cols <- setdiff(colnames(meta), "sample_id")
    for (cc in color_cols) {
      if (!cc %in% colnames(meta)) next
      df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                       group = meta[[cc]], sample = rownames(meta))
      p <- ggplot(df, aes(PC1, PC2, color = group, label = sample)) +
        geom_point(size = 3) + theme_bw() +
        labs(title = paste0(name, " ssGSEA PCA"), color = cc)
      print(p)
    }
    dev.off()
  }

  # ---- limma on score matrix ----
  if (!is.null(args$formula) && nzchar(args$formula)) {
    lm_dir <- file.path(out_dir, "limma")
    dir.create(lm_dir, showWarnings = FALSE)
    design <- tryCatch(model.matrix(as.formula(args$formula), data = meta),
                       error = function(e) { message("design error: ", e$message); NULL })
    if (is.null(design)) return(invisible())
    if (nrow(design) != ncol(scores)) {
      message("[", name, "] design rows != samples (NA in meta?); skipping limma")
      return(invisible())
    }
    fit <- eBayes(lmFit(scores[, rownames(design), drop = FALSE], design))

    # Render volcano/dotplot/boxplots/trajectories for one named term (coef or contrast).
    render_term <- function(label, tt) {
      safe <- gsub("[^A-Za-z0-9_.-]", "_", label)
      term_dir <- file.path(lm_dir, safe)
      dir.create(term_dir, showWarnings = FALSE, recursive = TRUE)
      write.csv(tt, file.path(term_dir, "toptable.csv"), row.names = FALSE)

      tt_v <- tt
      tt_v$neglogp <- -log10(pmax(tt_v$adj.P.Val, 1e-300))
      tt_v$sig <- ifelse(tt_v$adj.P.Val < 0.05,
                         ifelse(tt_v$logFC > 0, "up", "down"), "ns")
      vol_n <- if (args$volcano_label_n > 0) args$volcano_label_n else args$top_pathways
      lab_n <- min(vol_n, nrow(tt_v))
      tt_v$label <- ifelse(seq_len(nrow(tt_v)) <= lab_n, tt_v$pathway, "")
      vp <- ggplot(tt_v, aes(logFC, neglogp, color = sig, label = label)) +
        geom_point(alpha = 0.6, size = 1.6) +
        scale_color_manual(values = c(up = "#c0392b", down = "#2980b9", ns = "grey70")) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
        ggrepel::geom_text_repel(size = 2.6, max.overlaps = 20, segment.size = 0.2) +
        theme_bw() + labs(title = paste0(name, " — ", label),
                          x = "logFC (Δ ssGSEA score)", y = "-log10 adj.P")
      ggsave(file.path(term_dir, "volcano.pdf"), vp, width = 7, height = 6)

      dp_n  <- if (args$dotplot_top_n > 0) args$dotplot_top_n else args$top_pathways
      top_n <- min(dp_n, nrow(tt))
      if (top_n >= 1) {
        td <- tt[seq_len(top_n), ]
        td$pathway <- factor(td$pathway, levels = rev(td$pathway))
        td$neglogp <- -log10(pmax(td$adj.P.Val, 1e-300))
        dp <- ggplot(td, aes(x = logFC, y = pathway, size = neglogp, color = logFC)) +
          geom_point() +
          scale_color_gradient2(low = "#2980b9", mid = "grey85", high = "#c0392b", midpoint = 0) +
          geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
          theme_bw() + labs(title = paste0(name, " — ", label, " (top ", top_n, ")"),
                            x = "logFC (Δ ssGSEA score)", size = "-log10 adj.P")
        ggsave(file.path(term_dir, "dotplot.pdf"), dp,
               width = 8, height = max(4, 0.3 * top_n + 2))
      }

      top_paths <- head(tt$pathway, args$top_pathways)
      if (!length(top_paths)) return(invisible())
      long <- data.frame(
        pathway = rep(top_paths, each = ncol(scores)),
        sample  = rep(colnames(scores), times = length(top_paths)),
        score   = as.vector(t(scores[top_paths, , drop = FALSE]))
      )
      long <- merge(long, cbind(sample = rownames(meta), meta), by = "sample")
      long$pathway <- factor(long$pathway, levels = top_paths)

      box_cols <- split_csv(args$box_x); box_cols <- box_cols[box_cols %in% colnames(meta)]
      if (length(box_cols)) {
        long$box_x <- do.call(paste, c(lapply(box_cols, function(c) long[[c]]), sep = "_"))
        color_col <- box_cols[1]
        bp <- ggplot(long, aes(x = box_x, y = score)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.6,
                       aes(fill = .data[[color_col]])) +
          geom_jitter(width = 0.15, size = 0.9, alpha = 0.8) +
          facet_wrap(~ pathway, scales = "free_y", ncol = 3) +
          theme_bw() + theme(axis.text.x = element_text(angle = 35, hjust = 1),
                             legend.position = "bottom") +
          labs(title = paste0(name, " — ", label, " (top pathways)"),
               x = paste(box_cols, collapse = " × "), y = "ssGSEA score", fill = color_col)
        ggsave(file.path(term_dir, "boxplots.pdf"), bp,
               width = max(8, length(unique(long$box_x)) * 0.6 + 4),
               height = max(5, ceiling(length(top_paths) / 3) * 2.5))
      }

      tx <- args$trajectory_x; tg <- args$trajectory_group
      if (nzchar(tx) && tx %in% colnames(meta)) {
        long$x_val <- long[[tx]]
        group_col <- if (nzchar(tg) && tg %in% colnames(meta)) tg else NULL
        agg_cols <- c("pathway", "x_val", group_col)
        agg <- aggregate(long$score, by = long[, agg_cols, drop = FALSE],
                         FUN = function(v) c(mean = mean(v), se = sd(v) / sqrt(length(v))))
        agg <- do.call(data.frame, agg)
        colnames(agg)[(ncol(agg) - 1):ncol(agg)] <- c("mean", "se")
        aes_map <- if (!is.null(group_col))
          aes(x = x_val, y = mean, color = .data[[group_col]], group = .data[[group_col]])
        else aes(x = x_val, y = mean, group = 1)
        tp <- ggplot(agg, aes_map) +
          geom_line(linewidth = 0.7) + geom_point(size = 2) +
          geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15, alpha = 0.7) +
          geom_jitter(data = long, aes(x = x_val, y = score,
                                       color = if (!is.null(group_col)) .data[[group_col]] else NULL),
                      inherit.aes = FALSE, width = 0.1, alpha = 0.3, size = 0.8) +
          facet_wrap(~ pathway, scales = "free_y", ncol = 3) +
          theme_bw() + labs(title = paste0(name, " — ", label, " trajectories"),
                            x = tx, y = "ssGSEA score (mean ± SE)", color = tg)
        ggsave(file.path(term_dir, "trajectories.pdf"), tp,
               width = max(8, length(unique(long$x_val)) * 1.2 + 4),
               height = max(5, ceiling(length(top_paths) / 3) * 2.5))
      }

      # Gene-driver bar plots: for the top pathways from this term, show each
      # pathway member's DESeq2 log2FC for the matching DE contrast. Answers
      # "which genes drove the pathway-level change?"
      if (!is.null(de_dir)) {
        de_name <- de_mapping[[label]] %||% label
        de_path <- file.path(de_dir, paste0(de_name, ".csv"))
        if (!file.exists(de_path)) {
          message("[", name, "/", label, "] no DE CSV at ", de_path,
                  " (set ssgsea.de_mapping or rename); skipping gene drivers")
        } else {
          de_df <- read.csv(de_path, stringsAsFactors = FALSE)
          if (!"gene_id" %in% colnames(de_df)) {
            message("[", name, "/", label, "] DE CSV missing gene_id column; skipping")
          } else {
            long_de <- do.call(rbind, lapply(as.character(top_paths), function(p) {
              members <- sets[[p]]
              if (!length(members)) return(NULL)
              sub <- de_df[de_df$gene_id %in% members & !is.na(de_df$log2FoldChange), ,
                           drop = FALSE]
              if (!nrow(sub)) return(NULL)
              sub <- sub[order(-abs(sub$log2FoldChange)), ]
              sub <- head(sub, args$driver_top_n)
              sub$pathway <- p
              sub
            }))
            if (!is.null(long_de) && nrow(long_de)) {
              # Display label: prefer gene_name (symbol), fall back to gene_id
              long_de$display <- if ("gene_name" %in% colnames(long_de))
                ifelse(is.na(long_de$gene_name) | !nzchar(long_de$gene_name),
                       long_de$gene_id, long_de$gene_name) else long_de$gene_id
              long_de$sig <- ifelse(is.na(long_de$padj) | long_de$padj >= 0.05, "ns",
                              ifelse(long_de$log2FoldChange > 0, "up", "down"))
              long_de$pathway <- factor(long_de$pathway, levels = top_paths)

              gp <- ggplot(long_de,
                           aes(x = .reorder_within(display, log2FoldChange, pathway),
                               y = log2FoldChange, fill = sig)) +
                geom_col() +
                scale_fill_manual(values = c(up = "#c0392b", down = "#2980b9", ns = "grey70")) +
                coord_flip() +
                facet_wrap(~ pathway, scales = "free_y", ncol = 3) +
                .scale_x_reordered() +
                geom_hline(yintercept = 0, color = "grey40", linewidth = 0.3) +
                theme_bw() +
                theme(legend.position = "bottom",
                      axis.text.y = element_text(size = 7)) +
                labs(title = paste0(name, " — ", label, " — gene drivers"),
                     subtitle = paste0("DE source: ", basename(de_path),
                                       "; top ", args$driver_top_n,
                                       " genes per pathway by |log2FC|"),
                     x = NULL, y = "log2FC (DESeq2)", fill = "padj < 0.05")
              ggsave(file.path(term_dir, "gene_drivers.pdf"), gp,
                     width = max(9, 8),
                     height = max(5, ceiling(length(top_paths) / 3) * 3))
            }
          }
        }
      }
    }

    # Helper: extract a toptable and tag the pathway column
    get_tt <- function(fit_obj, coef_arg) {
      tt <- topTable(fit_obj, coef = coef_arg, number = Inf, sort.by = "P")
      tt$pathway <- rownames(tt)
      tt[, c("pathway", setdiff(colnames(tt), "pathway"))]
    }

    # --- Raw coefficients ---
    coefs <- split_csv(args$coefficients)
    if (!length(coefs)) coefs <- setdiff(colnames(design), "(Intercept)")
    for (cf in coefs) {
      if (!cf %in% colnames(design)) {
        message("[", name, "] coef '", cf, "' not in design; skipping"); next
      }
      render_term(cf, get_tt(fit, cf))
    }

    # --- Named linear contrasts via limma::makeContrasts ---
    if (!is.null(args$contrasts_json) && nzchar(args$contrasts_json) &&
        file.exists(args$contrasts_json)) {
      contrasts_map <- jsonlite::fromJSON(args$contrasts_json, simplifyDataFrame = FALSE)
      if (length(contrasts_map)) {
        # Coefficient names containing ':' (interaction terms) collide with R's seq operator
        # when parsed by makeContrasts; backtick-quote them in the expression strings.
        needs_quote <- colnames(design)[grepl("[^A-Za-z0-9_.]", colnames(design))]
        quote_terms <- function(expr) {
          for (cn in needs_quote)
            expr <- gsub(cn, paste0("`", cn, "`"), expr, fixed = TRUE)
          expr
        }
        exprs <- vapply(contrasts_map, function(x) quote_terms(as.character(x)), character(1))
        cmat <- tryCatch(
          makeContrasts(contrasts = exprs, levels = design),
          error = function(e) { message("[", name, "] makeContrasts failed: ", e$message); NULL }
        )
        if (!is.null(cmat)) {
          colnames(cmat) <- names(contrasts_map)
          fit_c <- eBayes(contrasts.fit(fit, cmat))
          for (cn in colnames(cmat)) render_term(cn, get_tt(fit_c, cn))
        }
      }
    }
  }
}

for (n in names(collections)) run_one(n, collections[[n]])
message("ssGSEA complete.")
