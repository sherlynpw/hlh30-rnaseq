#!/usr/bin/env Rscript
# Render sample correlation heatmaps from VST count matrices.
# - Replicates within (genotype, stage) collapsed by per-gene mean.
# - Two heatmaps emitted: from corrected and raw VST.
# - Pearson by default; Spearman/Kendall via --method.

suppressPackageStartupMessages({
  library(optparse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(readr)
  library(dplyr)
  library(matrixStats)
})

# >= 10pt in print: see _common/typography.R for FIG_SCALE / FONT_MIN.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "typography.R"))

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "palette.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--vst-corrected", type = "character", dest = "vst_corr"),
  make_option("--vst-raw",       type = "character", dest = "vst_raw"),
  make_option("--log2-norm",     type = "character", dest = "log2_norm",
              help = "DESeq2 normalized counts CSV; will be log2(x+1)-transformed"),
  make_option("--meta",          type = "character"),
  make_option("--method",        type = "character", default = "pearson"),
  make_option("--top-n",         type = "integer",   default = 0, dest = "top_n",
              help = "0 = all genes; otherwise top-N rows by variance"),
  make_option("--outdir",        type = "character"),
  make_option("--font-scale",    type = "double", default = 1, dest = "font_scale",
              help = "Multiply the pheatmap font size. Leave width/height at auto."),
  make_option("--width",         type = "double", default = 0,
              help = "Canvas width in inches; 0 = auto-fit. Set to the placed width for 1:1 layout."),
  make_option("--height",        type = "double", default = 0),
  make_option("--legend",        type = "character", default = "horizontal",
              help = "Colour bar orientation under the row labels: horizontal or vertical."),
  make_option("--mut-label", type = "character", default = "hlh-30(tm1978)", dest = "mut_label",
              help = "mutant name in the row and column labels (shorter for narrow panels)"),
  make_option("--cell-lines",    type = "double", default = 1.5, dest = "cell_lines",
              help = "Cell size in text lines. Raise it to fill a wider placement without enlarging the type.")
)))

stopifnot(!is.null(opt$vst_corr), !is.null(opt$vst_raw),
          !is.null(opt$meta), !is.null(opt$outdir))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

meta <- read_tsv(opt$meta, show_col_types = FALSE)
stage_short <- c(L1 = "L1", L4 = "L4", Day3 = "D3", Day6 = "D6")

# Accepts both new (WT/tm1978) and legacy (N2/hlh-30) tokens.
geno_label <- c(N2 = "WT", WT = "WT",
                `hlh-30` = opt$mut_label, tm1978 = opt$mut_label)

relabel <- function(x, map) ifelse(x %in% names(map), map[x], x)

meta$group <- paste(relabel(meta$genotype, geno_label),
                    relabel(meta$stage,    stage_short))

geno_order  <- c("WT", opt$mut_label)
stage_order <- c("L1", "L4", "D3", "D6")
ord <- expand.grid(stage = stage_order, genotype = geno_order,
                   stringsAsFactors = FALSE)
group_order <- paste(ord$genotype, ord$stage)
group_order <- intersect(group_order, unique(meta$group))

aggregate_vst <- function(path, transform = identity) {
  m <- read_csv(path, show_col_types = FALSE)
  gene_col <- names(m)[1]
  genes <- m[[gene_col]]
  mat <- as.matrix(m[, -1, drop = FALSE])
  mat <- transform(mat)
  rownames(mat) <- genes

  # Collapse replicates: per-gene mean across samples in each group.
  groups <- meta$group[match(colnames(mat), meta$sample_id)]
  if (any(is.na(groups))) {
    warning("samples missing from meta: ",
            paste(colnames(mat)[is.na(groups)], collapse = ", "))
    mat <- mat[, !is.na(groups)]
    groups <- groups[!is.na(groups)]
  }
  uniq <- intersect(group_order, unique(groups))
  agg <- sapply(uniq, function(g) rowMeans(mat[, groups == g, drop = FALSE]))
  colnames(agg) <- uniq
  agg
}

render <- function(mat, label) {
  if (opt$top_n > 0 && opt$top_n < nrow(mat)) {
    vars <- rowVars(mat)
    keep <- order(vars, decreasing = TRUE)[seq_len(opt$top_n)]
    mat <- mat[keep, , drop = FALSE]
  }
  cm <- cor(mat, method = opt$method)

  # Diverging blue->black->yellow palette over full correlation range [-1, 1],
  # matching the reference figure style.
  breaks <- seq(-1, 1, length.out = 101)
  pal    <- colorRampPalette(c("#2C5F9E", "#000000", "#F2D33C"))(100)

  # pheatmap's own clustering orders leaves purely by the tree, which scrambles
  # the stage sequence. Build the same tree, then rotate branches toward
  # group_order (L1 -> L4 -> D3 -> D6); topology is unchanged, only branch
  # orientation, so the leaf order follows the stages as far as the tree allows.
  hc <- hclust(dist(cm), method = "average")
  hc <- as.hclust(reorder(as.dendrogram(hc),
                          wts = match(colnames(cm), group_order),
                          agglo.FUN = mean))

  # "<genotype> <stage>" -> gene name italic, stage roman. WT is not a gene.
  lab_expr <- parse(text = vapply(colnames(cm), function(s) {
    p     <- strsplit(s, " ")[[1]]
    stage <- p[length(p)]
    geno  <- paste(p[-length(p)], collapse = " ")
    if (geno == "WT") sprintf('"%s %s"', geno, stage)
    else              sprintf('italic("%s")~"%s"', geno, stage)
  }, character(1)))

  base <- file.path(opt$outdir, sprintf("corr_%s_%s", label, opt$method))
  n  <- nrow(cm)
  fs <- 10 * opt$font_scale                              # label size in R pt
  gp <- gpar(fontsize = fs)
  # Cells scale with the type (1.5 lines each) so the row labels never crowd,
  # whatever the font scale. Dendrograms top and left, labels right and bottom.
  cell  <- unit(opt$cell_lines * fs, "pt")
  col_fun <- colorRamp2(c(-1, 0, 1), unname(PAL$div))
  ht <- Heatmap(
    cm, name = "r", col = col_fun,
    cluster_rows = hc, cluster_columns = hc,
    row_dend_side = "left", column_dend_side = "top",
    row_dend_width = unit(2.2 * fs, "pt"), column_dend_height = unit(2.2 * fs, "pt"),
    row_labels = lab_expr, column_labels = lab_expr,
    row_names_side = "right", column_names_side = "bottom",
    row_names_gp = gp, column_names_gp = gp, column_names_rot = 90,
    width = n * cell, height = n * cell,
    rect_gp = gpar(col = NA), show_heatmap_legend = FALSE
  )
  # Colour bar goes UNDER the row labels, beside the column labels, so it does
  # not add a column to the panel width. Its width is the row-label width.
  lab_w <- max_text_width(lab_expr, gp = gp)
  vertical <- opt$legend == "vertical"
  # "Pearson r" with the statistic in italic, centred over the bar
  lgd_title <- bquote(.(tools::toTitleCase(opt$method)) ~ italic(r))
  # vertical: the title is drawn by hand, centred on the BAR (ComplexHeatmap
  # centres it on bar + labels)
  bar_w <- unit(0.7 * fs, "pt")
  lgd <- if (vertical)
    Legend(col_fun = col_fun, at = seq(-1, 1, 0.5), direction = "vertical",
           title = NULL, labels_gp = gp, legend_height = lab_w - unit(1.6 * fs, "pt"),
           grid_width = bar_w, border = NA)
  else
    Legend(col_fun = col_fun, at = c(-1, 0, 1), direction = "horizontal",
           title = lgd_title, title_gp = gp, title_position = "topcenter",
           labels_gp = gp, legend_width = lab_w,
           grid_height = unit(0.7 * fs, "pt"), border = NA)
  # auto canvas: body + dendrograms + labels, plus a little padding
  pad_in <- 0.15
  w_in <- if (opt$width  > 0) opt$width  else
    convertWidth(n * cell + unit(2.2 * fs, "pt") + lab_w + unit(2 * fs, "pt"), "in", valueOnly = TRUE) + pad_in
  h_in <- if (opt$height > 0) opt$height else
    convertHeight(n * cell + unit(2.2 * fs, "pt") + lab_w + unit(2 * fs, "pt"), "in", valueOnly = TRUE) + pad_in
  # draw the legend in the column-names strip, just right of the heatmap body
  draw_all <- function() {
    draw(ht, padding = unit(c(2, 2, 2, 2), "pt"))
    decorate_column_names("r", {
      pushViewport(viewport(x = unit(1, "npc") + unit(fs, "pt"), y = unit(1, "npc") - unit(0.5 * fs, "pt"),
                            width = lab_w, height = lab_w, just = c("left", "top"),
                            clip = "off"))
      if (vertical) {
        # title centred in the row-label column; the bar sits directly under
        # the title's centre, its tick labels to the right
        tw <- convertWidth(grobWidth(textGrob(lgd_title, gp = gp)), "pt")
        cx <- 0.5 * max(tw, unit(0, "pt")) + unit(0.5 * fs, "pt")
        grid.text(lgd_title, x = cx, y = unit(1, "npc"), just = c("centre", "top"), gp = gp)
        draw(lgd, x = cx - 0.5 * bar_w, y = unit(1, "npc") - unit(1.4 * fs, "pt"), just = c("left", "top"))
      } else {
        draw(lgd, x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"))
      }
      popViewport()
    })
  }
  pdf(paste0(base, ".pdf"), width = w_in, height = h_in); draw_all(); dev.off()
  png(paste0(base, ".png"), width = w_in, height = h_in, units = "in", res = 200); draw_all(); dev.off()
  message("[wrote] ", base, ".{pdf,png}  (", nrow(mat), " genes)")
}

# RUV-corrected VST, for comparison against the log2norm version.
render(aggregate_vst(opt$vst_corr), "corrected")
# render(aggregate_vst(opt$vst_raw),  "raw")
if (!is.null(opt$log2_norm)) {
  render(aggregate_vst(opt$log2_norm, transform = function(x) log2(x + 1)),
         "log2norm")
  # Linear-scale Pearson on normalized counts (no log) — closer to the
  # TMM+Pearson recipe Plasmidsaurus describes; high-expression genes dominate
  # so stage-vs-stage drops sharply rather than sitting at 0.98+.
  # render(aggregate_vst(opt$log2_norm), "linnorm")
}
