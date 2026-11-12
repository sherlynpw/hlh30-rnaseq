#!/usr/bin/env Rscript
# GO dot plots, matching the wb.rs.ht app's EnrichmentDotPlot component.
#
# Same targets.tsv driver as render-bars.R. Encoding is taken from
# src/lib/components/EnrichmentDotPlot.svelte:
#   x      = fold enrichment (linear, 0 .. max(2, foldMax), 4 intervals)
#   radius = 4 + sqrt(k / kMax) * 7  px, k = gene count
#   colour = -log10(q) / logQMax on a 5-stop viridis ramp
#   order  = by q ascending (the component's default sortBy)
# so the R output reads the same as the web view.

suppressPackageStartupMessages({
  library(optparse); library(ggplot2); library(dplyr); library(readr); library(stringr)
})
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "typography.R"))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "go-plot-common.R"))

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "palette.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--bucket", type = "character", default = "s3://celegans-pipeline-rayandrew"),
  make_option("--cache",  type = "character", default = "cache"),
  make_option("--top-n",  type = "integer", default = 20, dest = "top_n",
              help = "When the terms column is empty, plot the top N most significant terms"),
  make_option("--sort-by", type = "character", default = "q", dest = "sort_by",
              help = "Row order: q (default, as the component), fold, or k"),
  make_option("--row-height", type = "double", default = 0.55, dest = "row_height",
              help = paste("Inches of canvas per GO term at font-scale 1.",
                           "0.55 is airy; ~0.30 is tight but still clears 10pt.")),
  make_option("--font-scale", type = "double", default = 1, dest = "font_scale",
              help = "Multiply every font size. Render at the ORIGINAL canvas and raise this."),
  make_option("--target-pt", type = "double", default = 0, dest = "target_pt",
              help = paste("Printed size of the term labels once placed. Needs a",
                           "placed_in column (inches) in the config; overrides",
                           "--font-scale per row. 0 = off.")),
  make_option("--width",  type = "double", default = 7.0,
              help = "Panel width in inches. Set to the placed width for 1:1 layout."),
  make_option("--height", type = "double", default = 0,
              help = "Panel height in inches; 0 = auto from the number of terms."),
  make_option("--legend", type = "character", default = "right",
              help = paste("right (default) or bottom. At paper width the",
                           "vertical colourbar costs ~25mm that the 10pt term",
                           "labels need, so put it under the plot.")),
  make_option("--x-breaks", type = "integer", default = 5, dest = "x_breaks",
              help = "Number of x-axis breaks; 5 collides at paper width."),
  make_option("--wrap",   type = "integer", default = 38,
              help = "Wrap GO term labels at this many characters. Narrow panels need a smaller value at 10pt."),
  make_option("--label-frac", type = "double", default = 0, dest = "label_frac",
              help = paste("Width of the term-label column as a share of the placed",
                           "width, so every panel of one width has the same label and",
                           "plot columns; labels wrap to fit. Needs --target-pt. 0 = off.")),
  make_option("--ramp",   type = "character", default = "viridis",
              help = "Colour ramp: viridis (as the component) or purple (brand tokens)")
)))
stopifnot(!is.null(opt$config), !is.null(opt$outdir))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(opt$cache,  showWarnings = FALSE, recursive = TRUE)

targets <- read_tsv(opt$config, show_col_types = FALSE)

# px -> mm diameter, so the app's radii carry over to ggplot point sizes
PX <- 25.4 / 96
r_for <- function(k, kmax) 4 + sqrt(k / max(1, kmax)) * (11 - 4)

# --- fixed-width layout helpers (--label-in) -----------------------------------
# Smallest round step (1, 1.5, 2, 2.5, 3, 4, 5, 6, 8 x 10^k) at or above x.
nice_step <- function(x) {
  m <- c(1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10) * 10^floor(log10(x))
  m[m >= x - 1e-9][1]
}
# Set the label column to label_w and give the panel what is left of total_w,
# so every panel has the same label and plot widths whatever its terms.
fix_widths <- function(p, label_w, total_w) {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  g  <- ggplotGrob(p)
  lc <- unique(g$layout$l[g$layout$name == "axis-l"])
  pc <- unique(g$layout$l[g$layout$name == "panel"])
  g$widths[lc] <- unit(label_w, "in")
  other <- sum(grid::convertWidth(g$widths[-pc], "in", TRUE))
  g$widths[pc] <- unit(total_w - other, "in")
  g
}

render_one <- function(row) {
  label <- sprintf("%s | %s | %s", row$run, row$path, row$ontology)
  # Term labels are axis.text = rel(0.8) of base_size 12, so a panel rendered
  # `width` in wide and placed `placed_in` wide prints them at
  # 9.6 * font_scale * placed_in / width pt. Solve for font_scale.
  if (opt$target_pt > 0) {
    stopifnot("placed_in" %in% names(row), is.finite(row$placed_in), row$placed_in > 0)
    opt$font_scale <- opt$target_pt / (9.6 * row$placed_in / opt$width)
    # pdf() writes whole-point sizes, so 19.5pt prints as 19: round the term
    # label size up to the next whole point
    opt$font_scale <- ceiling(9.6 * opt$font_scale) / 9.6
    message(sprintf("[%s] placed %.2f in -> font_scale %.3f", label, row$placed_in, opt$font_scale))
  }
  csv <- ensure_csv(row$run, row$path, row$ontology, opt$bucket, opt$cache)
  if (is.null(csv)) return(invisible(NULL))

  go  <- read_csv(csv, show_col_types = FALSE)
  fold_col <- if ("FoldEnrichment" %in% names(go)) "FoldEnrichment" else "fold_enrichment"
  if (!fold_col %in% names(go)) {
    warning(sprintf("[%s] no fold-enrichment column, skipping", label)); return(invisible(NULL))
  }

  hit <- select_terms(go, row$terms, opt$top_n, label)
  if (nrow(hit) == 0) {
    warning(sprintf("[%s] no terms matched, skipping figure", label)); return(invisible(NULL))
  }

  # fixed layout: canvas inches per printed inch, and the label column on the canvas
  fixed   <- opt$label_frac > 0 && opt$target_pt > 0
  # on the canvas, the label column is the same share of the width
  label_w <- opt$label_frac * opt$width
  label_pt <- 9.6 * opt$font_scale
  hit <- hit %>%
    mutate(fold      = .data[[fold_col]],
           neglog10q = -log10(p.adjust),
           # full GO name when it fits in 2 lines, else the short name from
           # targets.tsv; 8% slack because the column also holds the
           # label-to-panel margin and png sets text slightly wider than pdf
           display_wrapped = if (fixed) paper_label(Description, display, 0.92 * label_w, label_pt)
                             else str_wrap(display, width = opt$wrap)) %>%
    # ggplot puts the first level at the bottom, so order ascending to place the
    # most significant term at the top - matching the component's sort.
    arrange(dplyr::case_when(opt$sort_by == "fold" ~ fold,
                             opt$sort_by == "k"    ~ as.numeric(Count),
                             TRUE                  ~ neglog10q)) %>%
    mutate(display_wrapped = factor(display_wrapped, levels = unique(display_wrapped)))

  # scale bounds exactly as the component derives them
  fold_max <- max(2, hit$fold, na.rm = TRUE)
  k_max    <- max(1, hit$Count, na.rm = TRUE)
  logq_max <- max(1, hit$neglog10q, na.rm = TRUE)
  # panels of one scale_group share the colour range, so one bar serves all
  grp <- if ("scale_group" %in% names(row)) row$scale_group else NA
  if (!is.na(grp) && nzchar(grp)) logq_max <- GROUP_QMAX[[grp]]
  show_legend <- !("legend" %in% names(row)) || is.na(row$legend) || row$legend != "no"
  hit$size_mm <- 2 * r_for(hit$Count, k_max) * PX

  # Round tick positions (pretty), axis ends on the last tick
  if (fixed) {
    # every panel gets the same grid: 0, a round step and twice the step
    x_breaks <- c(0, 1, 2) * nice_step(fold_max / 2)
  } else {
    x_breaks <- scales::breaks_pretty(n = opt$x_breaks)(c(0, fold_max))
    x_breaks <- x_breaks[x_breaks >= 0]
    if (max(x_breaks) < fold_max) x_breaks <- c(x_breaks, max(x_breaks) + diff(x_breaks)[1])
  }
  x_max <- max(x_breaks)
  # colour-bar ticks at both ends and the middle. The top is the real maximum
  # rounded UP to a whole number, so the midpoint lands on .0 or .5 and the
  # most significant term still gets (nearly) the top colour.
  q_top    <- ceiling(logq_max)
  q_breaks <- c(0, q_top / 2, q_top)
  q_labels <- sub("\\.0$", "", sprintf("%.1f", q_breaks))

  # Gene counts: with 3-5 genes per term the app's dot-size encoding is not
  # readable, so every dot is the same size and its count is printed beside it
  dot_mm <- 2 * r_for(k_max, k_max) * PX * 0.8
  # Paper panels: dot area follows the gene count on one scale for every
  # panel (1.0 mm printed at 0 genes to 4.0 mm at 40 or more), and the count
  # is still printed beside the dot. Canvas size = printed size x canvas /
  # placed width.
  hit$dot_size <- if (fixed) (1.0 + 3.0 * sqrt(pmin(hit$Count, 40) / 40)) * opt$width / row$placed_in
                  else dot_mm
  p <- ggplot(hit, aes(x = fold, y = display_wrapped)) +
    geom_point(aes(colour = neglog10q, size = dot_size)) +
    scale_size_identity() +
    # leading spaces = a gap that scales with the type, not with the x range
    geom_text(aes(label = paste0("    n = ", Count)), hjust = 0, colour = WB$muted,
              size = 9.6 * opt$font_scale / .pt) +
    scale_colour_gradientn(colours = if (opt$ramp == "purple") WB$purple_seq else PAL$go_ramp,
                           # bar runs exactly from the first to the last tick
                           limits = c(0, q_top),
                           breaks = q_breaks, labels = q_labels,
                           # two lines so the title does not set the legend width
                           name = if (opt$legend == "bottom") expression(-log[10] * "(adj. " * italic(p) * ")")
                                  else expression(atop(-log[10], "(adj. " * italic(p) * ")"))) +
    scale_x_continuous(limits = c(0, x_max), breaks = x_breaks,
                       labels = function(v) sub("\\.0$", "", sprintf("%.1f", v)),
                       # room for the "n = k" text right of the top dot: the
                       # text keeps its printed size, so a narrower plot needs
                       # a larger share
                       expand = expansion(mult = c(0.06, if (fixed) 0.30 * 3.75 / row$placed_in * 1.5 else 0.30))) +
    # a count label beside a dot at the x maximum may run past the panel edge
    coord_cartesian(clip = "off") +
    # size legend, drawn separately since the radii are identity-scaled
    guides(colour = if (opt$legend == "bottom")
             guide_colourbar(barwidth = unit(52 * opt$font_scale, "pt"), barheight = unit(6 * opt$font_scale, "pt"),
                             ticks.colour = NA, frame.colour = NA, order = 1,
                             title.position = "left", title.vjust = 1)
           else
             guide_colourbar(barwidth = unit(6, "pt"), barheight = unit(40 * opt$font_scale, "pt"),
                             ticks.colour = NA, frame.colour = NA, order = 1,
                             title.hjust = 0.5)) +
    labs(x = "Fold enrichment", y = NULL) +
    theme_minimal(base_size = 12 * opt$font_scale, base_family = WB$font) +
    theme(
      plot.background    = element_rect(fill = WB$background, colour = NA),
      panel.background   = element_rect(fill = WB$background, colour = NA),
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(colour = WB$border, linewidth = 0.4),
      axis.text.y        = element_text(colour = WB$foreground, lineheight = 0.9),
      axis.text.x        = element_text(colour = WB$muted),
      axis.title.x       = element_text(colour = WB$muted, size = 10 * opt$font_scale),
      axis.ticks         = element_blank(),
      legend.title       = element_text(colour = WB$muted, size = 9.6 * opt$font_scale),
      legend.text        = element_text(colour = WB$muted, size = 9.6 * opt$font_scale),
      legend.position    = opt$legend,
      # a bottom legend lines up with the right edge of the whole panel, so
      # a colour bar wider than a narrow plot area runs left, under the labels
      legend.location    = if (opt$legend == "bottom") "plot" else "panel",
      legend.justification = if (opt$legend == "bottom") "right" else "center",
      plot.caption       = element_text(colour = WB$muted, size = 8 * opt$font_scale, hjust = 0),
      plot.margin        = margin(8, 12, 8, 8)
    )

  # `outfile` (optional column): explicit path under outdir, no extension, may
  # contain a subdirectory. When set it fully controls the arrangement/filename
  # (e.g. "compare_DE/[Fig 2E] GO-CC ...") instead of the run/path/ontology tree.
  has_outfile <- "outfile" %in% names(row) && !is.na(row$outfile) && nzchar(str_trim(row$outfile))
  if (has_outfile) {
    base <- file.path(opt$outdir, str_trim(row$outfile))
    dir.create(dirname(base), showWarnings = FALSE, recursive = TRUE)
  } else {
    out_dir <- file.path(opt$outdir, row$run, row$path)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    base <- file.path(out_dir, sprintf("go_%s_dot%s", row$ontology,
                                       if (opt$ramp == "viridis") "" else paste0("_", opt$ramp)))
  }
  # Panels are placed at ~3.4-4.1 in; rendering at 7.0 and scaling down halves
  # every label. Pass --width/--height to render at the placed size instead.
  # row spacing is driven by type size, so the auto height has to scale with
  # font_scale or the term labels collide vertically
  # ... and never shorter than the legend column (~8 lines), or the
  # colour-bar title is cropped on 1-2 term panels
  # a wrapped label needs a taller row: rows are evenly spaced, so the
  # label with the most lines sets the row height of the whole panel
  n_lines <- max(str_count(as.character(hit$display_wrapped), "\n") + 1)
  row_h   <- max(opt$row_height, 0.11 * n_lines + 0.10)
  legend_h <- if (opt$legend == "bottom" && show_legend) 0.45 else 0
  height <- if (opt$height > 0) opt$height
            else max(2.2, (row_h * nrow(hit) + 1.4 + legend_h) * opt$font_scale,
                     8 * 1.3 * 9.6 * opt$font_scale / 72)
  if (!show_legend) p <- p + guides(colour = "none")
  if (fixed) p <- fix_widths(p, label_w, opt$width)
  ggsave(paste0(base, ".pdf"), p, width = opt$width, height = height)
  ggsave(paste0(base, ".png"), p, width = opt$width, height = height, dpi = 200)
  message("[wrote] ", base, ".{pdf,png}  (", nrow(hit), " terms)")
}

# Largest -log10(adj. p) of the plotted terms in each scale_group (optional
# targets column): those panels share one colour range.
GROUP_QMAX <- list()
if ("scale_group" %in% names(targets)) {
  for (i in seq_len(nrow(targets))) {
    g <- targets$scale_group[i]
    if (is.na(g) || !nzchar(g)) next
    csv <- ensure_csv(targets$run[i], targets$path[i], targets$ontology[i], opt$bucket, opt$cache)
    hit <- select_terms(read_csv(csv, show_col_types = FALSE), targets$terms[i], opt$top_n, "group scan")
    GROUP_QMAX[[g]] <- max(c(GROUP_QMAX[[g]], 1, -log10(hit$p.adjust)), na.rm = TRUE)
  }
}

for (i in seq_len(nrow(targets))) {
  row <- targets[i, ]
  tryCatch(render_one(row),
           error = function(e) warning(sprintf("[row %d] %s", i, conditionMessage(e))))
}
