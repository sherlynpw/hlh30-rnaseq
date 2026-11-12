#!/usr/bin/env Rscript
# Re-render the stage-PCA scatter from pipeline pca_coordinates.csv.
# - Color: Genotype_Age (6 levels: WT / hlh-30(tm1978) x Day 0 / Day 3 / Day 6)
# - Shape: dots only, no per-sample text labels
# - Ellipses: one per stage (Day 0 / Day 3 / Day 6), grey

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(ggforce)
  library(dplyr)
})

# >= 10pt in print: see _common/typography.R for FIG_SCALE / FONT_MIN.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "typography.R"))

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "palette.R"))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "enclose.R"))


opt <- parse_args(OptionParser(option_list = list(
  make_option("--coords",   type = "character"),
  make_option("--variance", type = "character"),
  make_option("--outdir",   type = "character"),
  make_option("--short-stages", action = "store_true", default = FALSE, dest = "short_stages",
              help = "Use the short stage labels L1/L4/D3/D6, matching the rest of the figure set."),
  make_option("--font-scale", type = "double", default = 1, dest = "font_scale",
              help = paste("Multiply every font size by this. Render at the",
                           "ORIGINAL canvas size and raise this until the type",
                           "is legible at the placed size -- that keeps the",
                           "layout identical instead of rescaling the canvas.")),
  make_option("--ypad",     type = "double", default = 15,
              help = paste("paper style: empty band below the data reserved for",
                           "the in-panel legend, in PC2 units. Bigger type needs",
                           "a bigger band to clear the Day 0 ellipse.")),
  make_option("--legend",   type = "character", default = "inside",
              help = paste("inside (default) or bottom. At the panel's true",
                           "placed width (95mm) no in-panel corner is free at",
                           "10pt, so 'bottom' moves both keys under the plot",
                           "and drops the now-redundant ellipse leader labels.")),
  make_option("--style",    type = "character", default = "current",
              help = paste("current (default) or paper. 'paper' reproduces the",
                           "pre-2f6ed52 look used for the published panel:",
                           "grey stage ellipses with grey leader labels, roman",
                           "genotype in the legend, Stage above Genotype.")),
  make_option("--width",    type = "double", default = 7.5,
              help = "Panel width in inches. Set to the placed width for 1:1 layout."),
  make_option("--height",   type = "double", default = 5.5)
)))

stopifnot(!is.null(opt$coords), !is.null(opt$variance), !is.null(opt$outdir))
stopifnot(opt$style %in% c("current", "paper"))
stopifnot(opt$legend %in% c("inside", "inside-left", "inside-genotype", "bottom", "inside-top"))
legend_bottom <- opt$legend == "bottom"
# one-third width: only the genotype key, untitled, inside the top of the plot
# (stages are labelled at their outlines)
legend_top    <- opt$legend == "inside-top"
legend_left   <- opt$legend == "inside-left"
# The Stage key repeats the on-plot Day 0/Day 3/Day 6 leader labels. Dropping it
# halves the legend, which is what lets the published in-panel layout survive
# type larger than its original ~4.8pt without covering a cluster.
legend_geno_only <- opt$legend == "inside-genotype"
paper_style <- opt$style == "paper"
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# --- variance line for axes ---
var_lines <- readLines(opt$variance)
pct <- function(pc) {
  m <- regmatches(var_lines, regexec(paste0("^", pc, ":\\s*([0-9.]+)%"), var_lines))
  v <- Filter(function(x) length(x) == 2, m)[[1]][2]
  as.numeric(v)
}
pc1_pct <- pct("PC1")
pc2_pct <- pct("PC2")

# --- coordinates + label remapping ---
df <- read.csv(opt$coords, stringsAsFactors = FALSE, check.names = FALSE)

stage_all  <- if (opt$short_stages) c("L1","L4","D3","D6") else c("L1","Day 0","Day 3","Day 6")
df$stage_label_chr <- dplyr::recode(df$stage,
                                    "L1" = "L1", "L4" = if (opt$short_stages) "L4" else "Day 0",
                                    "Day3" = if (opt$short_stages) "D3" else "Day 3",
                                    "Day6" = if (opt$short_stages) "D6" else "Day 6")
stage_lvls <- intersect(stage_all, unique(df$stage_label_chr))  # drop unused
df$stage_label <- factor(df$stage_label_chr, levels = stage_lvls)

# Accepts both new (WT/tm1978) and legacy (N2/hlh-30) tokens.
geno_all <- c("WT", "hlh-30(tm1978)")
df$geno_label_chr <- dplyr::recode(df$genotype,
                                   "N2" = "WT", "WT" = "WT",
                                   "hlh-30" = "hlh-30(tm1978)",
                                   "tm1978" = "hlh-30(tm1978)")
geno_lvls <- intersect(geno_all, unique(df$geno_label_chr))
df$geno_label <- factor(df$geno_label_chr, levels = geno_lvls)

# Color encodes stage, shape encodes genotype — two short legends.
palette_stage <- unname(PAL$stage)
# names follow --short-stages so the palette still keys on the display label
names(palette_stage) <- if (opt$short_stages) c("L1","L4","D3","D6") else c("L1","Day 0","Day 3","Day 6")
shape_geno <- c(
  "WT"             = 16,   # filled circle
  "hlh-30(tm1978)" = 17    # filled triangle
)

# Gene names set in italic; WT is not a gene name so it stays roman.
geno_lab_expr <- parse(text = vapply(geno_lvls, function(g) {
  if (g == "WT") sprintf('"%s"', g) else sprintf('italic("%s")', g)
}, character(1)))

# keyed off stage_all so it follows --short-stages (was hardcoded "Day 0" etc,
# which silently produced NA lookups under the short labels)
label_buffer_mm <- setNames(c(6, 6, 6, 20), stage_all)

has_l4 <- stage_all[2] %in% stage_lvls   # the L4/"Day 0" label, long or short
# L1 sits far out on PC1 and occupies the bottom-right, so the legend moves to
# the top-right and the extra bottom padding is dropped.
has_l1 <- "L1" %in% stage_lvls

# Cluster outline: the smallest ellipse that encloses every replicate of a
# stage, padded a little. Smooth like a confidence ellipse (the author asked
# for round outlines), but no outlier falls outside it, which a fitted
# covariance ellipse cannot promise. (A padded hull was used until 10/2026.)
HULL_EXPAND <- unit(3.5, "mm")
HULL_RADIUS <- unit(2, "mm")
hull_outline <- enclosing_outline(df, "PC1", "PC2", "stage_label", opt$width, opt$height,
                                  opt$font_scale, pad_in = 0.22,
                                  legend_in = if (legend_bottom) 0.9 * opt$font_scale else 0)
outline_layers <- if (paper_style) {
  # Verbatim from pca-stage.R at 7142b97, which rendered the published panel:
  # one grey stat_ellipse per stage, plus an invisible geom_mark_ellipse that
  # contributes only the bracketed leader and its grey label.
  list(
    hull_outline,
    # The leader labels duplicate the Stage key, so they are dropped when the
    # legend moves out of the panel -- at 95mm ggforce cannot fit all three and
    # silently renders only one, which reads as an error.
    if (legend_bottom) NULL else geom_mark_ellipse(
      aes(group = stage_label, label = stage_label),
      color = NA, fill = NA,
      # keep each label against its own outline; ggforce otherwise drops the
      # Day 0 label far enough away that it lands behind the in-panel legend
      label.buffer = unit(2, "mm"),
      expand = HULL_EXPAND, label.fontsize = 9.6 * opt$font_scale,
      label.colour = "grey25", con.colour = "grey55", con.size = 0.3
    )
  )
} else if (has_l1) {
  # Day 3 and Day 6 are nearly co-located, so ggforce stacks both labels above
  # the same cluster and the upper connector crosses the lower one. It offers
  # no control over label direction, so anchors are placed by hand here: Day 6
  # up-left, Day 3 right, which keeps the two leaders from ever meeting.
  cent <- aggregate(cbind(PC1, PC2) ~ stage_label, data = df, FUN = mean)
  anchor <- data.frame(
    stage_label = stage_all,
    lx          = c( 44,      3,       4,      -14),
    ly          = c(  5,     14,      -1,        7)
  )
  a <- merge(cent, anchor, by = "stage_label")
  a$col <- palette_stage[as.character(a$stage_label)]
  # Stop the leader at the ellipse edge rather than running it over the points.
  a$xend <- a$PC1 + 0.18 * (a$lx - a$PC1)
  a$yend <- a$PC2 + 0.18 * (a$ly - a$PC2)

  list(
    hull_outline,
    geom_segment(
      data = a, inherit.aes = FALSE,
      aes(x = lx, y = ly, xend = xend, yend = yend),
      color = a$col, linewidth = 0.5
    ),
    geom_label(
      data = a, inherit.aes = FALSE,
      aes(x = lx, y = ly, label = stage_label),
      color = a$col, fill = "white", label.size = NA,
      fontface = "bold", size = 3.4 * opt$font_scale
    )
  )
} else if (has_l4) {
  c(
    list(hull_outline),
    # One layer per stage so each leader + label takes that stage's colour
    # (con.colour is not an aesthetic, so it cannot vary within one layer).
    # con.cap = 0 makes the leader actually touch its ellipse - the default
    # 3mm gap left the lines ending in mid-air over the Day3/Day6 overlap.
    lapply(stage_lvls, function(s) {
      geom_mark_ellipse(
        data = df[df$stage_label == s, ],
        aes(group = stage_label, label = stage_label),
        color = NA, fill = NA,
        expand = HULL_EXPAND, label.fontsize = 9.6 * opt$font_scale,
        label.colour = palette_stage[[s]],
        label.fontface = "bold",
        # Day 3 and Day 6 ellipses overlap, so their labels compete for the same
        # spot and one gets clipped; push Day 6's further out to separate them.
        label.buffer = unit(label_buffer_mm[[s]], "mm"),
        con.colour = palette_stage[[s]], con.size = 0.5,
        con.cap = unit(0, "mm"), con.type = "straight"
      )
    })
  )
} else {
  list(
    geom_mark_ellipse(
      aes(group = stage_label, label = stage_label),
      color = "grey55", fill = NA, linewidth = 0.4,
      expand = unit(3.5, "mm"),
      label.fontsize = 9.6 * opt$font_scale, label.colour = "grey25",
      con.colour = "grey55", con.size = 0.3
    )
  )
}

p <- ggplot(df, aes(x = PC1, y = PC2)) +
  outline_layers +
  geom_point(aes(color = stage_label, shape = geno_label),
             size = 3, alpha = 0.95) +
  scale_color_manual(values = palette_stage, name = "Stage", drop = FALSE) +
  scale_shape_manual(values = shape_geno,    name = "Genotype", drop = FALSE,
                     labels = geno_lab_expr) +   # gene name italic in every PCA panel
  guides(shape = guide_legend(order = if (paper_style) 2 else 1),
         # with L1, every stage is labelled at its outline, so a bottom legend
         # (one-third width) keeps only Genotype
         color = if (legend_geno_only || legend_top || (legend_bottom && has_l1)) "none"
                 else guide_legend(order = if (paper_style) 1 else 2)) +
  # paper: the bottom band holds the in-panel legend. The original reserved 15
  # for 9pt type; at 10pt the legend is taller and rode up into the Day 0
  # ellipse and its leader label, so the band is widened to 22.
  coord_cartesian(ylim = c(min(df$PC2) - if (legend_bottom) 3 else if (paper_style) opt$ypad else if (has_l1) 3 else 15,
                           max(df$PC2) + if (legend_top) 20 else if (has_l1) 10 else 2), clip = "off") +
  labs(
    x = sprintf("PC1 (%.2f%%)", pc1_pct),
    y = sprintf("PC2 (%.2f%%)", pc2_pct)
  ) +
  theme_bw(base_size = 12 * opt$font_scale) +
  theme(
    panel.grid.minor    = element_blank(),
    # paper: bottom-left, not the original's bottom-right -- at 10pt the legend
    # covers Day 0's leader label, which ggforce auto-places into that corner.
    legend.position      = if (legend_top) c(0.5, 0.99) else if (legend_bottom) "bottom"
                           else if (legend_left) c(0.02, 0.02)
                           else if (paper_style) c(0.98, 0.02) else if (has_l1) c(0.98, 0.98) else if (has_l4) c(0.98, 0.02) else c(0.02, 0.02),
    legend.justification = if (legend_top) c(0.5, 1) else if (legend_bottom) "center"
                           else if (legend_left) c(0, 0)
                           else if (paper_style) c(1, 0)       else if (has_l1) c(1, 1)       else if (has_l4) c(1, 0)       else c(0, 0),
    # legend title at the same 10pt as its keys
    legend.title        = if (legend_top) element_blank() else element_text(size = 9.6 * opt$font_scale),
    legend.background   = element_rect(fill = alpha("white", 0.85), color = NA),
    legend.box          = "vertical",
    legend.direction    = if (legend_top) "horizontal" else "vertical",
    legend.spacing.y    = unit(2, "pt"),
    legend.key          = element_blank()
  )

out_pdf <- file.path(opt$outdir, "pca_stage.pdf")
out_png <- file.path(opt$outdir, "pca_stage.png")
ggsave(out_pdf, p, width = opt$width, height = opt$height)
ggsave(out_png, p, width = opt$width, height = opt$height, dpi = 200)

message("wrote: ", out_pdf)
message("wrote: ", out_png)
