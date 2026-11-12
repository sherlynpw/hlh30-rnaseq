#!/usr/bin/env Rscript
# PCA scatter from pipeline pca_coordinates.csv.
# Ellipses grouped by genotype x stage (4 ellipses for no_l4 runs).
# Color encodes stage, shape encodes genotype (same legend style as pca-stage).

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
              help = "Multiply every font size. Render at the ORIGINAL canvas and raise this."),
  make_option("--legend",   type = "character", default = "inside",
              help = paste("inside (default) or bottom. At the placed width",
                           "(95mm) the 'WT | Day 3' ellipse labels cannot fit",
                           "at 10pt; 'bottom' moves the keys under the plot and",
                           "drops those labels, which the keys already encode.")),
  make_option("--group-shape", type = "character", default = "hull", dest = "group_shape",
              help = paste("Outline per genotype x stage: ellipse (stat_ellipse level 0.85,",
                           "as fig 1B), mark (smallest enclosing ellipse) or hull.")),
  make_option("--width",    type = "double", default = 7.5,
              help = "Panel width in inches. Set to the placed width for 1:1 layout."),
  make_option("--height",   type = "double", default = 5.5),
  make_option("--color-by", type = "character", default = "stage", dest = "color_by",
              help = "stage (default) or cohort: the adult growth temperature (fig S1C)"),
  make_option("--outname",  type = "character", default = "pca_genotype")
)))
stopifnot(opt$color_by %in% c("stage", "cohort"))
stopifnot(opt$legend %in% c("inside", "inside-h", "bottom"))
legend_bottom <- opt$legend == "bottom"
# Keys inside the panel, laid out as two horizontal rows in the empty
# bottom-left corner. A vertical block is far too tall to sit inside at this
# panel size, and putting it under the plot wastes ~40% of the panel.
legend_inside_h <- opt$legend == "inside-h"

stopifnot(!is.null(opt$coords), !is.null(opt$variance), !is.null(opt$outdir))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

var_lines <- readLines(opt$variance)
pct <- function(pc) {
  m <- regmatches(var_lines, regexec(paste0("^", pc, ":\\s*([0-9.]+)%"), var_lines))
  v <- Filter(function(x) length(x) == 2, m)[[1]][2]
  as.numeric(v)
}
pc1_pct <- pct("PC1")
pc2_pct <- pct("PC2")

df <- read.csv(opt$coords, stringsAsFactors = FALSE, check.names = FALSE)

stage_chr  <- dplyr::recode(df$stage, "L1" = "L1",
                            "L4" = if (opt$short_stages) "L4" else "Day 0",
                            "Day3" = if (opt$short_stages) "D3" else "Day 3",
                            "Day6" = if (opt$short_stages) "D6" else "Day 6")
stage_lvls <- intersect(if (opt$short_stages) c("L1","L4","D3","D6") else c("L1","Day 0","Day 3","Day 6"), unique(stage_chr))
df$stage_label <- factor(stage_chr, levels = stage_lvls)

# Accepts both new (WT/tm1978) and legacy (N2/hlh-30) tokens.
geno_chr  <- dplyr::recode(df$genotype, "N2" = "WT", "WT" = "WT",
                           "hlh-30" = "hlh-30(tm1978)", "tm1978" = "hlh-30(tm1978)")
geno_lvls <- intersect(c("WT", "hlh-30(tm1978)"), unique(geno_chr))
df$geno_label <- factor(geno_chr, levels = geno_lvls)
df$geno_stage <- interaction(df$geno_label, df$stage_label, drop = TRUE, sep = " | ")

has_l1 <- "L1" %in% stage_lvls

palette_stage <- unname(PAL$stage)
# names follow --short-stages so the palette still keys on the display label
names(palette_stage) <- if (opt$short_stages) c("L1","L4","D3","D6") else c("L1","Day 0","Day 3","Day 6")
# fig S1C: the same points coloured by the adult cohort (20 C or room
# temperature, 22 C), cool and warm Okabe-Ito colours that no stage uses
df$cohort_label <- factor(dplyr::recode(df$temperature, "20C" = "20 \u00b0C", "RT" = "22 \u00b0C"),
                          levels = c("20 \u00b0C", "22 \u00b0C"))
palette_cohort <- c("20 \u00b0C" = PAL$okabe[5], "22 \u00b0C" = PAL$okabe[6])
colour_var   <- if (opt$color_by == "cohort") "cohort_label" else "stage_label"
colour_vals  <- if (opt$color_by == "cohort") palette_cohort else palette_stage
colour_title <- if (opt$color_by == "cohort") "Cohort" else "Stage"
shape_geno <- c(
  "WT"             = 16,
  "hlh-30(tm1978)" = 17
)

# Gene names set in italic; WT is not a gene name so it stays roman.
geno_lab_expr <- parse(text = vapply(geno_lvls, function(g) {
  if (g == "WT") sprintf('"%s"', g) else sprintf('italic("%s")', g)
}, character(1)))

stopifnot(opt$group_shape %in% c("ellipse", "mark", "hull"))
ellipse_layer <- if (opt$group_shape == "ellipse") {
  # Same fit as pca-stage (fig 1B): the smallest enclosing ellipse stretches
  # to any one outlying replicate, a covariance ellipse follows the cloud.
  stat_ellipse(aes(group = geno_stage), type = "t", level = 0.85,
               color = "grey55", linewidth = 0.4)
} else if (opt$group_shape == "hull") {
  # A padded, rounded hull: follows the replicates instead of fitting a shape
  geom_mark_hull(aes(group = geno_stage), concavity = 5,
                 color = "grey55", fill = NA, linewidth = 0.4,
                 expand = unit(3.5, "mm"), radius = unit(2, "mm"))
} else if (legend_bottom || legend_inside_h) {
  # Unlabelled outlines: colour (stage) and shape (genotype) in the bottom keys
  # already identify every ellipse, and the "WT | Day 3" labels are far too
  # long to place at 10pt in a 95mm panel. Our own enclosing ellipse: ggforce's
  # left points outside here.
  enclosing_outline(df, "PC1", "PC2", "geno_stage", opt$width, opt$height,
                    opt$font_scale, pad_in = 0.22)
} else if (nlevels(droplevels(df$geno_stage)) <= 4) {
  geom_mark_ellipse(
    aes(group = geno_stage, label = as.character(geno_stage)),
    color = "grey55", fill = NA, linewidth = 0.4,
    expand = unit(2.5, "mm"), label.fontsize = 9.6 * opt$font_scale,
    label.colour = "grey25", con.colour = "grey55", con.size = 0.3
  )
} else {
  geom_mark_ellipse(
    aes(group = geno_stage),
    color = "grey55", fill = NA, linewidth = 0.4,
    expand = unit(2.5, "mm")
  )
}

p <- ggplot(df, aes(x = PC1, y = PC2)) +
  # One outline per (genotype, stage). geom_mark_ellipse rather than
  # stat_ellipse: the n=3 groups (L1, Day 0) are below stat_ellipse's minimum.
  # Leader labels only for <=4 groups; beyond that they collide and the
  # colour/shape legends already identify every ellipse.
  ellipse_layer +
  geom_point(aes(color = .data[[colour_var]], shape = geno_label),
             size = 3, alpha = 0.95) +
  scale_color_manual(values = colour_vals, name = colour_title, drop = FALSE) +
  scale_shape_manual(values = shape_geno,    name = "Genotype", drop = FALSE,
                     labels = geno_lab_expr) +
  # Stage above Genotype, matching the published panel and pca-stage.
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  coord_cartesian(ylim = c(min(df$PC2) - if (legend_bottom) 3 else if (legend_inside_h) 13 else if (has_l1) 3 else 15,
                           max(df$PC2) + if (has_l1) 6 else 2), clip = "off") +
  labs(
    x = sprintf("PC1 (%.2f%%)", pc1_pct),
    y = sprintf("PC2 (%.2f%%)", pc2_pct)
  ) +
  theme_bw(base_size = 12 * opt$font_scale) +
  theme(
    panel.grid.minor    = element_blank(),
    legend.position     = if (legend_bottom) "bottom" else if (legend_inside_h) c(0.01, 0.01) else if (has_l1) c(0.98, 0.98) else c(0.98, 0.02),
    legend.justification = if (legend_bottom) "center" else if (legend_inside_h) c(0, 0) else if (has_l1) c(1, 1)      else c(1, 0),
    legend.direction    = if (legend_inside_h) "horizontal" else "vertical",
    # legend title at the same 10pt as its keys
    legend.title        = element_text(size = 9.6 * opt$font_scale),
    legend.background   = element_rect(fill = alpha("white", 0.85), color = NA),
    legend.box          = "vertical",
    legend.spacing.y    = unit(2, "pt"),
    legend.key          = element_blank()
  )

out_pdf <- file.path(opt$outdir, paste0(opt$outname, ".pdf"))
out_png <- file.path(opt$outdir, paste0(opt$outname, ".png"))
ggsave(out_pdf, p, width = opt$width, height = opt$height)
ggsave(out_png, p, width = opt$width, height = opt$height, dpi = 200)

message("wrote: ", out_pdf)
message("wrote: ", out_png)
