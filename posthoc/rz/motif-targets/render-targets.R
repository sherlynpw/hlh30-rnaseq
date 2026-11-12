#!/usr/bin/env Rscript
# Two figures for "what % of each DESeq2 DEG list carries an HLH-30 motif in
# the promoter", from motif_targets_summary.tsv (DEG family only), one column
# per motif definition (motif-targets.R):
#   dumbbell - background% vs list%, grouped by timepoint, all/up/down stacked
#   dotgrid  - timepoint x direction matrix, colour = fold, size = n
# Each motif has its own background line; stars are BH-adjusted Fisher tests.

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(ggplot2); library(stringr); library(ggtext); library(patchwork)
})

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "palette.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--summary", type = "character", default = "out/hlh30-all-ruv/motif_targets_summary.tsv"),
  make_option("--font-scale", type = "double", default = 1, dest = "font_scale",
              help = "Multiply every font size; render at the ORIGINAL canvas."),
  make_option("--width",   type = "double", default = 14.18, help = "canvas width (in)"),
  make_option("--height",  type = "double", default = 9, help = "canvas height (in)"),
  make_option("--geom-scale", type = "double", default = 2, dest = "geom_scale",
              help = "multiply point and line sizes; canvas / placed width keeps them at print size"),
  make_option("--outdir",  type = "character", default = "out/hlh30-all-ruv"),
  make_option("--accent",  type = "character", default = "#7B3294"),
  make_option("--wt-label", type = "character", default = "WT", dest = "wt_label",
              help = "wild-type name in the group label (WT or N2)"),
  make_option("--notes", action = "store_true", default = FALSE,
              help = "draw the title, subtitle and caption (off for the paper panel; the caption text goes to *_legend.txt)"),
  make_option("--motifs", type = "character", default = "hlh30",
              help = "motif sets to draw, comma-separated, of any,hlh30,chip (fig S5C: hlh30 only)")
)))
ACC <- opt$accent
GREY <- "#9A9AA0"
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

res <- read_tsv(opt$summary, show_col_types = FALSE) %>% filter(family == "deg")
MOT <- c(any = "A  Any HLH-30 motif",
         hlh30 = "B  Canonical E-box",
         chip = "C  ChIP-seq motif")
# the motif each panel scans for, shown under its title
MONO <- function(x) sprintf("<span style='font-family:mono'><b>%s</b></span>", x)
SEQ <- c(any = "a site of B or C<br>&nbsp;",
         hlh30 = paste0(MONO("TCACGTGA"), " (JASPAR MA1449.1)"),
         chip = paste0(MONO("CACGTGRT"), "<br>de novo, L4 ChIP peaks"))
names(SEQ) <- MOT[names(SEQ)]
SHOW <- strsplit(opt$motifs, ",")[[1]]
MOT <- MOT[SHOW]
res <- res %>% filter(motif %in% SHOW) %>% mutate(motif = factor(MOT[motif], levels = MOT))
# BH over the tests drawn, not over every set the scan wrote
res$padj <- p.adjust(res$p, "BH")
ONE <- length(SHOW) == 1
bgs <- res %>% distinct(motif, bg_pct)
res$star <- ifelse(is.na(res$padj), "", ifelse(res$padj < 0.001, "***",
            ifelse(res$padj < 0.01, "**", ifelse(res$padj < 0.05, "*", ""))))

# timepoint display + order (L4 is reported as "Day 0")
stage_disp <- c(L1 = "L1", `Day 0` = "L4", `Day 3` = "D3", `Day 6` = "D6")
res <- res %>%
  mutate(stage = factor(stage_disp[group], levels = c("L1","L4","D3","D6")),
         dir_lab = factor(c(all = "all DEGs", up = "up DEGs", down = "down DEGs")[direction],
                          levels = c("all DEGs","up DEGs","down DEGs")))

# per-timepoint strip label: "WT vs. hlh-30(tm1978)" over the timepoint
# (the comparison itself is named once, in the figure subtitle)
strip_md <- setNames(sprintf("<b>%s</b>", levels(res$stage)), levels(res$stage))

cap <- paste(
  "Promoter = -2000/+200 bp of the TSS.",
  if (ONE && SHOW == "hlh30") "Sites are matches to the JASPAR MA1449.1 matrix at p < 1e-4, so they include close variants such as TCATGTGA."
  else "Motif sites called at p < 1e-4. R = A or G.",
  "Dashed line (BG) = background, the share of expressed genes with the motif.",
  sprintf("Stars = Fisher's exact test vs background, BH-adjusted across all %d tests.", nrow(res)))

base_theme <- theme_minimal(base_size = 11 * opt$font_scale) +
  theme(plot.title.position = "plot", plot.caption.position = "plot",
        plot.margin = margin(4, 16, 4, 4),
panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12 * opt$font_scale),
        # the subtitle's bottom margin is the gap the in-panel "background %"
        # label (drawn above the first row, clip off) needs to clear it
        plot.subtitle = element_text(size = 9 * opt$font_scale, colour = "grey40",
                                     margin = margin(b = 16 * opt$font_scale)),
        plot.caption = element_text(size = 8.8 * opt$font_scale, colour = "grey45", hjust = 0, lineheight = 1.1))
GS <- opt$geom_scale

# ================================================================= dumbbell ===
# ggplot draws the first y level at the bottom, so order down/up/all to get
# all on top, up middle, down bottom within each timepoint.
du <- res %>% mutate(dir_lab = factor(dir_lab, levels = c("down DEGs","up DEGs","all DEGs")))
# one panel per motif set in a row; timepoint and list labels on the left only
panel <- function(m, left) {
  d  <- du %>% filter(motif == m)
  bv <- d$bg_pct[1]
  ggplot(d) +
    geom_vline(xintercept = bv, linetype = "dashed", colour = "grey55", linewidth = 0.4 * GS) +
    geom_segment(aes(x = bg_pct, xend = pct, y = dir_lab, yend = dir_lab),
                 colour = "grey80", linewidth = 0.7 * GS) +
    geom_point(aes(bg_pct, dir_lab), colour = GREY, size = 1.5 * GS) +
    geom_point(aes(pct, dir_lab), colour = ACC, size = 2.3 * GS) +
    geom_text(aes(x = pmax(pct, bg_pct), y = dir_lab, label = sprintf("  %.0f%%%s", pct, star)),
              hjust = 0, size = 3.1 * opt$font_scale, colour = "grey20") +
    # the background value, once, just above the dashed line
    geom_text(data = d %>% filter(stage == levels(stage)[1]) %>% slice(1),
              aes(x = bv, y = Inf, label = sprintf("BG: %.1f%%", bv)),
              vjust = -0.4, size = 3.1 * opt$font_scale, colour = "grey40") +
    facet_grid(stage ~ ., switch = "y", labeller = labeller(stage = strip_md)) +
    coord_cartesian(clip = "off") +
    scale_x_continuous(expand = expansion(mult = c(0.06, if (ONE) 0.3 else 1.45)),
                       breaks = scales::breaks_pretty(n = 3)) +
    labs(x = if (!left) " " else if (ONE && SHOW == "hlh30") "promoters with a TCACGTGA-like site (%)"
             else "promoters with the motif (%)",
         y = NULL, title = if (ONE) NULL else m,
         subtitle = SEQ[[as.character(m)]]) +
    base_theme +
    theme(plot.title = element_text(face = "bold", size = 10 * opt$font_scale),
          plot.subtitle = element_markdown(size = 9 * opt$font_scale, colour = "grey40", lineheight = 1.2,
                                           margin = margin(b = 14 * opt$font_scale)),
          axis.title.x = element_text(hjust = 0),
          axis.text.y = if (left) element_text(colour = "grey20") else element_blank(),
          panel.grid.major.y = element_line(colour = "grey93"),
          panel.spacing.y = unit(4 * opt$font_scale, "pt"),
          strip.placement = "outside",
          strip.text.y.left = if (left) element_markdown(angle = 0, hjust = 1, lineheight = 1.35,
                                                         size = 10 * opt$font_scale, colour = "grey15",
                                                         margin = margin(r = 4 * opt$font_scale))
                              else element_blank(),
          plot.margin = margin(4, 10, 4, 2))
}
ps  <- Map(panel, levels(du$motif), seq_along(MOT) == 1)
# the paper panel is the plot alone; its notes go to the figure legend
writeLines(c("Grey dot = expressed-gene background, purple dot = the DEG list.", stringr::str_wrap(cap, 100)),
           file.path(opt$outdir, "targets_dumbbell_legend.txt"))
pdb <- wrap_plots(ps, nrow = 1)
if (opt$notes) pdb <- pdb +
  plot_annotation(title = if (ONE && SHOW == "hlh30") "Canonical HLH-30 E-box in DEG promoters" else "HLH-30 motifs in DEG promoters",
                  subtitle = sprintf("DEGs of %s vs. <i>hlh-30(tm1978)</i>.%sGrey dot = expressed-gene background, purple dot = the DEG list.",
                                     opt$wt_label, if (ONE) "<br>" else " "),
                  caption = stringr::str_wrap(cap, round(105 * opt$width / 14.18)),
                  theme = base_theme + theme(plot.subtitle = element_markdown(size = 9 * opt$font_scale,
                                                                              colour = "grey40",
                                                                              margin = margin(b = (if (ONE) 2 else 16) * opt$font_scale))))
ggsave(file.path(opt$outdir, "targets_dumbbell.pdf"), pdb, width = opt$width, height = opt$height)
ggsave(file.path(opt$outdir, "targets_dumbbell.png"), pdb, width = opt$width, height = opt$height, dpi = 200)

# ================================================================== dot grid ===
gd <- res %>%
  mutate(direction = factor(c(all="all DEGs", up="up DEGs", down="down DEGs")[direction],
                            levels = c("all DEGs","up DEGs","down DEGs")),
         stage = factor(stage, levels = rev(c("L1","L4","D3","D6"))))  # L1 on top
flim <- max(abs(gd$fold - 1), na.rm = TRUE)
flim <- min(flim, 1)
pdg <- ggplot(gd, aes(direction, stage)) +
  geom_point(aes(fill = fold, size = n), shape = 21, colour = "grey55", stroke = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%%s", pct, star)), size = 2.7 * opt$font_scale, colour = "grey15",
            nudge_y = 0.32) +
  facet_wrap(~ motif, nrow = 1) +
  scale_fill_gradient2(low = PAL$div[["low"]], mid = "white", high = ACC, midpoint = 1,
                       limits = c(1 - flim, 1 + flim), oob = scales::squish,
                       name = "fold vs\nbackground") +
  scale_size(range = c(3, 14), name = "genes in list",
             breaks = c(50, 150, 300, 500)) +
  guides(size = guide_legend(override.aes = list(fill = "grey70", colour = "grey40")),
         fill = guide_colourbar(barwidth = unit(7, "pt"), barheight = unit(60, "pt"))) +
  labs(x = NULL, y = NULL,
       title = "HLH-30 motifs across DEG lists",
       subtitle = "Fill = fold vs background (1 = background; purple = enriched, blue = depleted); size = list size.",
       caption = cap) +
  base_theme +
  theme(panel.grid.major = element_line(colour = "grey94"),
        axis.text = element_text(size = 10 * opt$font_scale, colour = "grey20"))
ggsave(file.path(opt$outdir, "targets_dotgrid.pdf"), pdg, width = 19, height = 9)
ggsave(file.path(opt$outdir, "targets_dotgrid.png"), pdg, width = 19, height = 9, dpi = 200)

message("wrote targets_{dumbbell,dotgrid}.{pdf,png} to ", opt$outdir)
