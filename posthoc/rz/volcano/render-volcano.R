#!/usr/bin/env Rscript
# Volcano plots: mutant vs WT at each stage.
# x = log2FoldChange, y = -log10(padj), coloured by DEG direction.
# Genes named in --highlight are drawn on top in a third colour and labelled.

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggrepel)
  library(patchwork)
  library(ggtext)
})

# >= 10pt in print: see _common/typography.R for FIG_SCALE / FONT_MIN.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "typography.R"))

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "palette.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--cache",     type = "character"),
  make_option("--stages",    type = "character",
              help = "space-separated stage names, e.g. 'L1 L4 Day3 Day6'"),
  make_option("--mut",       type = "character", default = "tm1978",
              help = "mutant token used in the contrast filename"),
  make_option("--wt",        type = "character", default = "WT",
              help = "WT token used in the contrast filename"),
  make_option("--highlight", type = "character", default = "",
              help = "comma-separated gene names to highlight and label"),
  make_option("--padj-cut",  type = "double", default = 0.05, dest = "padj_cut"),
  make_option("--lfc-cut",   type = "double", default = 1,    dest = "lfc_cut"),
  make_option("--y-max",     type = "double", default = 30,   dest = "y_max",
              help = "Cap -log10(padj) at this value; outliers drawn as triangles. 0 = no cap, every dot at its real position."),
  make_option("--label-top", type = "integer", default = 0,   dest = "label_top",
              help = "Label this many top genes (combined rank). 0 = no labels."),
  make_option("--outdir",    type = "character"),
  make_option("--short-stages", action = "store_true", default = FALSE, dest = "short_stages",
              help = "Use the short stage labels L1/L4/D3/D6, matching the rest of the figure set."),
  make_option("--short-counts", action = "store_true", default = FALSE, dest = "short_counts",
              help = "Label the DEG counts 'up: N' / 'down: N' without the mutant name, for narrow panels."),
  make_option("--count-size", type = "double", default = 3.4, dest = "count_size",
              help = "DEG count annotation size in mm (3.4 = same as 10pt axis text at the paper scale)."),
  make_option("--one-line-counts", action = "store_true", default = FALSE, dest = "one_line_counts",
              help = "Each count on ONE full-width line: 'N up in <mut>' right-aligned on row 1, 'N down in <mut>' left-aligned on row 2."),
  make_option("--plain-counts", action = "store_true", default = FALSE, dest = "plain_counts",
              help = "Counts as 'N up' / 'N down' on one line each, no mutant name."),
  make_option("--counts-in-title", action = "store_true", default = FALSE, dest = "counts_in_title",
              help = "Put the counts beside the stage name above the plot, so the plot area is free for the gene labels."),
  make_option("--stage-inside", action = "store_true", default = FALSE, dest = "stage_inside",
              help = paste("Stage name inside the plot, top left; 'N down' and 'N up' above the",
                           "plot at its left and right edges. Implies --plain-counts.")),
  make_option("--label-gap", type = "double", default = 0.3, dest = "label_gap",
              help = "Gap (in lines) between a highlighted dot and its label (ggrepel point.padding)."),
  make_option("--stack-counts", action = "store_true", default = FALSE, dest = "stack_counts",
              help = "Put the 'down' block under the 'up' block instead of beside it (narrow panels)."),
  make_option("--font-scale", type = "double", default = 1, dest = "font_scale",
              help = "Multiply every font size. Render at the ORIGINAL canvas and raise this."),
  make_option("--width",     type = "double", default = 5.0,
              help = "Panel width in inches. Set to the placed width for 1:1 layout."),
  make_option("--height",    type = "double", default = 4.6)
)))

stopifnot(!is.null(opt$cache), !is.null(opt$stages), !is.null(opt$outdir))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

stages <- strsplit(trimws(opt$stages), "\\s+")[[1]]
highlight <- trimws(strsplit(opt$highlight, ",")[[1]])
highlight <- highlight[nzchar(highlight)]

up_col   <- PAL$up      # up in the mutant
down_col <- PAL$down    # down in the mutant
ns_col   <- PAL$ns
hl_col   <- PAL$hlh30   # highlighted HLH-30 target genes

# Display name for the mutant, used in the on-plot counts.
mut_display <- if (opt$mut %in% c("tm1978", "hlh-30")) "hlh-30(tm1978)" else opt$mut

stage_label <- function(s) {
  m <- if (opt$short_stages) {
    c(L1 = "L1", L4 = "L4", Day3 = "D3", Day6 = "D6")
  } else {
    c(L1 = "L1", L4 = "Day 0", Day3 = "Day 3", Day6 = "Day 6")
  }
  if (s %in% names(m)) m[[s]] else s
}

build_one <- function(stage) {
  csv <- file.path(opt$cache, sprintf("%s_%s_vs_%s_%s.csv", opt$mut, stage, opt$wt, stage))
  if (!file.exists(csv)) {
    warning("missing contrast CSV: ", csv); return(NULL)
  }
  d <- read_csv(csv, show_col_types = FALSE) %>%
    filter(!is.na(padj), !is.na(log2FoldChange))
  # no cap: the axis runs to the largest value, rounded up to a multiple of 10
  y_cap <- if (opt$y_max > 0) opt$y_max else 10 * ceiling(max(-log10(d$padj)) / 10)
  d <- d %>%
    mutate(
      neglog10 = pmin(-log10(padj), y_cap),
      clipped  = -log10(padj) > y_cap,
      sig      = padj < opt$padj_cut & abs(log2FoldChange) >= opt$lfc_cut,
      dir      = case_when(
        sig & log2FoldChange >=  opt$lfc_cut ~ "up",
        sig & log2FoldChange <= -opt$lfc_cut ~ "down",
        TRUE                                 ~ "ns"
      )
    )

  n_up   <- sum(d$dir == "up")
  n_down <- sum(d$dir == "down")

  hl <- d %>% filter(gene_name %in% highlight)

  x_lim <- max(3, ceiling(max(abs(d$log2FoldChange), na.rm = TRUE)))
  # Headroom above the y cap for the two count annotations: "N up in" / mutant
  # over the red half (right-aligned), "N down in" / mutant over the blue half
  # (left-aligned), two rows each so neither line crosses the centre at 80mm.
  # Row height grows with the type.
  # One text line in data units, times 1.4 for leading, so rows never touch
  # at any scale. The y-range depends on the rows and the rows on the
  # y-range, so iterate; the plot area is the canvas minus ~1.4 in of title,
  # axis text and axis title at paper size.
  line_pt <- opt$count_size * opt$font_scale * .pt
  n_rows  <- if (opt$short_counts || opt$one_line_counts || opt$plain_counts) 1 else 2
  head_rows <- if (opt$stack_counts || opt$one_line_counts) 2 * n_rows + 0.4 else n_rows + 0.2
  if (opt$stage_inside) { n_rows <- 1; head_rows <- 1.2 }   # one row for the stage name
  y_top <- y_cap
  if (!opt$counts_in_title || opt$stage_inside) for (k in 1:5) {
    row_h <- 1.4 * line_pt * y_top / ((opt$height - 1.4) * 72)
    y_top <- y_cap + head_rows * row_h
  }
  count_rows <- function(dir, n, x, hjust, col, block = 0) {
    lines <- if (opt$plain_counts) sprintf('"%d %s"', n, dir)
             else if (opt$short_counts) sprintf('"%s: %d"', dir, n)
             else if (opt$one_line_counts) sprintf('"%d %s in "*italic("%s")', n, dir, mut_display)
             else c(sprintf('"%d %s in"', n, dir), sprintf('italic("%s")', mut_display))
    lapply(seq_along(lines), function(i)
      annotate("text", x = x, y = y_top - (block * (n_rows + 0.2) + i - 1) * row_h, hjust = hjust, vjust = 1,
               parse = TRUE, label = lines[i], color = col, fontface = "bold",
               size = opt$count_size * opt$font_scale))
  }

  p <- ggplot(d, aes(x = log2FoldChange, y = neglog10)) +
    geom_vline(xintercept = c(-opt$lfc_cut, opt$lfc_cut),
               linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_hline(yintercept = -log10(opt$padj_cut),
               linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_point(aes(color = dir, shape = clipped), size = 1.1, alpha = 0.75) +
    scale_color_manual(values = c(up = up_col, down = down_col, ns = ns_col),
                       guide = "none") +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    (if (!opt$counts_in_title && !opt$stage_inside) count_rows("up", n_up, x_lim, 1, up_col)) +
    (if (!opt$counts_in_title && !opt$stage_inside) count_rows("down", n_down, -x_lim, 0, down_col, block = if (opt$stack_counts || opt$one_line_counts) 1 else 0)) +
    # stage inside, top left; counts just above the plot, flush with its edges
    (if (opt$stage_inside) list(
      annotate("text", x = -x_lim, y = y_top, hjust = 0, vjust = 1, label = stage_label(stage),
               fontface = "bold", colour = "grey15", size = 3.4 * opt$font_scale * 1.2),
      # the panel edges are 5% of the range beyond the limits (default expansion)
      annotate("text", x = -1.1 * x_lim, y = 1.05 * y_top, hjust = 0, vjust = -0.5,
               label = sprintf("%d down", n_down), fontface = "bold", colour = down_col,
               size = 3.4 * opt$font_scale * 1.2),
      annotate("text", x = 1.1 * x_lim, y = 1.05 * y_top, hjust = 1, vjust = -0.5,
               label = sprintf("%d up", n_up), fontface = "bold", colour = up_col,
               size = 3.4 * opt$font_scale * 1.2))) +
    coord_cartesian(xlim = c(-x_lim, x_lim), ylim = c(0, y_top), clip = if (opt$stage_inside) "off" else "on") +
    labs(
      x = expression(log[2] ~ "fold change"),
      y = expression(-log[10] ~ "adj. " * italic(p)),
      subtitle = if (opt$stage_inside) " "
        else if (opt$counts_in_title)
        sprintf("%s\u2003<span style='color:%s'>%d down</span>\u2002<span style='color:%s'>%d up</span>",
                stage_label(stage), down_col, n_down, up_col, n_up)
        else stage_label(stage),
      # explain the capped points, only when a non-highlighted gene is capped
      # (highlighted genes are drawn as labelled dots whatever their p)
      caption = if (any(d$clipped & !(d$gene_name %in% highlight)))
        bquote("Triangles: adj." ~ italic(p) < 10^-.(y_cap))
    ) +
    theme_bw(base_size = 12 * opt$font_scale) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92"),
      # with --stage-inside the (blank) subtitle only holds the room for the counts
      plot.subtitle    = element_markdown(face = "bold", color = "grey15",
                                          size = (if (opt$stage_inside) 16 else 12) * opt$font_scale),
      # 9pt in print: axis text is 9.6 x font_scale -> 10pt
      plot.caption     = element_text(size = 9.1 * opt$font_scale, colour = "grey40", hjust = 0),
      plot.caption.position = "plot",
      plot.margin = margin(if (opt$stage_inside) 8 * opt$font_scale else 5.5, 5.5, 5.5, 5.5)
    )

  if (nrow(hl) > 0) {
    # Gene labels: each starts up and outward of its dot (up-left for down
    # genes, up-right for up genes), where a volcano is empty. The DEG dots
    # join with empty labels, so the gene labels repel from them too.
    lab_df <- bind_rows(
      hl %>% mutate(lab = sprintf('italic("%s")', gene_name), seg = hl_col,
                    nx = sign(log2FoldChange) * 0.12 * x_lim,
                    ny = pmin(0.3 * y_cap, 0.92 * y_top - neglog10)),
      d %>% filter(dir != "ns", !gene_name %in% highlight) %>%
        mutate(lab = '""', seg = NA, nx = 0, ny = 0),
      # invisible points over the stage name, so no label lands on it
      if (opt$stage_inside)
        tibble(log2FoldChange = -x_lim + c(0, 0.12, 0.24) * x_lim, neglog10 = 0.97 * y_top,
               lab = '""', seg = NA, nx = 0, ny = 0))
    p <- p +
      scale_discrete_identity(aesthetics = "segment.colour") +
      geom_point(data = hl, shape = 21, size = 2.9, stroke = 0.6,
                 fill = hl_col, colour = "grey15") +
      geom_text_repel(
        data = lab_df, parse = TRUE, aes(label = lab, segment.colour = seg),
        size = 3.4 * opt$font_scale, colour = hl_col, fontface = "bold",
        min.segment.length = 0, segment.size = 0.4,
        box.padding = 1.0, point.padding = opt$label_gap, max.overlaps = Inf, seed = 1,
        # weak pull back to the dot, so a label keeps its nudged spot
        force = 2, force_pull = 0.1, max.iter = 1e5, max.time = 3,
        nudge_x = lab_df$nx, nudge_y = lab_df$ny,
        # keep gene labels out of the stage-name row
        ylim = if (opt$stage_inside) c(0, y_top) else c(-Inf, Inf)
      )
  }

  if (opt$label_top > 0) {
    top <- d %>%
      filter(sig, !gene_name %in% highlight) %>%
      mutate(score = neglog10 * abs(log2FoldChange)) %>%
      arrange(desc(score)) %>%
      slice_head(n = opt$label_top)
    if (nrow(top) > 0) {
      p <- p + geom_text_repel(
        data = top, aes(label = gene_name),
        size = 3.4 * opt$font_scale, color = "grey15", min.segment.length = 0,
        box.padding = 0.4, max.overlaps = Inf, seed = 1
      )
    }
  }

  list(plot = p, n_up = n_up, n_down = n_down, stage = stage, n_hl = nrow(hl))
}

results <- lapply(stages, build_one)
results <- results[!sapply(results, is.null)]

for (r in results) {
  base <- file.path(opt$outdir, sprintf("volcano_%s", r$stage))
  ggsave(paste0(base, ".pdf"), r$plot, width = opt$width, height = opt$height)
  ggsave(paste0(base, ".png"), r$plot, width = opt$width, height = opt$height, dpi = 200)
  message(sprintf("[wrote] %s.{pdf,png}  up=%d  down=%d  highlighted=%d",
                  base, r$n_up, r$n_down, r$n_hl))
}

if (length(results) >= 2) {
  combo <- Reduce(`+`, lapply(results, `[[`, "plot")) +
    plot_layout(nrow = 1)
  base <- file.path(opt$outdir, "volcano_combined")
  ggsave(paste0(base, ".pdf"), combo,
         width = opt$width * length(results), height = opt$height)
  ggsave(paste0(base, ".png"), combo,
         width = opt$width * length(results), height = opt$height, dpi = 200)
  message("[wrote] ", base, ".{pdf,png}  (combined ", length(results), " stages)")
}
