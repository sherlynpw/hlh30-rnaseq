#!/usr/bin/env Rscript
# K-means cluster profiles across conditions, one panel per cluster.
# Line = cluster mean of per-gene z-scores, band = interquartile range,
# split by genotype so WT/mutant separation (or lack of it) is visible.

suppressPackageStartupMessages({
  library(optparse); library(ggplot2); library(dplyr); library(tidyr)
})

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "palette.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--vst",      type = "character"),
  make_option("--meta",     type = "character"),
  make_option("--clusters", type = "character"),
  make_option("--labels",   type = "character", default = NULL,
              help = "Optional TSV (cluster, label) annotating each panel"),
  make_option("--indiv-dir", type = "character", default = NULL, dest = "indiv_dir",
              help = "If set, also write one figure per cluster here (cluster_N_profile.{pdf,png})"),
  make_option("--short-stages", action = "store_true", default = FALSE, dest = "short_stages",
              help = paste("Abbreviate the x tick labels to L1/L4/D3/D6. At 10pt",
                           "the full 'Day 0' labels are wider than a facet and",
                           "collide -- this is horizontal crowding at a fixed",
                           "placed width, so a bigger canvas cannot fix it.")),
  make_option("--ncol", type = "integer", default = 3,
              help = "Facet columns; 2 for k = 4 so no facet hangs alone on a second row."),
  make_option("--width", type = "double", default = 9.5, help = "Canvas width in inches."),
  make_option("--font-scale", type = "double", default = 1, dest = "font_scale",
              help = "Multiply every font size; render at the ORIGINAL canvas."),
  make_option("--outdir",   type = "character")
)))
stopifnot(!is.null(opt$vst), !is.null(opt$meta), !is.null(opt$clusters), !is.null(opt$outdir))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

vst  <- read.csv(opt$vst, row.names = 1, check.names = FALSE)
meta <- read.delim(opt$meta, stringsAsFactors = FALSE)
km   <- read.csv(opt$clusters, stringsAsFactors = FALSE)

# Accepts both new (WT/tm1978) and legacy (N2/hlh-30) tokens.
geno_label  <- c(N2 = "WT", WT = "WT",
                 `hlh-30` = "hlh-30(tm1978)", tm1978 = "hlh-30(tm1978)")
stage_label <- if (opt$short_stages) {
  c(L1 = "L1", L4 = "L4", Day3 = "D3", Day6 = "D6")
} else {
  c(L1 = "L1", L4 = "Day 0", Day3 = "Day 3", Day6 = "Day 6")
}
relabel <- function(x, map) ifelse(x %in% names(map), map[x], x)

meta$cond <- paste(meta$genotype, meta$stage, sep = "_")
conds <- unique(meta$cond)

m    <- as.matrix(vst)
prof <- sapply(conds, function(c) rowMeans(m[, meta$sample_id[meta$cond == c], drop = FALSE]))
z    <- t(scale(t(prof)))
z    <- z[complete.cases(z), , drop = FALSE]

km <- km[km$gene_id %in% rownames(z), ]
z  <- z[km$gene_id, , drop = FALSE]

cond_meta <- meta[!duplicated(meta$cond), c("cond", "genotype", "stage")]

long <- as.data.frame(z) |>
  mutate(gene = rownames(z), cluster = km$cluster) |>
  pivot_longer(all_of(conds), names_to = "cond", values_to = "zscore") |>
  left_join(cond_meta, by = "cond") |>
  mutate(
    geno  = factor(relabel(genotype, geno_label),
                   levels = c("WT", "hlh-30(tm1978)")),
    stg   = relabel(stage, stage_label),
    # levels come from stage_label so they follow --short-stages; hardcoding
    # the long forms here made every row NA when the labels were abbreviated
    stg   = factor(stg, levels = intersect(unname(stage_label), unique(stg)))
  )

summ <- long |>
  group_by(cluster, geno, stg) |>
  summarise(mid = mean(zscore), lo = quantile(zscore, .25), hi = quantile(zscore, .75),
            .groups = "drop")

sizes <- km |> count(cluster, name = "n")
ann <- if (!is.null(opt$labels) && file.exists(opt$labels)) {
  read.delim(opt$labels, stringsAsFactors = FALSE)
} else {
  data.frame(cluster = sizes$cluster, label = "")
}
sizes <- left_join(sizes, ann, by = "cluster")
sizes$label[is.na(sizes$label)] <- ""
# "#k Label" wrapped to the facet width, then "(n = ...)" on its own line.
# The label is wrapped here, not in the labeller, which would re-flow the
# line break too.
lab_wrap <- if (opt$ncol <= 2) 28 else 16
sizes$panel <- ifelse(
  nzchar(sizes$label),
  sprintf("%s\n(n = %s)", stringr::str_wrap(sprintf("#%s %s", sizes$cluster, sizes$label), lab_wrap),
          format(sizes$n, big.mark = ",", trim = TRUE)),
  sprintf("#%s\n(n = %s)", sizes$cluster, format(sizes$n, big.mark = ",", trim = TRUE))
)

# Panels in numerical cluster order.
sizes <- sizes[order(sizes$cluster), ]
summ$panel <- factor(sizes$panel[match(summ$cluster, sizes$cluster)], levels = sizes$panel)

# hlh-30 in a lighter purple than PAL$hlh30, so it stays distinct from the
# dark WT line where the two genotypes track each other
pal <- c("WT" = PAL$wt, "hlh-30(tm1978)" = "#A77BBF")
geno_expr <- parse(text = c('"WT"', 'italic("hlh-30(tm1978)")'))

# shared layers/scales so the combined and per-cluster figures match exactly
build <- function(df, faceted, title = NULL) {
  g <- ggplot(df, aes(stg, mid, group = geno)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = geno), alpha = 0.13, colour = NA) +
    geom_line(aes(colour = geno, linetype = geno), linewidth = 0.9) +
    geom_point(aes(colour = geno, shape = geno), size = 2.4) +
    # explicit breaks: labels are positional
    scale_colour_manual(values = pal, name = NULL, breaks = names(pal), labels = geno_expr) +
    scale_fill_manual(values = pal,  name = NULL, breaks = names(pal), labels = geno_expr) +
    scale_linetype_manual(values = c("WT" = "solid", "hlh-30(tm1978)" = "22"),
                          name = NULL, breaks = names(pal), labels = geno_expr) +
    scale_shape_manual(values = c("WT" = 16, "hlh-30(tm1978)" = 17),
                       name = NULL, breaks = names(pal), labels = geno_expr) +
    # whole-number ticks only: at 10pt the half steps collide in a short facet
    scale_y_continuous(breaks = scales::breaks_width(1)) +
    labs(x = NULL, y = "z-scored expression", title = title) +
    theme_bw(base_size = 11 * opt$font_scale) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.border = element_rect(colour = "grey80"),
          strip.background = element_rect(fill = "grey96", colour = "grey80"),
          strip.text = element_text(size = 8.5 * opt$font_scale, lineheight = 1.0, colour = "grey20",
                                    margin = margin(4.5, 3, 4.5, 3)),   # air above and below the title
          plot.title = element_text(face = "bold", size = 11 * opt$font_scale),
          legend.position = "top", legend.justification = "left",
          legend.key.width = unit(22, "pt"),
          # no dead band between the legend and the first row of facets
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
          plot.margin = margin(4, 6, 2, 4),
          axis.title.y = element_text(size = 9 * opt$font_scale, colour = "grey30"))
  if (faceted) g + facet_wrap(~ panel, ncol = opt$ncol,
                              labeller = label_value)
  else g
}

# combined (all clusters)
p <- build(summ, faceted = TRUE)
# three-line strips need the extra height on the 3-column layout
cp_h <- if (opt$ncol >= 3) 7.0 else 6.2
ggsave(file.path(opt$outdir, "cluster_profiles.pdf"), p, width = opt$width, height = cp_h)
ggsave(file.path(opt$outdir, "cluster_profiles.png"), p, width = opt$width, height = cp_h, dpi = 200)
message("wrote: ", file.path(opt$outdir, "cluster_profiles.{pdf,png}"))

# per-cluster individual figures (same design), if requested
if (!is.null(opt$indiv_dir)) {
  dir.create(opt$indiv_dir, showWarnings = FALSE, recursive = TRUE)
  for (cl in sort(unique(summ$cluster))) {
    ttl <- as.character(sizes$panel[sizes$cluster == cl])
    pi <- build(summ[summ$cluster == cl, ], faceted = FALSE, title = gsub("\n", " - ", ttl))
    base <- file.path(opt$indiv_dir, sprintf("cluster_%s_profile", cl))
    ggsave(paste0(base, ".pdf"), pi, width = 5.0, height = 4.2)
    ggsave(paste0(base, ".png"), pi, width = 5.0, height = 4.2, dpi = 200)
  }
  message("wrote ", length(unique(summ$cluster)), " per-cluster figures to ", opt$indiv_dir)
}
