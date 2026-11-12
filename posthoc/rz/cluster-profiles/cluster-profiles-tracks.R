#!/usr/bin/env Rscript
# K-means cluster profiles, two-genotype style (WT vs hlh-30(tm1978)) but with
# the individual per-gene tracks kept behind the genotype mean lines - a merge
# of the pipeline figure (grey gene tracks) and the posthoc figure (two lines).
# Each gene contributes two faint lines (one per genotype); bold lines are the
# per-genotype cluster mean of the z-scores.

suppressPackageStartupMessages({
  library(optparse); library(ggplot2); library(dplyr); library(tidyr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--vst",      type = "character", default = "cache/vst_counts.csv"),
  make_option("--meta",     type = "character", default = "cache/sample_meta.tsv"),
  make_option("--clusters", type = "character", default = "cache/cluster_assignments.csv"),
  make_option("--labels",   type = "character", default = "labels-hlh30-all-ruv.tsv"),
  make_option("--track-alpha", type = "double", default = 0.015, dest = "track_alpha"),
  make_option("--outdir",   type = "character", default = "out"),
  make_option("--indiv-dir", type = "character", default = "out/individual", dest = "indiv_dir",
              help = "per-cluster figures go here (destined for the results/ folder)")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(opt$indiv_dir, showWarnings = FALSE, recursive = TRUE)

vst  <- read.csv(opt$vst, row.names = 1, check.names = FALSE)
meta <- read.delim(opt$meta, stringsAsFactors = FALSE)
km   <- read.csv(opt$clusters, stringsAsFactors = FALSE)

geno_label  <- c(N2 = "WT", WT = "WT", `hlh-30` = "hlh-30(tm1978)", tm1978 = "hlh-30(tm1978)")
stage_label <- c(L1 = "L1", L4 = "Day 0", Day3 = "Day 3", Day6 = "Day 6")
relabel <- function(x, map) ifelse(x %in% names(map), map[x], x)

meta$cond <- paste(meta$genotype, meta$stage, sep = "_")
conds <- unique(meta$cond)

# per-gene mean per condition, z-scored across the 8 conditions
m    <- as.matrix(vst)
prof <- sapply(conds, function(c) rowMeans(m[, meta$sample_id[meta$cond == c], drop = FALSE]))
z    <- t(scale(t(prof)))
z    <- z[complete.cases(z), , drop = FALSE]
km   <- km[km$gene_id %in% rownames(z), ]
z    <- z[km$gene_id, , drop = FALSE]

cond_meta <- meta[!duplicated(meta$cond), c("cond", "genotype", "stage")]
long <- as.data.frame(z) |>
  mutate(gene = rownames(z), cluster = km$cluster) |>
  pivot_longer(all_of(conds), names_to = "cond", values_to = "zscore") |>
  left_join(cond_meta, by = "cond") |>
  mutate(geno = factor(relabel(genotype, geno_label), levels = c("WT", "hlh-30(tm1978)")),
         stg  = relabel(stage, stage_label),
         stg  = factor(stg, levels = intersect(c("L1","Day 0","Day 3","Day 6"), unique(stg))))

summ <- long |> group_by(cluster, geno, stg) |>
  summarise(mid = mean(zscore), .groups = "drop")

# panel labels: "Cluster N (n = ...)" + optional annotation
sizes <- km |> count(cluster, name = "n")
ann <- if (!is.null(opt$labels) && file.exists(opt$labels))
  read.delim(opt$labels, stringsAsFactors = FALSE) else data.frame(cluster = sizes$cluster, label = "")
sizes <- left_join(sizes, ann, by = "cluster"); sizes$label[is.na(sizes$label)] <- ""
sizes$panel <- ifelse(nzchar(sizes$label),
  sprintf("Cluster %s  (n = %s)\n%s", sizes$cluster, format(sizes$n, big.mark = ","), sizes$label),
  sprintf("Cluster %s  (n = %s)", sizes$cluster, format(sizes$n, big.mark = ",")))
sizes <- sizes[order(sizes$cluster), ]
pf <- function(cl) factor(sizes$panel[match(cl, sizes$cluster)], levels = sizes$panel)
long$panel <- pf(long$cluster); summ$panel <- pf(summ$cluster)

pal <- c("WT" = "#1F78B4", "hlh-30(tm1978)" = "#E66101")
geno_expr <- parse(text = c('"WT"', 'italic("hlh-30(tm1978)")'))

# shared layers/scales so the combined and per-cluster figures match exactly
build <- function(long_df, summ_df, faceted, title = NULL) {
  g <- ggplot() +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_line(data = long_df, aes(stg, zscore, colour = geno, group = interaction(gene, geno)),
              linewidth = 0.12, alpha = opt$track_alpha) +
    geom_line(data = summ_df, aes(stg, mid, colour = geno, linetype = geno, group = geno),
              linewidth = 1.1) +
    geom_point(data = summ_df, aes(stg, mid, colour = geno, shape = geno, group = geno), size = 2.3) +
    scale_colour_manual(values = pal, name = NULL, labels = geno_expr) +
    scale_linetype_manual(values = c("WT" = "solid", "hlh-30(tm1978)" = "22"),
                          name = NULL, labels = geno_expr) +
    scale_shape_manual(values = c("WT" = 16, "hlh-30(tm1978)" = 17), name = NULL, labels = geno_expr) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1.1))) +
    labs(x = NULL, y = "z-scored expression", title = title) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.border = element_rect(colour = "grey80"),
          strip.background = element_rect(fill = "grey96", colour = "grey80"),
          strip.text = element_text(size = 8.5, lineheight = 1.15, colour = "grey20"),
          plot.title = element_text(face = "bold", size = 11.5),
          legend.position = "inside",
          legend.position.inside = c(0.012, 0.986),
          legend.justification.inside = c(0, 1),
          legend.direction = "horizontal",
          legend.background = element_rect(fill = scales::alpha("white", 0.72), colour = NA),
          legend.key = element_blank(),
          legend.key.width = unit(20, "pt"),
          legend.key.height = unit(12, "pt"),
          legend.text = element_text(size = 8.3),
          legend.margin = margin(2, 5, 3, 3),
          axis.title.y = element_text(size = 9, colour = "grey30"))
  if (faceted) g + facet_wrap(~ panel, ncol = 3) else g
}

# --- combined 6-panel (posthoc) ---
pc <- build(long, summ, faceted = TRUE,
            title = "K-means cluster profiles: WT vs hlh-30(tm1978) with per-gene tracks")
ggsave(file.path(opt$outdir, "cluster_profiles_tracks.pdf"), pc, width = 9.5, height = 6.4, device = cairo_pdf)
ggsave(file.path(opt$outdir, "cluster_profiles_tracks.png"), pc, width = 9.5, height = 6.4, dpi = 200)
message("wrote cluster_profiles_tracks.{pdf,png} to ", opt$outdir)

# --- per-cluster individual figures (destined for the results/ folder) ---
for (cl in sort(unique(km$cluster))) {
  ld <- long[long$cluster == cl, ]; sd <- summ[summ$cluster == cl, ]
  ttl <- as.character(sizes$panel[sizes$cluster == cl])
  pi <- build(ld, sd, faceted = FALSE, title = gsub("\n", " - ", ttl))
  base <- file.path(opt$indiv_dir, sprintf("cluster_%s_profile_bygenotype", cl))
  ggsave(paste0(base, ".pdf"), pi, width = 5.2, height = 4.4, device = cairo_pdf)
  ggsave(paste0(base, ".png"), pi, width = 5.2, height = 4.4, dpi = 200)
}
message("wrote ", length(unique(km$cluster)), " per-cluster figures to ", opt$indiv_dir)
