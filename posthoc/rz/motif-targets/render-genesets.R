#!/usr/bin/env Rscript
# Dumbbells of "% direct HLH-30 targets" (E-box in promoter) for the compare-DE
# category lists and the k-means GO-enriched gene sets, from the same cached
# scan / summary table as the DEG dumbbell. Same background line + Fisher stars.

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(ggplot2); library(stringr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--summary", type = "character", default = "out/hlh30-all-ruv/motif_targets_summary.tsv"),
  make_option("--labels",  type = "character", default = "../cluster-profiles/labels-hlh30-all-ruv.tsv"),
  make_option("--outdir",  type = "character", default = "out/hlh30-all-ruv"),
  make_option("--accent",  type = "character", default = "#7B3294")
)))
ACC <- opt$accent; GREY <- "#9A9AA0"
res <- read_tsv(opt$summary, show_col_types = FALSE)
bg  <- res$bg_pct[1]
res$star[is.na(res$star)] <- ""

cap <- str_wrap(sprintf(
  "Background = %.1f%% of expressed genes carry ~CAYRTG (JASPAR MA1449.1, p<1e-4, promoter -2000/+200). Grey dot = background; purple = the gene set. Stars = BH-adjusted Fisher vs background.",
  bg), width = 100)
base_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "grey93"),
        panel.spacing.y = unit(8, "pt"),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8.5, colour = "grey25", lineheight = 1.1),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.5, colour = "grey40"),
        plot.caption = element_text(size = 7.5, colour = "grey55", hjust = 0))

# grouped dumbbell: facet by `group`, sub-rows = `sub`
dumbbell <- function(df, title, subtitle, w, h, file) {
  df <- df %>% mutate(bar_end = pmax(pct, bg))
  p <- ggplot(df) +
    geom_vline(xintercept = bg, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_segment(aes(x = bg_pct, xend = pct, y = sub, yend = sub), colour = "grey80", linewidth = 0.9) +
    geom_point(aes(bg_pct, sub), colour = GREY, size = 2.1) +
    geom_point(aes(pct, sub), colour = ACC, size = 3.2) +
    geom_text(aes(bar_end, sub, label = sprintf("  %.0f%%%s", pct, star)),
              hjust = 0, size = 2.8, colour = "grey20") +
    facet_grid(grp_lab ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.16))) +
    labs(x = "genes with ~CAYRTG in promoter (%)", y = NULL,
         title = title, subtitle = subtitle, caption = cap) +
    base_theme
  ggsave(file.path(opt$outdir, paste0(file, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(opt$outdir, paste0(file, ".png")), p, width = w, height = h, dpi = 200)
  message("wrote ", file, ".{pdf,png}")
}

# =================================================== compare-DE dumbbell ======
cat_lab <- c(tm1978_L4_vs_tm1978_L1_exclusive = "hlh-30(tm1978)-only\n(L1 -> L4)",
             WT_L4_vs_WT_L1_exclusive = "WT-only\n(L1 -> L4)",
             shared = "shared\n(both L1 -> L4)",
             differential = "differential\n(L1 -> L4)")
cmp <- res %>% filter(family == "compare") %>%
  mutate(grp_lab = factor(cat_lab[group], levels = unname(cat_lab)),
         sub = factor(c(up = "up", down = "down")[direction], levels = c("down","up")))
dumbbell(cmp,
  "Direct HLH-30 targets among compare-DE gene lists",
  "Genes split by category (WT-only / hlh-30-only / shared / differential across the L1->L4 change).",
  7.4, 5.2, "targets_compareDE_dumbbell")

# ================================================= GO (per cluster) dumbbell ==
lab <- read_tsv(opt$labels, show_col_types = FALSE)
clab <- setNames(lab$label, paste0("cluster_", lab$cluster, "_genes"))
go <- res %>% filter(family == "go") %>%
  mutate(cl_key = sub(" (BP|CC|MF)$", "", group),
         ont    = sub(".* ", "", group),
         cl_num = as.integer(sub("cluster_(\\d+)_genes", "\\1", cl_key)),
         grp_lab = factor(sprintf("Cluster %d", cl_num),
                          levels = sprintf("Cluster %d", sort(unique(cl_num)))),
         sub = factor(ont, levels = c("MF","CC","BP")))
dumbbell(go,
  "Direct HLH-30 targets among GO-enriched cluster genes",
  "Per k-means cluster: union of genes across ALL its significantly enriched GO terms (BP / CC / MF), not any single function.",
  7.4, 6.6, "targets_GO_dumbbell")
