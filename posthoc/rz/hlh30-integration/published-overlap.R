#!/usr/bin/env Rscript
# Our hlh-30(tm1978) vs WT DEGs against published HLH-30 gene lists
# (published/, published/new_degsets/, published/new_autophagy/; see the
# README.md / NOTES.md in each for sources and caveats).
#
#   published_overlap   every published list (rows, in four groups) among our
#                       up / down DEGs per stage: % of the list that is one of
#                       our DEGs, Fisher's exact test (universe = genes with a
#                       padj at that stage; lists reduced to it), BH over all
#                       cells.
#   published_autophagy log2FC at our four stages for (B1) core autophagy genes
#                       that papers show to need hlh-30 for expression and
#                       (B2) every gene of the Zhang & Melendez 2025 core
#                       autophagy table, regardless of hlh-30.
#
#   Rscript published-overlap.R          # -> out/published/
#   Rscript published-overlap.R --paper  # -> out/published/paper/03-*, 04-*
#                                          (7.09 in, 10 pt, no titles or notes)

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "palette.R"))
opt <- parse_args(OptionParser(option_list = list(
  make_option("--deg-dir", type = "character", default = "../volcano/cache/hlh30-all-ruv", dest = "deg_dir"),
  make_option("--pub",     type = "character", default = "published"),
  make_option("--padj",    type = "double", default = 0.05),
  make_option("--lfc",     type = "double", default = 1),
  make_option("--outdir",  type = "character", default = "out/published"),
  make_option("--paper",   action = "store_true", default = FALSE,
              help = "paper versions only: full page width, 10 pt, no titles, subtitles or captions")
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
STAGES <- c(L1 = "L1", L4 = "L4", D3 = "Day3", D6 = "Day6")

# ---- list registry -----------------------------------------------------------
# group, label, file (under --pub), direction values kept. "lower" / "higher"
# = in the hlh-30 sample; group 3 lists are induced in WT and need hlh-30.
G <- c("hlh-30 loss vs WT, whole animal",
       "hlh-30 loss in a modified background or one tissue",
       "Induced by a stimulus in WT, needs hlh-30",
       "Binding, motif or protein evidence")
R <- function(g, lab, file, dir = NULL) data.frame(group = G[g], label = lab, file = file,
                                                    dir = if (is.null(dir)) NA else paste(dir, collapse = ","))
N <- "new_degsets/"
reg <- rbind(
  R(1, "Grove 2009: embryo-L2, fed; lower in hlh-30",        paste0(N, "grove2009_hlh30_microarray_de.tsv"), "down"),
  R(1, "Grove 2009: embryo-L2, fed; higher in hlh-30",       paste0(N, "grove2009_hlh30_microarray_de.tsv"), "up"),
  R(1, "Gerisch 2020: fed adults; lower in hlh-30",          paste0(N, "gerisch2020_hlh30_vs_wt_al_de.tsv"), "down"),
  R(1, "Gerisch 2020: fed adults; higher in hlh-30",         paste0(N, "gerisch2020_hlh30_vs_wt_al_de.tsv"), "up"),
  R(1, "Gerisch 2020: starved adults (diapause); lower",     paste0(N, "gerisch2020_hlh30_vs_wt_ard_de.tsv"), "down"),
  R(1, "Gerisch 2020: starved adults (diapause); higher",    paste0(N, "gerisch2020_hlh30_vs_wt_ard_de.tsv"), "up"),
  R(1, "Munoz-Barrera 2026: starved L1; lower in hlh-30",    paste0(N, "munozbarrera2026_hlh30_vs_wt_starved_l1_de.tsv"), "down"),
  R(1, "Munoz-Barrera 2026: starved L1; higher in hlh-30",   paste0(N, "munozbarrera2026_hlh30_vs_wt_starved_l1_de.tsv"), "up"),
  R(1, "Murphy 2019: starved L1 (KEGG genes); lower",        paste0(N, "murphy2019_hlh30_starved_l1_kegg_clusters.tsv"), "down"),
  R(1, "Murphy 2019: starved L1 (KEGG genes); higher",       paste0(N, "murphy2019_hlh30_starved_l1_kegg_clusters.tsv"), "up_in_hlh30"),
  R(1, "Sohn 2025: fed adult males (male genes); lower",     paste0(N, "sohn2025_hlh30_down_in_males.tsv"), "down"),
  R(2, "Shioda 2023: neuronal hlh-30 RNAi, glp-1; lower",    paste0(N, "shioda2023_glp1_neuronal_hlh30_rnai_de.tsv"), "down"),
  R(2, "Shioda 2023: neuronal hlh-30 RNAi, glp-1; higher",   paste0(N, "shioda2023_glp1_neuronal_hlh30_rnai_de.tsv"), "up"),
  R(2, "Lin 2024: Abeta strain GMC101; lower",               paste0(N, "lin2024_gmc101_hlh30_ko_de.tsv"), "down"),
  R(2, "Lin 2024: Abeta strain GMC101; higher",              paste0(N, "lin2024_gmc101_hlh30_ko_de.tsv"), "up"),
  R(2, "Metcalf 2024: glial XBP-1s; lower",                  paste0(N, "metcalf2024_glial_xbp1s_hlh30_vs_xbp1s_de.tsv"), "down"),
  R(2, "Metcalf 2024: glial XBP-1s; higher",                 paste0(N, "metcalf2024_glial_xbp1s_hlh30_vs_xbp1s_de.tsv"), "up"),
  R(3, "Visvikis 2014: S. aureus",                           "visvikis2014_hlh30_dependent_saureus_induced.tsv"),
  R(3, "Wani 2021: S. aureus vs starved",                    paste0(N, "wani2021_hlh30_dependent_saureus_induced.tsv")),
  R(3, "Chen 2017: Cry5B toxin",                             paste0(N, "chen2017_hlh30_dependent_cry5b_induced.tsv")),
  R(3, "Hsieh 2026: Cry5B toxin (HLH-30 R38K)",              paste0(N, "hsieh2026_hlh30_r38_dependent_cry5b_induced.tsv")),
  R(3, "Possik 2023: pgph-2 overexpression",                 paste0(N, "possik2023_hlh30_dependent_pgph2_induced.tsv")),
  R(3, "Murphy 2019: starved L1 (KEGG genes)",               paste0(N, "murphy2019_hlh30_starved_l1_kegg_clusters.tsv"), "up"),
  R(3, "Visvikis 2014: S. aureus, NOT hlh-30-dependent",     "visvikis2014_hlh30_independent_saureus_induced.tsv"),
  R(3, "Wani 2021: S. aureus, NOT hlh-30-dependent",         paste0(N, "wani2021_hlh30_independent_saureus_induced.tsv")),
  R(4, "ENCODE: HLH-30::GFP ChIP, L4, promoter",             "encode_modern_hlh30_l4_chip_promoter_bound.tsv"),
  R(4, "ENCODE: HLH-30::GFP ChIP, late embryo, promoter",    paste0(N, "encode_modern_hlh30_late_embryo_chip_promoter_bound.tsv")),
  R(4, "Grove 2009: predicted from HLH-30 binding (PBM)",    "grove2009_hlh30_pbm_predicted_targets.tsv"),
  R(4, "Harvald 2017: starvation proteins, CLEAR motif",     "harvald2017_starvation_hlh30_protein_highB_clear.tsv"))
reg$key <- sprintf("r%02d", seq_len(nrow(reg)))
# direction of the published change moves out of the label into its own
# glyph: "; lower ..." / "; higher ..." = lower / higher in the hlh-30 sample
reg$pdir <- ifelse(grepl("; lower", reg$label), "down", ifelse(grepl("; higher", reg$label), "up", NA))
reg$label <- sub(";\\s*(lower|higher).*$", "", reg$label)
read_list <- function(i) {
  d <- read.delim(file.path(opt$pub, reg$file[i]), stringsAsFactors = FALSE, quote = "")
  if (!is.na(reg$dir[i])) d <- d[d$direction %in% strsplit(reg$dir[i], ",")[[1]], ]
  unique(d$gene_id[!is.na(d$gene_id) & nzchar(d$gene_id)])
}
pub <- setNames(lapply(seq_len(nrow(reg)), read_list), reg$key)

dge <- bind_rows(lapply(names(STAGES), function(st) {
  read.csv(file.path(opt$deg_dir, sprintf("tm1978_%s_vs_WT_%s.csv", STAGES[[st]], STAGES[[st]])),
           stringsAsFactors = FALSE) %>%
    filter(!is.na(padj)) %>%
    transmute(gene_id, gene_name, stage = st, lfc = log2FoldChange, padj,
              deg = case_when(padj < opt$padj & log2FoldChange >=  opt$lfc ~ "up",
                              padj < opt$padj & log2FoldChange <= -opt$lfc ~ "down", TRUE ~ "ns"))
}))

res <- bind_rows(lapply(names(STAGES), function(st) {
  d <- dge[dge$stage == st, ]
  bind_rows(lapply(names(pub), function(k) {
    inl <- d$gene_id %in% pub[[k]]
    bind_rows(lapply(c("up", "down"), function(dir) {
      isd <- d$deg == dir
      tb <- table(factor(inl, c(TRUE, FALSE)), factor(isd, c(TRUE, FALSE)))
      ft <- fisher.test(tb)
      data.frame(stage = st, key = k, direction = dir, n_list = sum(inl), n_deg = sum(isd),
                 overlap = tb[1, 1], expected = sum(inl) * sum(isd) / nrow(d),
                 or = unname(ft$estimate), p = ft$p.value)
    }))
  }))
})) %>% mutate(q = p.adjust(p, "BH")) %>% left_join(reg %>% select(key, group, label, file, dir), by = "key")
write.table(res, file.path(opt$outdir, "published_overlap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- A: heatmap ---------------------------------------------------------------
# Columns: our DEGs, one strip per stage with down / up under it, and the
# chance rate (% of all expressed genes that are DEGs in that set, the
# reference each % in the column is read against) under each. Rows: the
# published lists, with list size and a direction glyph next to the heatmap.
nl <- res %>% group_by(key) %>% summarise(n = max(n_list))
GLYPH <- c(down = sprintf("<span style='color:%s'>&#9660;</span>", PAL$down),
           up   = sprintf("<span style='color:%s'>&#9650;</span>", PAL$up))
# rows without a direction get an invisible glyph, so the sizes stay aligned
reg$row <- sprintf("%s <span style='color:grey45'>%s</span> %s", reg$label,
                   format(nl$n[match(reg$key, nl$key)], big.mark = ",", trim = TRUE),
                   ifelse(is.na(reg$pdir), "<span style='color:white'>&#9660;</span>", GLYPH[reg$pdir]))
base <- dge %>% group_by(stage) %>%
  summarise(down = 100 * mean(deg == "down"), up = 100 * mean(deg == "up"), .groups = "drop") %>%
  pivot_longer(c(down, up), names_to = "direction", values_to = "bg")
# header: our direction as the same glyph, the chance rate under it
COLLAB <- with(base, setNames(sprintf("%s<br><span style='color:grey40'>%s%%</span>", GLYPH[direction],
                                      ifelse(bg < 1, sprintf("%.1f", bg), sprintf("%.0f", bg))),
                              paste(stage, direction)))
hm <- res %>%
  mutate(col = factor(paste(stage, direction), levels = paste(rep(names(STAGES), each = 2), c("down", "up"))),
         stage = factor(stage, levels = names(STAGES)),
         row = factor(reg$row[match(key, reg$key)], levels = rev(reg$row)),
         group = factor(group, levels = G),
         # colour only significant cells: a 0 overlap with < 1 gene expected is
         # not evidence of depletion
         l2or = ifelse(q < 0.05, pmax(pmin(log2(pmax(or, 1e-3)), 4), -4), NA),
         star = case_when(q < 0.001 ~ "***", q < 0.01 ~ "**", q < 0.05 ~ "*", TRUE ~ ""),
         pct = 100 * overlap / pmax(n_list, 1),
         # only significant cells carry a number, so the eye goes to them; the
         # full table is published_overlap.tsv
         lab = ifelse(q < 0.05, case_when(pct < 0.1 ~ "<0.1%", pct < 1 ~ sprintf("%.1f%%", pct),
                                          TRUE ~ sprintf("%.0f%%", pct)), ""),
         txt = ifelse(!is.na(l2or) & abs(l2or) > 2.5, "white", "grey10"))
# one block per group, the group name as a left-aligned heading over its
# rows; the stage and direction header is its own block on top
pg <- function(g, first = FALSE) {
  ggplot(hm %>% filter(group == g) %>% droplevels(), aes(col, row, fill = l2or)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = lab, colour = I(txt)), size = 10 / .pt) +
    facet_grid(. ~ stage, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = PAL$div[["low"]], mid = "white", high = PAL$div[["high"]], midpoint = 0,
                         limits = c(-4, 4), breaks = c(-4, -2, 0, 2, 4),
                         labels = c("1/16", "1/4", "1", "4", "16"), name = "odds ratio\nvs chance",
                         na.value = "grey95") +
    scale_x_discrete(position = "top", labels = function(v) COLLAB[v]) +
    labs(x = NULL, y = NULL, title = g) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "right",
          plot.title.position = "plot",
          plot.title = element_text(size = 10, face = "bold", colour = "grey15", hjust = 0,
                                    margin = margin(t = 6, b = 2)),
          axis.text.y.left = ggtext::element_markdown(size = 10, colour = "grey10"),
          axis.text.x.top = element_blank(), strip.text.x = element_blank(),
          panel.spacing.x = unit(6, "pt"),
          legend.title = element_text(size = 10), legend.text = element_text(size = 10),
          plot.margin = margin(0, 2, 0, 2))
}
# header block: stage strips above the direction glyphs and chance rates
hdr <- ggplot(distinct(hm, col, stage), aes(col, "")) + geom_blank() +
  facet_grid(. ~ stage, scales = "free_x", space = "free_x") +
  scale_x_discrete(position = "top", labels = function(v) COLLAB[v]) +
  labs(x = NULL, y = NULL) + theme_void(base_size = 10) +
  theme(strip.placement = "outside", strip.text.x = element_text(size = 10, face = "bold", margin = margin(b = 3)),
        axis.text.x.top = ggtext::element_markdown(size = 10, lineheight = 1.1),
        panel.spacing.x = unit(6, "pt"), plot.margin = margin(0, 2, 0, 2))
blocks <- c(list(hdr), lapply(G, pg))
nrows <- sapply(G, function(g) sum(reg$group == g))
# heights in rows; each group block also holds its heading
pA <- wrap_plots(blocks, ncol = 1, heights = c(0.1, nrows + 1)) + plot_layout(guides = "collect")
capA <- paste0("Our DEGs: hlh-30(tm1978) vs WT at each stage (padj < 0.05, |log2FC| >= 1). Cell = % of the published list ",
               "that is in our DEG set (Grove 2009 up, D3: 2 of its 11 genes = 18%), shown only when the overlap differs from chance (Fisher's exact test, BH over all ",
               "cells, q < 0.05); the % under each column glyph is the chance rate (our DEGs as a % of all expressed genes). Universe: genes with a DESeq2 padj at that ",
               "stage; lists reduced to it. Grey number = list size. Glyph next to a list = lower (down) or higher (up) in ",
               "the hlh-30 sample of that study; no glyph = the list has no direction. The two 'NOT hlh-30-dependent' rows are the controls ",
               "for the stimulus lists above them. Some cut-offs were set by us where the paper gives only a full table ",
               "(published/new_degsets/NOTES.md). Nonninger 2025 is left out (fold changes, no p values). ENCODE promoter ",
               "assignment (summit -1000/+200 of a TSS) is ours.")
figA <- pA + plot_annotation(title = "Our hlh-30 DEGs vs published HLH-30 gene lists (draft for discussion)",
                             subtitle = "Cell: % of the published list that is one of our DEGs (not a gene count). Under each column: the same % expected by chance.",
                             caption = stringr::str_wrap(capA, 150),
                             theme = theme(plot.title = element_text(face = "bold", size = 12),
                                           plot.subtitle = element_text(size = 10, colour = "grey20", margin = margin(b = 6)),
                                           plot.caption = element_text(colour = "grey35", hjust = 0, size = 10,
                                                                       lineheight = 1.1)))
hA <- 3.1 + 0.24 * nrow(reg)
if (!opt$paper) {
  ggsave(file.path(opt$outdir, "published_overlap.pdf"), figA, width = 10, height = hA, device = cairo_pdf)
  ggsave(file.path(opt$outdir, "published_overlap.png"), figA, width = 10, height = hA, dpi = 150)
}

# ---- B: autophagy genes gene by gene --------------------------------------------
AU <- file.path(opt$pub, "new_autophagy")
ROLE <- c(core_machinery = "core machinery", maturation_fusion = "maturation / fusion",
          selective_autophagy = "selective autophagy", transcriptional_regulator = "regulate autophagy genes",
          autophagy_related_or_regulator = "other autophagy-related")
lnk <- read.delim(file.path(AU, "hlh30_autophagy_genes_summary.tsv"), stringsAsFactors = FALSE, quote = "")
core <- read.delim(file.path(AU, "core_autophagy_zhang_melendez2025.tsv"), stringsAsFactors = FALSE, quote = "") %>%
  filter(nzchar(gene_id)) %>% distinct(gene_id, .keep_all = TRUE)
# B1: core autophagy genes (in the review table) whose expression papers show to need hlh-30
b1 <- lnk %>% filter(n_papers_expr_reduced_hlh30 > 0, in_core_list %in% c("True", "TRUE", TRUE)) %>%
  left_join(core %>% select(gene_id, role_class), by = "gene_id") %>%
  transmute(gene_id, role_class, n_pap = n_papers_expr_reduced_hlh30,
            lab = sprintf("%s (%d)", gene_name, n_papers_expr_reduced_hlh30))
b2 <- core %>% transmute(gene_id, role_class,
                         lab = ifelse(gene_name == "hlh-30", "hlh-30 (deleted)", gene_name))
tile_df <- function(genes) {
  genes %>% tidyr::crossing(stage = names(STAGES)) %>%
    left_join(dge %>% select(gene_id, stage, lfc, padj, deg), by = c("gene_id", "stage")) %>%
    mutate(stage = factor(stage, levels = names(STAGES)),
           role = factor(ROLE[role_class], levels = ROLE),
           sig = case_when(is.na(padj) ~ "", deg != "ns" ~ "**", padj < 0.05 ~ "*", TRUE ~ ""))
}
ord_lab <- function(d) {   # within a role, by mean log2FC (most down at the bottom)
  o <- d %>% group_by(role, lab) %>% summarise(m = mean(lfc, na.rm = TRUE), .groups = "drop") %>%
    arrange(role, desc(m))
  d %>% mutate(lab = factor(lab, levels = rev(unique(o$lab))))
}
tiles <- function(d, title, subtitle) {
  ggplot(d, aes(stage, lab, fill = pmax(pmin(lfc, 1.5), -1.5))) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sig), size = 3.2, vjust = 0.75) +
    facet_grid(role ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(14)) +
    scale_fill_gradient2(low = PAL$div[["low"]], mid = "white", high = PAL$div[["high"]], midpoint = 0,
                         limits = c(-1.5, 1.5), na.value = "grey85", name = expression(log[2] ~ "FC")) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text.y = element_text(face = "italic"),
          strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),
          strip.background = element_rect(fill = "grey96", colour = NA), panel.spacing.y = unit(4, "pt"),
          plot.title = element_text(face = "bold"))
}
d1 <- ord_lab(tile_df(b1)); d2 <- ord_lab(tile_df(b2))
p1 <- tiles(d1, "B1  Autophagy genes known to need hlh-30",
            sprintf("%d core autophagy genes shown to be lower\nin hlh-30 (mutant or RNAi) in at least one paper;\n(n) = number of papers", nrow(b1)))
p2 <- tiles(d2, "B2  All core autophagy genes",
            sprintf("%d genes of Table 1, Zhang & Melendez 2025\n(Genetics), regardless of hlh-30", nrow(b2)))
capB <- paste0("log2FC hlh-30(tm1978) vs WT at each stage. * padj < 0.05; ** DEG (padj < 0.05 and |log2FC| >= 1). ",
               "Grey = no padj at that stage (not expressed or filtered). Colour capped at +/-1.5. ",
               "B1 sources: 18 papers that perturb hlh-30 and test single genes (published/new_autophagy/hlh30_autophagy_evidence.tsv, ",
               "each row with its figure panel). Of B1's genes, ",
               sum(b1$gene_id %in% core$gene_id), " are also in B2.")
# same row height in both: B1 sits on top of an empty block
r1 <- nrow(b1) + 3 * n_distinct(b1$role_class); r2 <- nrow(b2) + 3 * n_distinct(b2$role_class)
figB <- ((p1 / plot_spacer() + plot_layout(heights = c(r1, r2 - r1))) | p2) +
  plot_layout(widths = c(1, 1), guides = "collect") +
  plot_annotation(title = "Autophagy genes in our hlh-30 data (draft for discussion)",
                  caption = stringr::str_wrap(capB, 150),
                  theme = theme(plot.title = element_text(face = "bold", size = 13),
                                plot.caption = element_text(colour = "grey40", hjust = 0, size = 8)))
hB <- 2.4 + 0.17 * nrow(b2)
if (!opt$paper) {
  ggsave(file.path(opt$outdir, "published_autophagy.pdf"), figB, width = 10, height = hB, device = cairo_pdf)
  ggsave(file.path(opt$outdir, "published_autophagy.png"), figB, width = 10, height = hB, dpi = 150)
}

# ---- paper versions ------------------------------------------------------------
# Same data, placed at full page width (7.09 in) with every label at 10 pt.
# The explanations move to the figure legend, so no titles or captions here.
if (opt$paper) {
  PW <- 7.09
  pdir <- file.path(opt$outdir, "paper")
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  ital <- function(x) gsub("hlh-30", "<i>hlh-30</i>", x, fixed = TRUE)
  paper_theme <- theme_minimal(base_size = 10) +
    theme(text = element_text(size = 10), panel.grid = element_blank(),
          legend.title = element_text(size = 10), legend.text = element_text(size = 10),
          legend.position = "bottom")

  # 03: one block per group under its heading, the stage header on top
  pgp <- function(g) {
    pg(g) + labs(title = ital(g)) + paper_theme +
      theme(plot.title.position = "plot",
            plot.title = ggtext::element_markdown(size = 10, face = "bold", colour = "grey15", hjust = 0,
                                                  margin = margin(t = 5, b = 2)),
            axis.text.y.left = ggtext::element_markdown(size = 10, colour = "grey10"),
            axis.text.x.top = element_blank(), strip.text.x = element_blank(),
            panel.spacing.x = unit(4, "pt"), plot.margin = margin(0, 2, 0, 2)) +
      guides(fill = guide_colourbar(title.position = "left", title.vjust = 0.9,
                                    barwidth = unit(1.6, "in"), barheight = unit(8, "pt")))
  }
  short <- function(x) {
    for (r in list(c(" (diapause)", ""), c(" (KEGG genes)", ""), c(" (male genes)", ""),
                   c(", promoter", ""), c("late embryo", "embryo"), c("HLH-30::GFP ChIP", "ChIP"),
                   c("predicted from HLH-30 binding (PBM)", "HLH-30 binding (PBM)"),
                   c("starvation proteins, CLEAR motif", "CLEAR motif, starvation"),
                   c("S. aureus", "<i>S. aureus</i>")))
      x <- gsub(r[1], r[2], x, fixed = TRUE)
    ital(x)
  }
  # colour = fold vs chance (observed / expected genes, i.e. the cell % over
  # the chance % under the glyph), easier to read than an odds ratio
  hm <- hm %>% mutate(l2or = ifelse(q < 0.05, pmax(pmin(log2(pmax(overlap, 1e-3) / expected), 4), -4), NA),
                      txt = ifelse(!is.na(l2or) & abs(l2or) > 2.5, "white", "grey10"))
  hm$lab <- sub("%", "", hm$lab, fixed = TRUE)
  COLLAB <- sub("%", "", COLLAB, fixed = TRUE)
  hdrp <- hdr + theme(strip.text.x = element_text(size = 10, face = "bold", margin = margin(b = 2)),
                      axis.text.x.top = ggtext::element_markdown(size = 10, lineheight = 1.05),
                      panel.spacing.x = unit(4, "pt"))
  # Row labels are their own panel (name, list size, direction glyph), so a
  # study given in both directions is named once, between its two rows.
  # x runs in inches; y is the row position counted from the bottom of the
  # group, the same positions the heatmap's discrete axis uses.
  rl <- reg %>% mutate(name = short(label), size = format(nl$n[match(key, nl$key)], big.mark = ",", trim = TRUE)) %>%
    group_by(group) %>% mutate(pos = rev(seq_along(key))) %>% ungroup()
  rn <- rl %>% group_by(group, name) %>% summarise(pos = mean(pos), .groups = "drop")
  cairo_pdf(tempfile()); grid::pushViewport(grid::viewport(gp = grid::gpar(fontsize = 10)))
  wid <- function(s) max(sapply(s, function(x) grid::convertWidth(grid::stringWidth(x), "in", TRUE)))
  Wn <- wid(gsub("<[^>]+>", "", rn$name)); Ws <- wid(rl$size); dev.off()
  XG <- Wn + Ws + 0.28; XS <- Wn + Ws + 0.14; LW <- XG + 0.08
  labp <- function(g) {
    d <- rl %>% filter(group == g); dn <- rn %>% filter(group == g)
    ggplot() +
      ggtext::geom_richtext(data = dn, aes(x = Wn, y = pos, label = name), hjust = 1, size = 10 / .pt,
                            colour = "grey10", fill = NA, label.colour = NA, label.padding = unit(0, "pt")) +
      geom_text(data = d, aes(x = XS, y = pos, label = size), hjust = 1, size = 10 / .pt, colour = "grey45") +
      geom_text(data = d %>% filter(!is.na(pdir)) %>% mutate(gc = ifelse(pdir == "up", PAL$up, PAL$down)),
                aes(x = XG, y = pos, label = ifelse(pdir == "up", "▲", "▼"), colour = I(gc)),
                size = 10 / .pt) +
      scale_x_continuous(limits = c(0, LW), expand = c(0, 0)) +
      scale_y_continuous(limits = c(1, nrow(d)), expand = expansion(add = 0.6)) +
      labs(title = ital(g)) + theme_void(base_size = 10) +
      theme(plot.title.position = "plot",
            plot.title = ggtext::element_markdown(size = 10, face = "bold", colour = "grey15", hjust = 0,
                                                  margin = margin(t = 5, b = 2)),
            plot.margin = margin(0, 0, 0, 2)) +
      coord_cartesian(clip = "off")
  }
  heatp <- function(g) pgp(g) + labs(title = NULL) + theme(axis.text.y.left = element_blank())
  # the first group's heading sits in the header row, level with the chance
  # rates, so there is no empty band between the header and the first rows
  blockp <- function(g, first = FALSE) {
    l <- labp(g); if (first) l <- l + labs(title = NULL)
    l + heatp(g) + plot_layout(widths = unit(c(LW, 1), c("in", "null")))
  }
  head1 <- ggplot() +
    ggtext::geom_richtext(aes(x = 0, y = 0, label = ital(G[1])), hjust = 0, vjust = -0.15, size = 10 / .pt,
                          fontface = "bold", colour = "grey15", fill = NA, label.colour = NA,
                          label.padding = unit(0, "pt")) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void(base_size = 10) + theme(plot.margin = margin(0, 0, 0, 2)) + coord_cartesian(clip = "off")
  hdrb <- head1 + hdrp + plot_layout(widths = unit(c(LW, 1), c("in", "null")))
  p03 <- wrap_plots(c(list(hdrb), Map(blockp, G, seq_along(G) == 1)), ncol = 1,
                    heights = c(0.1, nrows + 1.2)) +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  # the fill title has a line break; keep it on one line beside the bar
  p03 <- p03 & scale_fill_gradient2(low = PAL$div[["low"]], mid = "white", high = PAL$div[["high"]],
                                    midpoint = 0, limits = c(-4, 4), breaks = c(-4, -2, 0, 2, 4),
                                    labels = c("16x\nfewer", "4x\nfewer", "chance", "4x\nmore", "16x\nmore"),
                                    name = "vs chance", na.value = "grey95") &
    guides(fill = guide_colourbar(title.position = "left", title.vjust = 0.9,
                                  barwidth = unit(2.6, "in"), barheight = unit(8, "pt")))
  h03 <- 1.15 + 0.205 * (nrow(reg) + 1.2 * length(G) - 1.2)
  ggsave(file.path(pdir, "03-published-lists.pdf"), p03, width = PW, height = h03, device = cairo_pdf)
  ggsave(file.path(pdir, "03-published-lists.png"), p03, width = PW, height = h03, dpi = 300)

  # 04: A = the genes papers show to need hlh-30, B = every core autophagy
  # gene, split over two columns. Role names head their rows (not side
  # strips), so three blocks fit side by side.
  tiles_p <- function(d, tag = NULL) {
    ggplot(d, aes(stage, lab, fill = pmax(pmin(lfc, 1.5), -1.5))) +
      geom_tile(colour = "white", linewidth = 0.6) +
      geom_text(aes(label = sig), size = 10 / .pt, vjust = 0.75) +
      facet_grid(role ~ ., scales = "free_y", space = "free_y") +
      scale_fill_gradient2(low = PAL$div[["low"]], mid = "white", high = PAL$div[["high"]], midpoint = 0,
                           limits = c(-1.5, 1.5), na.value = "grey85", name = expression(log[2] ~ "FC"),
                           breaks = c(-1.5, 0, 1.5), labels = c("-1.5", "0", "1.5")) +
      scale_x_discrete(position = "top") +
      scale_y_discrete(expand = expansion(add = c(0.5, 1.6))) +
      labs(x = NULL, y = NULL, tag = tag) + paper_theme +
      theme(axis.text.y = element_text(face = "italic", size = 10, colour = "grey10"),
            axis.text.x.top = element_text(size = 10, face = "bold", colour = "grey10"),
            # role name as a heading over its rows, left of the tiles
            strip.text.y = element_blank(),
            panel.spacing.y = unit(2, "pt"),
            plot.tag = element_text(size = 10, face = "bold"), plot.tag.position = c(0, 1),
            plot.margin = margin(2, 6, 2, 2))
  }
  # role headings: text in the room left above each facet's first gene,
  # right-aligned to the tiles (clip off), so they take no width
  with_heads <- function(p, d) {
    top <- d %>% group_by(role) %>% summarise(lab = levels(droplevels(lab))[nlevels(droplevels(lab))], .groups = "drop") %>%
      mutate(stage = names(STAGES)[1])
    p + geom_text(data = top, aes(x = 4.5, y = lab, label = role), inherit.aes = FALSE,
                  hjust = 1, vjust = 0.5, nudge_y = 1.15, size = 10 / .pt, fontface = "bold", colour = "grey15") +
      coord_cartesian(clip = "off")
  }
  sub_roles <- function(d, roles) {
    d %>% filter(role %in% roles) %>% mutate(lab = droplevels(lab), role = droplevels(role))
  }
  # 04a and 04b are two figures, each half the page wide: 04a the genes papers
  # show to need hlh-30, 04b every core autophagy gene over two columns
  HW <- 3.75
  units <- function(d) nlevels(droplevels(d$lab)) + 2.3 * nlevels(droplevels(d$role))
  save04 <- function(p, file, um, w = HW, pad = 0.9) {
    h <- pad + 0.19 * um
    ggsave(file.path(pdir, paste0(file, ".pdf")), p, width = w, height = h, device = cairo_pdf)
    ggsave(file.path(pdir, paste0(file, ".png")), p, width = w, height = h, dpi = 300)
    h
  }
  # 04a is one third of the page (21 genes need no more)
  TW <- 2.5
  p04a <- with_heads(tiles_p(d1), d1) + theme(legend.position = "bottom") +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                  barwidth = unit(1.3, "in"), barheight = unit(8, "pt")))
  h04a <- save04(p04a, "04a-autophagy-hlh30-linked", units(d1), w = TW)
  colB1 <- sub_roles(d2, ROLE[1:2]); colB2 <- sub_roles(d2, ROLE[3:5])
  u <- c(units(colB1), units(colB2)); um <- max(u)
  # the colour bar sits under the shorter column, so it adds no height
  k <- which.min(u)
  lay <- c(patchwork::area(1, 3 - k, 2, 3 - k), patchwork::area(1, k, 1, k), patchwork::area(2, k, 2, k))
  cols <- list(with_heads(tiles_p(colB1), colB1), with_heads(tiles_p(colB2), colB2))
  p04b <- wrap_plots(cols[[3 - k]], cols[[k]], guide_area()) +
    plot_layout(design = lay, heights = c(u[k], um - u[k]), guides = "collect") &
    theme(legend.position = "bottom", legend.justification = c(1, 1)) &
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                  barwidth = unit(1.3, "in"), barheight = unit(8, "pt")))
  h04b <- save04(p04b, "04b-autophagy-core", um, pad = 0.4)
  # Sun et al. 2020 (eLife, Supplementary file 1) defines the lysosome genes,
  # as Zhang & Melendez 2025 defines the autophagy genes for 04a/04b
  CAT <- c("Lysosomal membrane proteins" = "membrane proteins", "V0 subunits" = "V-ATPase V0",
           "V1 subunits" = "V-ATPase V1", "Protease (cathepsins)" = "cathepsins",
           "Non-protease hydrolases" = "other hydrolases")
  lys <- read.delim(file.path(opt$pub, "new_lysosome", "sun2020_lysosome_genes.tsv"),
                    stringsAsFactors = FALSE, quote = "") %>% filter(nzchar(gene_id))
  # 04c: the genes of that list that papers show to need hlh-30, in its categories
  b3 <- lnk %>% filter(n_papers_expr_reduced_hlh30 > 0, gene_id %in% lys$gene_id) %>%
    mutate(role_class = lys$category[match(gene_id, lys$gene_id)]) %>%
    transmute(gene_id, role_class, n_pap = n_papers_expr_reduced_hlh30,
              lab = sprintf("%s (%d)", gene_name, n_papers_expr_reduced_hlh30))
  d3 <- tile_df(b3) %>% mutate(role = factor(CAT[role_class], levels = CAT)) %>% droplevels() %>% ord_lab()
  p04c <- with_heads(tiles_p(d3), d3) + theme(legend.position = "bottom") +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                  barwidth = unit(1.3, "in"), barheight = unit(8, "pt")))
  h04c <- save04(p04c, "04c-lysosome-hlh30-linked", units(d3), w = TW)
  message("04c groups: ", paste(names(table(d3$role)), table(d3$role) / 4, collapse = ", "))
  # 04d: every lysosome gene of the Sun list, in its five categories. Two
  # columns would be 11.5 in tall at 10 pt, so the categories go in three
  # full-width columns, the colour bar under the shortest.
  b4 <- lys %>% transmute(gene_id, role_class = category, lab = gene_name)
  d4 <- tile_df(b4) %>% mutate(role = factor(CAT[role_class], levels = CAT)) %>% ord_lab()
  colD <- list(sub_roles(d4, CAT[1:3]), sub_roles(d4, CAT[4]), sub_roles(d4, CAT[5]))
  u <- sapply(colD, units); um <- max(u); k <- which.min(u)
  # layout rows cut at each column's length, so every column keeps the same
  # row height; the colour bar fills the room under the shortest column
  cuts <- sort(unique(u)); nr <- length(cuts)
  lay <- do.call(c, lapply(seq_along(colD), function(j)
    patchwork::area(1, j, match(u[j], cuts), j)))
  lay <- c(lay, patchwork::area(2, k, nr, k))
  p04d <- wrap_plots(c(lapply(colD, function(d) with_heads(tiles_p(d), d)), list(guide_area()))) +
    plot_layout(design = lay, heights = diff(c(0, cuts)), guides = "collect") &
    theme(legend.position = "bottom", legend.justification = c(1, 1)) &
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                  barwidth = unit(1.3, "in"), barheight = unit(8, "pt")))
  h04d <- save04(p04d, "04d-lysosome-core", um, w = PW, pad = 0.4)
  message("04d groups: ", paste(names(table(d4$role)), table(d4$role) / 4, collapse = ", "),
          sprintf("; %.2f x %.2f in", PW, h04d))
  unlink(file.path(pdir, c("04-autophagy-genes.pdf", "04-autophagy-genes.png")))
  message(sprintf("paper: 03 %.2f x %.2f in, 04a %.2f x %.2f in, 04b %.2f x %.2f in, 04c %.2f x %.2f in -> %s",
                  PW, h03, TW, h04a, HW, h04b, TW, h04c, pdir))
}

print(res %>% filter(q < 0.05) %>% arrange(q) %>% select(stage, label, direction, overlap, expected, or, q), row.names = FALSE)
message("wrote ", opt$outdir, "/published_{overlap,autophagy}.{pdf,png}")
