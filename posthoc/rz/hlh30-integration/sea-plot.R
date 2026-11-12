#!/usr/bin/env Rscript
# Plot the SEA results from sea-dge.R as one motif x DEG-set grid, in the
# spirit of Chang et al. 2026 Table 1: rows = motifs significant (SEA
# q < 0.05) in at least one set, columns = stage x direction in the mutant,
# one block per promoter window (500 / 2,000 bp upstream). Motifs whose
# consensus holds an E-box core (CANNTG) are marked, since that is the
# HLH-30 site class.
#
#   Rscript sea-plot.R          # -> out/sea/sea_grid.{pdf,png}, sea_significant.tsv
#   Rscript sea-plot.R --paper  # -> out/sea/paper/01-sea-motif-grid.{pdf,png}:
#                               #    no title or notes, 7.09 in wide, all text 10 pt

suppressPackageStartupMessages({ library(optparse); library(dplyr); library(ggplot2) })
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "palette.R"))
opt <- parse_args(OptionParser(option_list = list(
  make_option("--in",     type = "character", default = "out/sea/sea_all.tsv", dest = "infile"),
  make_option("--top",    type = "integer", default = 25, help = "max motif rows (by best q)"),
  make_option("--deg-dir", type = "character", default = "../volcano/cache/hlh30-all-ruv", dest = "deg_dir"),
  make_option("--norm",   type = "character", default = "cache/normalized_counts.csv",
              help = "DESeq2 normalized counts of the run (results/hlh30-all-ruv/deseq2/)"),
  make_option("--meta",   type = "character", default = "cache/sample_meta.tsv"),
  make_option("--min-base", type = "double", default = 10, dest = "min_base",
              help = "a TF counts as expressed at a stage if its mean normalized count over that stage's samples is at least this"),
  make_option("--outdir", type = "character", default = "out/sea"),
  make_option("--paper",  action = "store_true", default = FALSE,
              help = "paper version only: no title or notes, full page width, 10 pt")
)))
# read every finished SEA run directly (sea_<stage>_<dir>_u<window>/sea.tsv),
# so the grid can be drawn while slower windows are still running
runs <- list.files(dirname(opt$infile), pattern = "^sea_.*_u[0-9]+$", full.names = TRUE)
runs <- runs[file.exists(file.path(runs, "sea.tsv")) & file.size(file.path(runs, "sea.tsv")) > 0]
d <- bind_rows(lapply(runs, function(r) {
  x <- read.delim(file.path(r, "sea.tsv"), comment.char = "#", quote = "", stringsAsFactors = FALSE)
  if (!nrow(x)) return(NULL)
  m <- regmatches(basename(r), regexec("^sea_(.+)_(up|down)_u([0-9]+)$", basename(r)))[[1]]
  # the SEA runs are named by the old short label D0; L4 is shown
  x$stage <- sub("^D0$", "L4", m[2]); x$direction <- m[3]; x$window <- as.integer(m[4]); x
}))
message(length(runs), " SEA runs read")
SETS <- paste(rep(c("L1", "L4", "D3", "D6", "any"), each = 2), c("down", "up"))

# one row per TF name: the best-q motif for that TF in each set (JASPAR and
# CIS-BP often carry the same factor)
# CIS-BP inferred motifs are named "(gene)_(source species)_(DBD_x)"; keep the
# worm gene and flag that the matrix was inferred from another species
d <- d %>% mutate(inferred = grepl("^\\(", ALT_ID),
                  tf = ifelse(is.na(ALT_ID) | ALT_ID == "", ID, sub("^\\(([^)]+)\\).*", "\\1", ALT_ID)),
                  tf = ifelse(inferred, paste0(tf, " (inf.)"), tf),
                  set_lab = factor(paste(stage, direction), levels = SETS),
                  ebox = grepl("CA..TG", CONSENSUS))
# Expression filter (as Chang et al. 2026 did for their candidates): a motif
# hit only counts if the TF's own gene is expressed at that stage in our data
# (DESeq2 baseMean >= min_base); for the pooled "any" columns, at any stage.
# Stage-specific expression: mean DESeq2-normalized count over every sample of
# that stage (both genotypes). The baseMean column of the contrast tables is
# NOT usable here: it averages all samples of the run, so it is the same at
# every stage.
STG <- c(L1 = "L1", L4 = "L4", D3 = "Day3", D6 = "Day6")
nc   <- read.csv(opt$norm, row.names = 1, check.names = FALSE)
meta <- read.delim(opt$meta, stringsAsFactors = FALSE)
nm <- read.csv(file.path(opt$deg_dir, "tm1978_L1_vs_WT_L1.csv"), stringsAsFactors = FALSE)[, c("gene_id", "gene_name")]
expr <- do.call(rbind, lapply(names(STG), function(st) {
  smp <- intersect(meta$sample_id[meta$stage == STG[[st]]], colnames(nc))
  m <- rowMeans(nc[, smp, drop = FALSE])
  data.frame(stage = st, gene = tolower(nm$gene_name[match(names(m), nm$gene_id)]), mean = m,
             on = m >= opt$min_base)
}))
expr <- expr[!is.na(expr$gene), ]
is_on <- function(gene, stage) {
  g <- tolower(sub(" \\(inf\\.\\)$", "", gene))
  mapply(function(gg, st) {
    e <- expr[expr$gene == gg & (st == "any" | expr$stage == st), "on"]
    length(e) > 0 && any(e)
  }, g, stage)
}
d$expressed <- is_on(d$tf, d$stage)
message(sprintf("%d of %d significant hits are for TFs not expressed at that stage",
                sum(d$QVALUE < 0.05 & !d$expressed), sum(d$QVALUE < 0.05)))
# rows: top factors by best q among hits that pass the expression filter
keep <- d %>% filter(QVALUE < 0.05, expressed) %>% group_by(tf) %>% summarise(best = min(QVALUE)) %>%
  arrange(best) %>% slice_head(n = opt$top)
g <- d %>% filter(tf %in% keep$tf) %>% group_by(tf, set_lab, window) %>%
  slice_min(QVALUE, n = 1, with_ties = FALSE) %>% ungroup() %>%
  mutate(sig = QVALUE < 0.05 & expressed,
         l2e = ifelse(sig, pmin(log2(ENR_RATIO), 3), NA),
         star = case_when(!expressed ~ "n.e.", QVALUE < 0.001 ~ "***", QVALUE < 0.01 ~ "**",
                          QVALUE < 0.05 ~ "*", TRUE ~ ""),
         window = factor(paste0(window, " bp upstream"),
                         levels = intersect(c("500 bp upstream", "2000 bp upstream"), paste0(sort(unique(window)), " bp upstream"))))
ordr <- g %>% group_by(tf) %>%
  summarise(e = any(ebox), n = sum(sig), b = suppressWarnings(min(QVALUE[sig]))) %>% arrange(e, n, -b)
g$tf <- factor(g$tf, levels = ordr$tf)
# the factors in the grid with the motif behind their best hit
write.table(d %>% filter(tf %in% keep$tf, QVALUE < 0.05, expressed) %>% group_by(tf) %>%
              slice_min(QVALUE, n = 1, with_ties = FALSE) %>% ungroup() %>%
              transmute(tf, DB, ID, CONSENSUS, best_q = QVALUE, ebox) %>% arrange(best_q),
            file.path(opt$outdir, "sea_grid_factors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(d %>% filter(QVALUE < 0.05) %>% arrange(QVALUE) %>%
              select(set = set_lab, window, DB, ID, tf, CONSENSUS, TP, TP., FP, FP., ENR_RATIO, QVALUE, expressed),
            file.path(opt$outdir, "sea_significant.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# every row x column x window, flagged where the TF itself is not expressed at
# that stage, whether or not SEA reported the motif there
full <- expand.grid(tf = levels(g$tf), set_lab = factor(SETS, levels = SETS), window = levels(g$window),
                    stringsAsFactors = FALSE) %>%
  mutate(stage = sub(" .*", "", set_lab), expressed = is_on(tf, stage),
         tf = factor(tf, levels = levels(g$tf)), set_lab = factor(set_lab, levels = SETS),
         window = factor(window, levels = levels(g$window)))
ne <- full %>% filter(!expressed) %>%
  mutate(tf = factor(tf, levels = levels(g$tf)), set_lab = factor(set_lab, levels = SETS),
         window = factor(window, levels = levels(g$window)))
# diagonal hatching for those cells, drawn as segments clipped to the cell box
# (x, y are the discrete positions; the box matches the outline tiles)
HW <- 0.46; HH <- 0.43; HSTEP <- 0.25
hatch <- if (nrow(ne)) do.call(rbind, lapply(seq_len(nrow(ne)), function(i) {
  cc <- seq(-1 + HSTEP / 2, 1 - HSTEP / 2, by = HSTEP)   # line u - v = c in the unit box
  data.frame(set_lab = ne$set_lab[i], tf = ne$tf[i], window = ne$window[i],
             u0 = pmax(0, cc), v0 = pmax(0, -cc), u1 = pmin(1, 1 + cc), v1 = pmin(1, 1 - cc))
})) else NULL
if (!is.null(hatch)) hatch <- hatch %>%
  mutate(x = as.integer(set_lab), y = as.integer(tf),
         x0 = x - HW + 2 * HW * u0, x1 = x - HW + 2 * HW * u1,
         y0 = y - HH + 2 * HH * v0, y1 = y - HH + 2 * HH * v1)
ghost <- g %>% filter(!expressed, QVALUE < 0.05) %>%
  mutate(star = case_when(QVALUE < 0.001 ~ "***", QVALUE < 0.01 ~ "**", TRUE ~ "*"))
g <- g %>% filter(expressed)
ytxt <- setNames(ifelse(ordr$e, paste0(ordr$tf, "  [E-box]"), ordr$tf), ordr$tf)
grid_layers <- function(txt_size) list(
  geom_tile(colour = "white", linewidth = 0.6),
  # every cell gets a barely-there outline so empty cells still read as a grid
  geom_tile(data = full, aes(set_lab, tf), inherit.aes = FALSE, fill = NA,
            colour = "grey90", linewidth = 0.25, width = 2 * HW, height = 2 * HH),
  # diagonal hatching = the TF is not expressed at that stage in our data
  geom_segment(data = hatch, aes(x = x0, y = y0, xend = x1, yend = y1), inherit.aes = FALSE,
               colour = "grey75", linewidth = 0.25),
  geom_text(aes(label = star), size = txt_size, vjust = 0.75, colour = "grey10"),
  # the hits the filter removed, in grey inside their outline
  geom_text(data = ghost, aes(set_lab, tf, label = star), inherit.aes = FALSE, size = txt_size, vjust = 0.75,
            colour = "grey60"))
fill_scale <- function(name) scale_fill_gradient(low = "#F3E3E0", high = PAL$up, na.value = "white", limits = c(0, 3),
                                                 breaks = 0:3, labels = c("1", "2", "4", "8+"), name = name)

if (opt$paper) {
  # one figure for the paper: 10 pt everywhere at the printed width
  pdir <- file.path(opt$outdir, "paper"); dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  win_lab <- setNames(sub("2000", "2,000", sub(" upstream$", " upstream", levels(g$window))), levels(g$window))
  # rows name the factor that binds the motif, so protein style: upper case, roman (NHR-25)
  prot <- sub("^(\\S+)", "\\U\\1", ordr$tf, perl = TRUE)
  ytxt_p <- setNames(ifelse(ordr$e, paste0(prot, " (E-box)"), prot), ordr$tf)
  pp <- ggplot(g, aes(set_lab, tf, fill = l2e)) + grid_layers(10 / .pt) +
    facet_wrap(~ window, nrow = 1, labeller = as_labeller(win_lab)) +
    fill_scale("fold enrichment") +
    scale_y_discrete(labels = ytxt_p) +
    labs(x = expression("DEGs of" ~ italic("hlh-30(tm1978)") ~ "vs. WT"), y = NULL) +
    guides(fill = guide_colourbar(title.vjust = 0.8, barwidth = unit(1.2, "in"), barheight = unit(0.12, "in"))) +
    theme_minimal(base_size = 10) +
    theme(text = element_text(size = 10), axis.text = element_text(size = 10, colour = "grey15"),
          axis.title = element_text(size = 10), strip.text = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 10), legend.title = element_text(size = 10),
          panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "bottom", legend.justification = "right",
          plot.margin = margin(4, 4, 4, 4))
  hp <- 1.9 + 0.27 * nrow(ordr)
  ggsave(file.path(pdir, "01-sea-motif-grid.pdf"), pp, width = 7.09, height = hp, device = cairo_pdf)
  ggsave(file.path(pdir, "01-sea-motif-grid.png"), pp, width = 7.09, height = hp, dpi = 300)
  message("wrote ", pdir, "/01-sea-motif-grid.{pdf,png}; ", nrow(ordr), " factors, ", round(hp, 2), " in tall")
  quit(save = "no")
}

p <- ggplot(g, aes(set_lab, tf, fill = l2e)) + grid_layers(2.8) +
  facet_wrap(~ window, nrow = 1) +
  fill_scale("enrichment\n(q < 0.05)") +
  scale_y_discrete(labels = ytxt) +
  labs(x = "promoters of DEGs, hlh-30(tm1978) vs WT  (stage, direction in the mutant)", y = NULL,
       title = "Motifs enriched in promoters of hlh-30 DEGs\n(SEA, all JASPAR nematode + CIS-BP C. elegans motifs)",
       subtitle = paste0("vs promoters of expressed non-DEGs at the same stage; * q < 0.05, ** < 0.01, *** < 0.001; empty box = not significant.\n",
                         "Hatched = the TF itself is not expressed at that stage (mean normalized count < ", opt$min_base,
                         "); grey stars inside = hits removed by that filter.\n'any' = union over stages. ",
                         "Top ", nrow(ordr), " factors by best q; (inf.) = CIS-BP matrix inferred from another species.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"))
h <- 2.4 + 0.2 * nrow(ordr); w <- 3.2 + 3.4 * nlevels(g$window)
ggsave(file.path(opt$outdir, "sea_grid.pdf"), p, width = w, height = h, device = cairo_pdf)
ggsave(file.path(opt$outdir, "sea_grid.png"), p, width = w, height = h, dpi = 150)
message("wrote ", opt$outdir, "/sea_grid.{pdf,png}; ", nrow(ordr), " factors")
