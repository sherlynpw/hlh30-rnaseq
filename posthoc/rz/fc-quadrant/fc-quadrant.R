#!/usr/bin/env Rscript
# Dot-plot companion to the DEG Venns: for one stage transition, each gene
# that is a DEG in WT or in hlh-30(tm1978) is a dot at (log2FC in WT, log2FC
# in hlh-30). The four quadrants show whether the two genotypes move genes the
# same way (top-right, bottom-left) or opposite ways (top-left, bottom-right).
# Colour = Venn class; a count sits in every quadrant corner.
#
# Same contrasts, orientation and DEG cut-offs as ../venn-degs/venn-degs.R;
# the DESeq2 tables come from that script's cache.

suppressPackageStartupMessages({
  library(optparse); library(ggplot2); library(dplyr); library(readr)
})
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "palette.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--cache",  type = "character", default = "../venn-degs/cache"),
  make_option("--outdir", type = "character", default = "out"),
  make_option("--padj",   type = "double",    default = 0.05),
  make_option("--lfc",    type = "double",    default = 1.0),
  make_option("--font-scale", type = "double", default = 1, dest = "font_scale",
              help = "Multiply every font size. Render at the ORIGINAL canvas and raise this."),
  make_option("--l4-name", type = "character", default = "L4", dest = "l4_name",
              help = "Short label for the L4 stage (paper panels: L4)."),
  make_option("--gtf",    type = "character", default = "../chip-tracks/cache/gtf.gtf",
              help = "Ensembl GTF, for sequence names of genes with no CGC name (optional)"),
  make_option("--legend-ids", type = "character", default = "", dest = "legend_ids",
              help = "comma-separated transition ids (A-D) that keep the legend; empty = all"),
  make_option("--legend-ncol", type = "integer", default = 2, dest = "legend_ncol",
              help = "legend columns; 1 for one-third width"),
  make_option("--width",  type = "double",    default = 4.6),
  make_option("--height", type = "double",    default = 4.6)
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

TR <- list(
  list(id="A", title=paste("L1 →", opt$l4_name), run="hlh30-all-ruv",
       wt="WT_L1_vs_WT_L4", wt_flip=TRUE, mut="tm1978_L4_vs_tm1978_L1", mut_flip=FALSE),
  list(id="B", title=paste(opt$l4_name, "→ D3"), run="VZRJXQ-ruv-rbe",
       wt="N2_Day3_vs_N2_L4", wt_flip=FALSE, mut="hlh-30_L4_vs_hlh-30_Day3", mut_flip=TRUE),
  list(id="C", title=paste(opt$l4_name, "→ D6"), run="VZRJXQ-ruv-rbe",
       wt="N2_Day6_vs_N2_L4", wt_flip=FALSE, mut="hlh-30_L4_vs_hlh-30_Day6", mut_flip=TRUE),
  list(id="D", title="D3 → D6", run="VZRJXQ-ruv-rbe",
       wt="N2_Day6_vs_N2_Day3", wt_flip=FALSE, mut="hlh-30_Day6_vs_hlh-30_Day3", mut_flip=FALSE)
)

load_fc <- function(run, name, flip) {
  d <- read.csv(file.path(opt$cache, paste0(run, "__", name, ".csv")), stringsAsFactors = FALSE)
  # genes with no padj (DESeq2 independent filtering) drop out here, and the
  # inner join below keeps only genes tested in both genotypes, as the Venn does
  d <- d[!is.na(d$padj) & !is.na(d$log2FoldChange), ]
  if (flip) d$log2FoldChange <- -d$log2FoldChange
  d$deg <- d$padj < opt$padj & abs(d$log2FoldChange) >= opt$lfc
  d[, c("gene_id", "gene_name", "concise_description", "log2FoldChange", "padj", "deg")]
}

# Colour = where the gene is significant, which is what the Venn shows.
# Opposite = a DEG in both genotypes with opposite signs (the Venn's rule). A
# DEG in one genotype only stays in that genotype's class, whatever the sign
# in the other. Non-DE genes are not drawn (15k of them inside the cut-off
# box pack into a solid grey patch at any size or opacity); they are counted
# in the quadrant totals.
CLASS_LEVELS <- c("DE in both", "DE in WT only", "DE in hlh-30 only", "opposite")
# shared genes light grey, so the genotype-specific genes stand out; WT grey
# and HLH-30 purple as in every panel; opposite in a bright orange (Okabe-Ito)
CLASS_COL <- c("DE in both" = "#E6E6E6", "DE in WT only" = "#8C8C8C",
               "DE in hlh-30 only" = PAL$hlh30, "opposite" = "#E69F00")
CLASS_LAB <- c("DE in both", "DE in WT only", expression("DE in" ~ italic("hlh-30") ~ "only"),
               "opposite (DE in both)")

draw <- function(t) {
  wt  <- load_fc(t$run, t$wt,  t$wt_flip)
  mut <- load_fc(t$run, t$mut, t$mut_flip)
  d <- inner_join(wt, mut %>% select(-gene_name, -concise_description), by = "gene_id", suffix = c("_wt", "_mut")) %>%
    mutate(deg = deg_wt | deg_mut,
           opp = deg_wt & deg_mut & sign(log2FoldChange_wt) != sign(log2FoldChange_mut),
           class = case_when(
             opp ~ "opposite",
             deg_wt & deg_mut ~ "DE in both",
             deg_wt ~ "DE in WT only",
             TRUE ~ "DE in hlh-30 only"),
           class = factor(class, levels = CLASS_LEVELS)) %>%
    arrange(deg, class)                            # shared genes first, the rest on top

  lim <- max(abs(c(d$log2FoldChange_wt[d$deg], d$log2FoldChange_mut[d$deg]))) * 1.08
  # quadrant labels "DEGs / all genes" by sign of the two fold changes, set
  # a little in from the corners
  q <- d %>%
    mutate(qx = log2FoldChange_wt > 0, qy = log2FoldChange_mut > 0) %>%
    group_by(qx, qy) %>% summarise(n_deg = sum(deg), n_all = n(), .groups = "drop") %>%
    mutate(pct = 100 * n_deg / n_all,
           # one decimal below 1%, so 11 / 2312 is not "0%"
           label = sprintf("%d / %d\n(%s%%)", n_deg, n_all,
                           ifelse(pct < 1, sprintf("%.1f", pct), sprintf("%.0f", pct))),
           x = ifelse(qx, lim, -lim) * 0.9, y = ifelse(qy, lim, -lim) * 0.9,
           hjust = ifelse(qx, 1, 0), vjust = ifelse(qy, 1, 0))

  mut_expr <- 'italic("hlh-30(tm1978)")'
  p <- ggplot(d, aes(log2FoldChange_wt, log2FoldChange_mut)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, colour = "grey85", linetype = "22", linewidth = 0.35) +
    geom_point(data = function(x) x[x$deg & !x$opp, ], aes(colour = class), size = 0.9, alpha = 0.9,
               show.legend = FALSE) +
    # opposite genes on top and a little larger, since there are only a few
    geom_point(data = function(x) x[x$opp, ], aes(colour = class), size = 1.6, show.legend = FALSE) +
    # invisible layer that carries the legend keys; one row per class, so a
    # class with no genes in this panel still gets its key
    geom_point(data = data.frame(log2FoldChange_wt = 0, log2FoldChange_mut = 0,
                                 class = factor(CLASS_LEVELS, levels = CLASS_LEVELS)),
               aes(colour = class), alpha = 0) +
    geom_text(data = q, aes(x, y, label = label, hjust = hjust, vjust = vjust),
              size = 3.4 * opt$font_scale, colour = "grey15", lineheight = 0.9) +
    scale_colour_manual(values = CLASS_COL, labels = CLASS_LAB, name = NULL, drop = FALSE) +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    labs(x = expression(log[2] ~ "FC, WT"),
         y = expression(log[2] ~ "FC," ~ italic("hlh-30(tm1978)")),
         subtitle = t$title) +
    theme_bw(base_size = 12 * opt$font_scale) +
    theme(panel.grid = element_blank(),
          plot.subtitle = element_text(face = "bold", size = 12 * opt$font_scale),
          legend.position = "bottom",
          legend.text = element_text(size = 9.6 * opt$font_scale),
          legend.key.spacing.x = unit(6 * opt$font_scale, "pt")) +
    guides(colour = guide_legend(ncol = opt$legend_ncol, override.aes = list(size = 2.5, alpha = 1)))

  # one legend per figure: panels not listed in --legend-ids drop it and
  # lose its height, so the plot keeps its size
  keep_ids <- trimws(strsplit(opt$legend_ids, ",")[[1]])
  h <- opt$height
  if (length(keep_ids) && !(t$id %in% keep_ids)) {
    p <- p + theme(legend.position = "none")
    h <- opt$height - 0.75 * opt$font_scale
  }
  base <- file.path(opt$outdir, sprintf("quadrant_%s_%s", t$id,
                                        gsub("[^A-Za-z0-9]+", "_", tolower(t$title))))
  ggsave(paste0(base, ".pdf"), p, width = opt$width, height = h, device = cairo_pdf)
  ggsave(paste0(base, ".png"), p, width = opt$width, height = h, dpi = 200)
  message("[wrote] ", base, ".{pdf,png}  DEGs=", sum(d$deg), "  ", paste(names(table(d$class[d$deg])), table(d$class[d$deg]), sep = "=", collapse = " "))
  invisible(d)
}
# the opposite genes of every transition, one row per gene and transition
opp <- bind_rows(lapply(TR, function(t) {
  draw(t) %>% filter(opp) %>%
    transmute(transition = sub("→", "to", t$title), gene_id, gene_name,
              log2FC_WT = round(log2FoldChange_wt, 3), padj_WT = signif(padj_wt, 3),
              log2FC_hlh30 = round(log2FoldChange_mut, 3), padj_hlh30 = signif(padj_mut, 3),
              description = concise_description)
}))
# genes with no CGC name carry their WormBase id as the name; use the
# sequence name from the GTF instead (C33A11.2), which readers can search
if (file.exists(opt$gtf) && nrow(opp)) {
  gl <- grep("\tgene\t", readLines(opt$gtf), value = TRUE)
  seqname <- setNames(sub('.*gene_name "([^"]+)".*', "\\1", gl), sub('.*gene_id "([^"]+)".*', "\\1", gl))
  fix <- opp$gene_name == opp$gene_id & opp$gene_id %in% names(seqname)
  opp$gene_name[fix] <- seqname[opp$gene_id[fix]]
}
write.csv(opp, file.path(opt$outdir, "opposite_genes.csv"), row.names = FALSE, na = "")
message("[wrote] ", file.path(opt$outdir, "opposite_genes.csv"), "  ", nrow(opp), " genes")
