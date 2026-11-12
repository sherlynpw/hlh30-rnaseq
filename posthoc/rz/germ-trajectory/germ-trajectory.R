#!/usr/bin/env Rscript
# Expression trajectory of the germ-cell genes that come out of the compare-DE
# analysis (genes in germ-cell-related enriched terms). For each gene, show its
# level across L1 -> L4 -> Day3 -> Day6 in WT vs hlh-30(tm1978), to test whether
# the two genotypes start together at L4 and diverge later, or differ at baseline.
#
# Inputs: DESeq2 VST counts (WBGene x sample), sample_meta, a gene-symbol list.

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(tidyr); library(readr); library(stringr)
  library(ggplot2); library(RSQLite)
})

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "palette.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--vst",     type = "character", default = "cache/vst.csv"),
  make_option("--meta",    type = "character", default = "cache/meta.tsv"),
  make_option("--genes",   type = "character", default = "cache/germ_genes.txt"),
  make_option("--gene-db", type = "character", dest = "gene_db", default = "../../../gene_db.sqlite"),
  make_option("--font-scale", type = "double", default = 1, dest = "font_scale",
              help = "Multiply every font size; render at the ORIGINAL canvas."),
  make_option("--outdir",  type = "character", default = "out"),
  make_option("--title",   type = "character", default = "Germ-cell genes",
              help = "gene set name, the title of the mean-line panel"),
  make_option("--prefix",  type = "character", default = "germ", help = "output file prefix"),
  make_option("--ylim",    type = "character", default = "",
              help = "mean-line y range as lo,hi, so paired panels share an axis (empty = auto)"),
  make_option("--width",   type = "double", default = 6.2, help = "mean-line canvas width (in)"),
  make_option("--height",  type = "double", default = 4.4, help = "mean-line canvas height (in)"),
  make_option("--hm-width", type = "double", default = 6.8, dest = "hm_width", help = "heatmap canvas width (in)"),
  make_option("--row-in",  type = "double", default = 0.2, dest = "row_in",
              help = "heatmap height per gene row (in), so paired heatmaps have equal rows"),
  make_option("--zlim",    type = "double", default = NA,
              help = "heatmap colour limit, so paired heatmaps share a scale (NA = auto)")
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
ACC <- PAL$hlh30; WTc <- PAL$wt; BLU <- PAL$div[["low"]]; RED <- PAL$div[["high"]]
STAGES <- c("L1","L4","Day3","Day6")
STAGE_DISP <- c(L1 = "L1", L4 = "L4", Day3 = "D3", Day6 = "D6")

# --- map gene symbols -> WBGene ----------------------------------------------
syms <- readLines(opt$genes); syms <- syms[nzchar(syms)]
con <- dbConnect(SQLite(), opt$gene_db); on.exit(dbDisconnect(con))
ann <- dbGetQuery(con, "SELECT wormbase_gene, gene_name, sequence_name FROM annotations")
sym2wb <- c(setNames(ann$wormbase_gene, ann$gene_name), setNames(ann$wormbase_gene, ann$sequence_name))
sym2wb <- sym2wb[!is.na(names(sym2wb)) & nzchar(names(sym2wb))]
wb2sym <- setNames(ifelse(!is.na(ann$gene_name) & nzchar(ann$gene_name), ann$gene_name, ann$sequence_name),
                   ann$wormbase_gene)
wb <- unname(sym2wb[syms]); names(wb) <- syms
wb <- wb[!is.na(wb)]
message(length(wb), "/", length(syms), " gene symbols mapped to WBGene")

# --- VST + meta --------------------------------------------------------------
vst  <- read.csv(opt$vst, row.names = 1, check.names = FALSE)
meta <- read.delim(opt$meta, stringsAsFactors = FALSE)
meta <- meta %>% mutate(stage = factor(stage, levels = STAGES),
                        genotype = factor(genotype, levels = c("WT","tm1978")))
keep_wb <- intersect(unname(wb), rownames(vst))
message(length(keep_wb), " genes present in the VST matrix")

long <- vst[keep_wb, , drop = FALSE] %>%
  tibble::rownames_to_column("wb") %>%
  pivot_longer(-wb, names_to = "sample_id", values_to = "vst") %>%
  inner_join(meta %>% select(sample_id, genotype, stage), by = "sample_id") %>%
  mutate(gene = wb2sym[wb])

# mean per gene x genotype x stage, then z-score per gene across the group means
grp <- long %>% group_by(gene, wb, genotype, stage) %>%
  summarise(mean_vst = mean(vst), .groups = "drop")
grp <- grp %>% group_by(gene) %>%
  mutate(z = as.numeric(scale(mean_vst))) %>% ungroup()

# =============================================================== heatmap ======
mat <- grp %>% mutate(col = paste(genotype, stage, sep = "_")) %>%
  select(gene, col, z) %>% pivot_wider(names_from = col, values_from = z) %>%
  tibble::column_to_rownames("gene")
gene_ord <- rownames(mat)[hclust(dist(as.matrix(mat)))$order]

hm <- grp %>% mutate(gene = factor(gene, levels = gene_ord))
geno_lab <- as_labeller(c(WT = "WT", tm1978 = "italic('hlh-30(tm1978)')"), default = label_parsed)
zlim <- if (is.na(opt$zlim)) max(abs(hm$z), na.rm = TRUE) else opt$zlim
ph <- ggplot(hm, aes(stage, gene, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  facet_grid(~ genotype, labeller = geno_lab) +
  scale_x_discrete(labels = STAGE_DISP) +
  scale_fill_gradient2(low = BLU, mid = "white", high = RED, midpoint = 0,
                       limits = c(-zlim, zlim), oob = scales::squish,
                       name = "z-score per gene") +
  labs(x = NULL, y = NULL, title = sprintf("%s (%d)", opt$title, nrow(mat)),
       subtitle = "VST expression, rows clustered") +
  theme_minimal(base_size = 8 * opt$font_scale) +
  theme(plot.title.position = "plot") +
  theme(panel.grid = element_blank(),
        panel.spacing.x = unit(5, "pt"),
        strip.text = element_text(face = "bold", size = 8 * opt$font_scale, colour = "grey15"),
        axis.text.x = element_text(size = 8 * opt$font_scale),
        # same sans italic as the other panels' gene names (was mono)
        axis.text.y = element_text(size = 8 * opt$font_scale, face = "italic"),
        # title and subtitle as in the mean-line panels (12.5pt, 9pt in print)
        plot.title = element_text(face = "bold", size = 10 * opt$font_scale),
        plot.subtitle = element_text(size = 7.2 * opt$font_scale, colour = "grey40"),
        legend.position = "bottom",
        legend.title = element_text(size = 8 * opt$font_scale, vjust = 0.8),
        legend.key.height = unit(8, "pt"), legend.key.width = unit(22, "pt"))
# fixed space for the title, strips, axis and legend, plus one row per gene,
# so two heatmaps with the same --row-in have rows of the same height
hm_h <- 1.75 + opt$row_in * nrow(mat)
ggsave(file.path(opt$outdir, paste0(opt$prefix, "_heatmap.pdf")), ph, width = opt$hm_width, height = hm_h)
ggsave(file.path(opt$outdir, paste0(opt$prefix, "_heatmap.png")), ph, width = opt$hm_width, height = hm_h, dpi = 200)

# ====================================================== group-mean trajectory =
gm <- grp %>% group_by(genotype, stage) %>%
  summarise(m = mean(z), se = sd(z)/sqrt(n()), .groups = "drop")
pg <- ggplot(gm, aes(stage, m, colour = genotype, group = genotype, fill = genotype)) +
  geom_ribbon(aes(ymin = m - se, ymax = m + se), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + geom_point(size = 2.6) +
  scale_x_discrete(labels = STAGE_DISP) +
  scale_colour_manual(values = c(WT = WTc, tm1978 = ACC), name = NULL,
                      labels = c(WT = "WT", tm1978 = expression(italic("hlh-30(tm1978)")))) +
  scale_fill_manual(values = c(WT = WTc, tm1978 = ACC), guide = "none") +
  labs(x = NULL, y = "mean z-scored\nexpression",
       title = sprintf("%s (%d)", opt$title, nrow(mat)),
       subtitle = "mean \u00b1 SEM across genes") +
  scale_y_continuous(breaks = scales::breaks_width(0.5)) +
  (if (nzchar(opt$ylim)) coord_cartesian(ylim = as.numeric(strsplit(opt$ylim, ",")[[1]])) else NULL) +
  theme_minimal(base_size = 12 * opt$font_scale) +
  theme(plot.title.position = "plot", plot.caption.position = "plot") +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12 * opt$font_scale),
        # 9pt in print: axis text is 9.6 * font_scale -> 10pt, so 9pt is 8.6x
        plot.subtitle = element_text(size = 8.6 * opt$font_scale, colour = "grey40"),
        plot.caption = element_text(size = 8.6 * opt$font_scale, colour = "grey55", hjust = 0),
        legend.position = "top", legend.justification = "left")
ggsave(file.path(opt$outdir, paste0(opt$prefix, "_meanline.pdf")), pg, width = opt$width, height = opt$height)
ggsave(file.path(opt$outdir, paste0(opt$prefix, "_meanline.png")), pg, width = opt$width, height = opt$height, dpi = 200)

# tidy table out
write_tsv(grp %>% arrange(gene, genotype, stage), file.path(opt$outdir, paste0(opt$prefix, "_trajectory.tsv")))
message("wrote ", opt$prefix, "_{heatmap,meanline}.{pdf,png} and ", opt$prefix, "_trajectory.tsv to ", opt$outdir)
