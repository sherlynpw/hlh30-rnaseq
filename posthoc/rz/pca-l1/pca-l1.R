#!/usr/bin/env Rscript
# L1-only PCA (WT vs hlh-30(tm1978)) in the posthoc pca-stage style.
# Colour encodes stage (all L1 green) and shape encodes genotype, the rule of
# every other PCA panel; one enclosing ellipse per genotype.

suppressPackageStartupMessages({
  library(optparse); library(ggplot2); library(ggforce); library(dplyr)
})

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "palette.R"))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "enclose.R"))


opt <- parse_args(OptionParser(option_list = list(
  make_option("--coords",   type = "character", default = "cache/pca_coordinates.csv"),
  make_option("--variance", type = "character", default = "cache/variance_explained.txt"),
  make_option("--font-scale", type = "double", default = 1, dest = "font_scale",
              help = "Multiply every font size; render at the ORIGINAL canvas."),
  make_option("--legend",   type = "character", default = "top",
              help = "top (inside, top centre), bottom (under the plot) or genotype (genotype key only, untitled, inside the top)"),
  make_option("--height",   type = "double", default = 5.5),
  make_option("--outdir",   type = "character", default = "out")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

var_lines <- readLines(opt$variance)
pct <- function(pc) {
  m <- regmatches(var_lines, regexec(paste0("^", pc, ":\\s*([0-9.]+)%"), var_lines))
  as.numeric(Filter(function(x) length(x) == 2, m)[[1]][2])
}
pc1_pct <- pct("PC1"); pc2_pct <- pct("PC2")

df <- read.csv(opt$coords, stringsAsFactors = FALSE, check.names = FALSE)
# accept both new (WT/tm1978) and legacy (N2/hlh-30) tokens
geno_chr  <- dplyr::recode(df$genotype, "N2" = "WT", "WT" = "WT",
                           "hlh-30" = "hlh-30(tm1978)", "tm1978" = "hlh-30(tm1978)")
geno_lvls <- intersect(c("WT", "hlh-30(tm1978)"), unique(geno_chr))
df$geno_label <- factor(geno_chr, levels = geno_lvls)

shape_geno   <- c("WT" = 16, "hlh-30(tm1978)" = 17)
# gene name italic; WT stays roman
geno_lab_expr <- parse(text = vapply(geno_lvls, function(g)
  if (g == "WT") sprintf('"%s"', g) else sprintf('italic("%s")', g), character(1)))

p <- ggplot(df, aes(x = PC1, y = PC2)) +
  # smallest ellipse around each genotype's replicates, as in the other PCA panels
  enclosing_outline(df, "PC1", "PC2", "geno_label", 7.5, opt$height, opt$font_scale,
                    pad_in = 0.22, legend_in = if (opt$legend == "bottom") 0.9 * opt$font_scale else 0) +
  geom_point(aes(color = "L1", shape = geno_label), size = 3, alpha = 0.95) +
  (if (opt$legend == "genotype") guides(color = "none", shape = guide_legend(title = NULL))) +
  scale_color_manual(values = PAL$stage["L1"], name = "Stage") +
  scale_shape_manual(values = shape_geno, name = "Genotype", drop = FALSE,
                     labels = geno_lab_expr) +
  # Stage above Genotype, as in pca-stage / pca-genotype
  (if (opt$legend != "genotype") guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))) +
  # room for the enclosing ellipses, which reach past the outer replicates
  coord_cartesian(xlim = range(df$PC1) + c(-6, 6),
                  ylim = c(min(df$PC2) - 12, max(df$PC2) + if (opt$legend == "bottom") 12 else 22), clip = "off") +
  labs(x = sprintf("PC1 (%.2f%%)", pc1_pct), y = sprintf("PC2 (%.2f%%)", pc2_pct)) +
  theme_bw(base_size = 12 * opt$font_scale) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12 * opt$font_scale),
        # top centre is the only empty region: hlh-30 left, WT right
        # genotype key: top left, above the hlh-30 group (the WT group reaches the top right)
        # genotype key: one line across the top
        legend.position = if (opt$legend == "bottom") "bottom" else c(0.5, 0.99),
        legend.justification = if (opt$legend == "bottom") "center" else c(0.5, 1),
        legend.title = element_text(size = 9.6 * opt$font_scale),   # same 10pt as its keys
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        legend.box = "horizontal", legend.box.just = "top",
        legend.direction = if (opt$legend == "genotype") "horizontal" else "vertical",
        legend.key = element_blank())

ggsave(file.path(opt$outdir, "pca_l1_genotype.pdf"), p, width = 7.5, height = opt$height)   # same canvas width as the other PCAs
ggsave(file.path(opt$outdir, "pca_l1_genotype.png"), p, width = 7.5, height = opt$height, dpi = 200)
message("wrote pca_l1_genotype.{pdf,png} to ", opt$outdir)
