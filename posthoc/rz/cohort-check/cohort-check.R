#!/usr/bin/env Rscript
# Are the two adult cohorts (room temperature, 22 C, and 20 C) different
# enough to matter? The paper model of the 28-sample set is
# ~ W_1 + replicate + condition, with no cohort term, and the replicate
# numbers 1 to 3 repeat in both cohorts, so replicate does not absorb cohort.
#
#   M0  ~ W_1 + replicate + condition              the paper model (28/84/253 DEGs; the files say 29/84/251)
#   M1  ~ W_1 + replicate + RT + condition         adds the cohort (RT = 1 for a room-temperature adult)
#   M1b ~ replicate + RT + condition               the cohort effect without RUV
#   M2  ~ W_1 + replicate + RT + RT_mut + condition RT_mut = RT for hlh-30 adults: a genotype effect
#                                                   that differs between cohorts
# plus a PCA of the 22 adult samples coloured by cohort, before and after
# removing W_1 and replicate (the correction of Fig 1D).
#
# Runs in the pipeline container (DESeq2 1.48.2, ashr 2.2.63), as the paper
# models did. DEG: padj < 0.05 and |ashr-shrunken log2FC| >= 1.

suppressPackageStartupMessages({
  library(optparse); library(DESeq2); library(limma); library(ggplot2)
})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--cache",   type = "character", default = "cache"),
  make_option("--outdir",  type = "character", default = "out"),
  make_option("--palette", type = "character", default = "../_common/palette.R")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
source(opt$palette)
PADJ <- 0.05; LFC <- 1

counts <- read.delim(file.path(opt$cache, "gene_counts.tsv"), row.names = 1, check.names = FALSE)
meta   <- read.delim(file.path(opt$cache, "sample_meta.tsv"), stringsAsFactors = FALSE)
W      <- read.delim(file.path(opt$cache, "w_factors.tsv"), stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id
meta <- meta[colnames(counts), ]
stopifnot(!anyNA(meta$sample_id), setequal(W$sample_id, meta$sample_id))
meta$W_1 <- W$W_1[match(meta$sample_id, W$sample_id)]
meta$replicate <- factor(meta$replicate)
meta$condition <- factor(paste(meta$genotype, meta$stage, sep = "_"))
meta$RT <- as.numeric(meta$temperature == "RT")
mut <- setdiff(unique(meta$genotype), c("N2", "WT"))
meta$RT_mut <- meta$RT * (meta$genotype == mut)
adult <- meta$stage != "L4"
message("cohorts: ", paste(names(table(meta$temperature)), table(meta$temperature), collapse = ", "))

keep <- rowSums(counts) >= 10
counts <- as.matrix(counts[keep, ])

fit <- function(design) {
  dds <- DESeqDataSetFromMatrix(counts, meta, design)
  suppressMessages(DESeq(dds, quiet = TRUE))
}
deg <- function(res) !is.na(res$padj) & res$padj < PADJ & abs(res$log2FoldChange) >= LFC
geno <- function(dds, st) {
  ct <- c("condition", paste(mut, st, sep = "_"), paste(meta$genotype[meta$genotype != mut][1], st, sep = "_"))
  r <- results(dds, contrast = ct, alpha = PADJ)
  s <- lfcShrink(dds, contrast = ct, type = "ashr", quiet = TRUE)
  s$padj <- r$padj
  s
}
coef_res <- function(dds, coef) {
  r <- results(dds, name = coef, alpha = PADJ)
  s <- lfcShrink(dds, coef = coef, type = "ashr", quiet = TRUE)
  s$padj <- r$padj
  s
}

m0  <- fit(~ W_1 + replicate + condition)
m1  <- fit(~ W_1 + replicate + RT + condition)
m1b <- fit(~ replicate + RT + condition)
m2  <- fit(~ W_1 + replicate + RT + RT_mut + condition)

stages <- c("L4", "Day3", "Day6")
g0 <- lapply(setNames(stages, stages), function(st) geno(m0, st))
g1 <- lapply(setNames(stages, stages), function(st) geno(m1, st))
n0 <- sapply(g0, function(r) sum(deg(r)))
# the paper's Fig 1 counts are 29/84/251; this refit gives 28/84/253 (genes at
# the cutoffs; at Day 6 one log2FC is 1.0016 vs 0.9998 and one padj sits at the
# independent-filtering edge), so M1 is compared with this refit, not the files
message("paper model DEGs: ", paste(n0, collapse = "/"), " (paper: 29/84/251)")
stopifnot(all(abs(n0 - c(29, 84, 251)) <= 3))

cmp <- do.call(rbind, lapply(stages, function(st) {
  a <- rownames(g0[[st]])[deg(g0[[st]])]; b <- rownames(g1[[st]])[deg(g1[[st]])]
  data.frame(stage = st, DEGs_paper_model = length(a), DEGs_with_cohort = length(b),
             in_both = length(intersect(a, b)), only_paper_model = length(setdiff(a, b)),
             only_with_cohort = length(setdiff(b, a)),
             # genes with a padj in both models; near-zero-count genes have unstable log2FC
             log2FC_pearson_r_tested_genes = with(list(k = !is.na(g0[[st]]$padj) & !is.na(g1[[st]]$padj)),
               round(cor(g0[[st]]$log2FoldChange[k], g1[[st]]$log2FoldChange[k]), 4)),
             padj_of_paper_DEGs_median_ratio = round(median(g1[[st]][a, "padj"] / g0[[st]][a, "padj"], na.rm = TRUE), 3))
}))
write.table(cmp, file.path(opt$outdir, "genotype_degs_with_and_without_cohort.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cohort  <- coef_res(m1, "RT")
cohortb <- coef_res(m1b, "RT")
inter   <- coef_res(m2, "RT_mut")
write_degs <- function(r, file) {
  d <- as.data.frame(r)[deg(r), c("baseMean", "log2FoldChange", "padj")]
  d <- d[order(d$padj), ]
  write.csv(cbind(gene_id = rownames(d), d), file.path(opt$outdir, file), row.names = FALSE)
  nrow(d)
}
summary <- data.frame(
  test = c("cohort effect, room temperature vs 20 C (paper model + cohort)",
           "cohort effect without RUV (replicate + cohort + condition)",
           "genotype effect that differs between cohorts (cohort x hlh-30)"),
  DEGs = c(write_degs(cohort, "cohort_degs.csv"), write_degs(cohortb, "cohort_degs_no_ruv.csv"),
           write_degs(inter, "cohort_by_genotype_degs.csv")),
  genes_padj_below_0.05_any_size = c(sum(cohort$padj < PADJ, na.rm = TRUE), sum(cohortb$padj < PADJ, na.rm = TRUE),
                                     sum(inter$padj < PADJ, na.rm = TRUE)))
summary$W_1_vs_cohort_r_adults <- round(cor(meta$W_1[adult], meta$RT[adult]), 3)
write.table(summary, file.path(opt$outdir, "cohort_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(cmp); print(summary)

# PCA of the adults (top 500 variable genes), before and after removing W_1 and replicate
vsd <- assay(vst(m0, blind = FALSE))[, adult]
md  <- meta[adult, ]
corr <- removeBatchEffect(vsd, batch = md$replicate, covariates = md$W_1,
                          design = model.matrix(~ condition, droplevels(md)))
pca_df <- function(x, label) {
  top <- head(order(apply(x, 1, var), decreasing = TRUE), 500)
  p <- prcomp(t(x[top, ]))
  v <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
  data.frame(PC1 = p$x[, 1], PC2 = p$x[, 2], genotype = ifelse(md$genotype == mut, "hlh-30(tm1978)", "WT"),
             stage = sub("Day", "Day ", md$stage), cohort = ifelse(md$RT == 1, "22 \u00b0C (room temp.)", "20 \u00b0C"),
             panel = sprintf("%s\nPC1 %.1f%%, PC2 %.1f%%", label, v[1], v[2]))
}
d <- rbind(pca_df(vsd, "VST, no correction"), pca_df(corr, "VST, W_1 and replicate removed (as Fig. 1D)"))
g <- ggplot(d, aes(PC1, PC2, colour = cohort, shape = interaction(genotype, stage, sep = ", "))) +
  geom_point(size = 2.6, stroke = 0.8) +
  facet_wrap(~ panel, scales = "free") +
  scale_colour_manual(values = c("22 \u00b0C (room temp.)" = PAL$okabe[6], "20 \u00b0C" = PAL$okabe[5])) +
  scale_shape_manual(values = c(16, 1, 17, 2)) +
  labs(shape = "Genotype, stage", colour = "Cohort") +
  theme_bw(base_size = 10) + theme(panel.grid = element_blank(), legend.position = "bottom", legend.box = "vertical")
ggsave(file.path(opt$outdir, "pca_cohort.pdf"), g, width = 7.09, height = 4.2)
ggsave(file.path(opt$outdir, "pca_cohort.png"), g, width = 7.09, height = 4.2, dpi = 200)
