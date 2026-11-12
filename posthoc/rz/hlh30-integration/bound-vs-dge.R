#!/usr/bin/env Rscript
# HLH-30 binding (modENCODE L4 ChIP, via ChIP-Atlas) against the
# expression change in hlh-30(tm1978) vs WT, per stage.
#
#   bound gene  = a peak summit reproducible in both L4 replicates (summits
#                 within 100 bp) lies within -1000/+500 bp of the gene's TSS
#   ChIP signal = mean log2(ChIP / input) over TSS +/- 500 bp, averaged over
#                 the two L4 replicates (each normalised by its genome mean)
#   DEG         = padj < 0.05 and |log2FC| >= 1 (as in the volcanoes)
#
# Panels: A ChIP signal vs log2FC at L4 (the ChIP stage); B log2FC
# distribution (ECDF) of bound vs unbound genes per stage; C overlap of bound
# genes with up / down DEGs per stage (Fisher odds ratio, BH); D Venn-style
# counts at L4. Also writes the per-gene table.
#
#   Rscript bound-vs-dge.R           # -> out/bound/
#   Rscript bound-vs-dge.R --paper   # -> out/bound/paper/02a-bound-vs-deg, 02b-bound-share-pairwise

suppressPackageStartupMessages({
  library(optparse); library(rtracklayer); library(GenomicRanges); library(dplyr)
  library(tidyr); library(ggplot2); library(patchwork)
})
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "..", "_common", "palette.R"))
opt <- parse_args(OptionParser(option_list = list(
  make_option("--deg-dir", type = "character", default = "../volcano/cache/hlh30-all-ruv", dest = "deg_dir"),
  make_option("--chip",    type = "character", default = "../chip-tracks/cache"),
  make_option("--gtf",     type = "character", default = "../chip-tracks/cache/gtf.gtf"),
  make_option("--padj",    type = "double", default = 0.05),
  make_option("--lfc",     type = "double", default = 1),
  make_option("--up",      type = "integer", default = 1000, help = "bound window, bp upstream of TSS"),
  make_option("--down",    type = "integer", default = 500,  help = "bound window, bp downstream of TSS"),
  make_option("--outdir",  type = "character", default = "out/bound"),
  make_option("--paper",   action = "store_true", default = FALSE,
              help = "write only the paper panels (<outdir>/paper/02a, 02b) at 10 pt")
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
STAGES <- c(L1 = "L1", L4 = "L4", D3 = "Day3", D6 = "Day6")
to_ucsc <- function(x) { x <- as.character(x); ifelse(x == "MtDNA", "chrM", paste0("chr", x)) }

# ---- genes and TSSs ---------------------------------------------------------
gtf <- import(opt$gtf); genes <- gtf[gtf$type == "gene"]; genes <- genes[!duplicated(genes$gene_id)]
minus <- as.character(strand(genes)) == "-"
tss <- GRanges(to_ucsc(seqnames(genes)), IRanges(ifelse(minus, end(genes), start(genes)), width = 1),
               strand = strand(genes), gene_id = genes$gene_id, gene_name = genes$gene_name)

# ---- reproducible L4 summits -> bound genes --------------------------------
read_np <- function(srx) {
  b <- read.table(file.path(opt$chip, paste0(srx, ".bed")), sep = "\t", stringsAsFactors = FALSE)
  GRanges(b$V1, IRanges(b$V2 + b$V10 + 1, width = 1), score = b$V5)
}
r1 <- read_np("SRX331310"); r2 <- read_np("SRX331312")
h <- distanceToNearest(r1, r2); summits <- r1[queryHits(h)[which(mcols(h)$distance <= 100)]]
win <- promoters(tss, upstream = opt$up, downstream = opt$down + 1)
bound <- unique(win$gene_id[queryHits(findOverlaps(win, summits))])
message(length(summits), " reproducible L4 summits; ", length(bound), " genes bound (-", opt$up, "/+", opt$down, ")")

# ---- ChIP signal at the TSS --------------------------------------------------
genome_mean <- function(bw) {
  si <- seqinfo(bw); gr <- GRanges(seqnames(si), IRanges(1, seqlengths(si)))
  s <- summary(bw, which = gr, size = 1, type = "mean", defaultValue = 0)
  v <- vapply(s, function(x) as.numeric(score(x))[1], numeric(1))
  sum(v * seqlengths(si), na.rm = TRUE) / sum(as.numeric(seqlengths(si)))
}
tss_mean <- function(srx, gr) {
  bw <- BigWigFile(file.path(opt$chip, paste0(srx, ".bw")))
  gr <- gr[as.character(seqnames(gr)) %in% seqnames(seqinfo(bw))]
  sl <- seqlengths(seqinfo(bw))[as.character(seqnames(gr))]
  start(gr) <- pmax(1, start(gr)); end(gr) <- pmin(sl, end(gr))
  v <- vapply(summary(bw, which = gr, size = 1, type = "mean", defaultValue = 0),
              function(x) as.numeric(score(x))[1], numeric(1))
  setNames(v / genome_mean(bw), gr$gene_id)
}
tw <- resize(tss, 1001, fix = "center"); strand(tw) <- "*"
pairs <- list(c("SRX331310", "SRX331309"), c("SRX331312", "SRX331311"))
lr <- sapply(pairs, function(p) {
  a <- tss_mean(p[1], tw); b <- tss_mean(p[2], tw)
  log2((a[names(a)] + 0.05) / (b[names(a)] + 0.05))
})
chip_sig <- data.frame(gene_id = rownames(lr), chip = rowMeans(lr), stringsAsFactors = FALSE)

# ---- expression -------------------------------------------------------------
dge <- bind_rows(lapply(names(STAGES), function(st) {
  d <- read.csv(file.path(opt$deg_dir, sprintf("tm1978_%s_vs_WT_%s.csv", STAGES[[st]], STAGES[[st]])),
                stringsAsFactors = FALSE)
  d %>% filter(!is.na(padj)) %>%
    transmute(gene_id, gene_name, stage = st, lfc = log2FoldChange, padj,
              deg = case_when(padj < opt$padj & log2FoldChange >=  opt$lfc ~ "up",
                              padj < opt$padj & log2FoldChange <= -opt$lfc ~ "down",
                              TRUE ~ "ns"))
})) %>% mutate(stage = factor(stage, levels = names(STAGES)),
               bound = gene_id %in% bound) %>%
  left_join(chip_sig, by = "gene_id")
write.table(dge, file.path(opt$outdir, "gene_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- statistics --------------------------------------------------------------
fish <- dge %>% group_by(stage) %>% group_modify(function(d, k) {
  bind_rows(lapply(c("up", "down"), function(dir) {
    tb <- table(factor(d$bound, c(TRUE, FALSE)), factor(d$deg == dir, c(TRUE, FALSE)))
    ft <- fisher.test(tb)
    data.frame(direction = dir, n_deg = sum(d$deg == dir), n_bound_deg = tb[1, 1],
               pct_bound_deg = 100 * tb[1, 1] / max(1, sum(d$deg == dir)),
               pct_bound_all = 100 * mean(d$bound),
               # the Fisher reference group: every other tested gene at this stage
               pct_bound_other = 100 * tb[1, 2] / max(1, sum(tb[, 2])),
               or = unname(ft$estimate), lo = ft$conf.int[1], hi = ft$conf.int[2], p = ft$p.value)
  }))
}) %>% ungroup() %>% mutate(q = p.adjust(p, "BH"))
wil <- dge %>% group_by(stage) %>%
  summarise(median_bound = median(lfc[bound]), median_unbound = median(lfc[!bound]),
            p = wilcox.test(lfc[bound], lfc[!bound])$p.value, .groups = "drop") %>%
  mutate(q = p.adjust(p, "BH"))
write.table(fish, file.path(opt$outdir, "fisher_bound_vs_deg.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wil,  file.path(opt$outdir, "wilcoxon_bound_lfc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(as.data.frame(fish)); print(as.data.frame(wil))

# ---- figures ------------------------------------------------------------------
# Wilson 95% CI for the bound share among each DEG set
wilson <- function(x, n) {
  if (n == 0) return(c(NA, NA)); z <- 1.96; ph <- x / n
  c((ph + z^2 / (2 * n) - z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / (1 + z^2 / n),
    (ph + z^2 / (2 * n) + z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)) * 100
}
fc <- fish %>% rowwise() %>% mutate(ci_lo = wilson(n_bound_deg, n_deg)[1], ci_hi = wilson(n_bound_deg, n_deg)[2]) %>%
  ungroup() %>% mutate(row = paste(stage, direction),
                      row = factor(row, levels = rev(paste(rep(names(STAGES), each = 2), c("down", "up")))),
                      or_p = pmax(or, 0.02), lo_p = pmax(lo, 0.02), hi_p = pmin(hi, 50),
                      lab = sprintf("%d/%d bound (%.0f%%)  q = %.1g", n_bound_deg, n_deg, pct_bound_deg, q))
# ---- pairwise: is the bound share different between two DEG sets? ------------
# Each stage's whole DEG list is compared (the question: does the bound share
# of one stage's DEGs differ from another's?). The lists overlap: a gene DE at
# both stages sits on both sides, so a 2x2 Fisher test (which assumes two
# independent groups) is wrong here. Instead an exact conditional test: under
# the null the bound labels are exchangeable over the genes in either list, so
# given which list each gene is in (A only, both, B only) and the number of
# bound genes, the bound counts in the three groups are multivariate
# hypergeometric. p = total probability of a difference in % bound at least
# as large as observed, both sides, by full enumeration (the same as shuffling
# the labels, without Monte Carlo error). BH runs within each direction (6
# pairs), since up and down are separate questions. Overall test per
# direction: shuffled labels with the chi-squared statistic over all stages.
# The Fisher test on stage-specific genes only (shared genes removed) is kept
# in the TSV as a check.
N_PERM <- 20000
set.seed(1)
deg_sets <- dge %>% filter(deg != "ns") %>% select(gene_id, stage, direction = deg, bound)
perm_p <- function(obs, null) (1 + sum(abs(null) >= abs(obs) - 1e-12)) / (1 + length(null))
exact_p <- function(nAo, nb, nBo, K, dobs) {
  nA <- nAo + nb; nB <- nb + nBo; lden <- lchoose(nAo + nb + nBo, K); tot <- 0
  for (a in 0:min(nAo, K)) for (b in 0:min(nb, K - a)) {
    c <- K - a - b
    if (c <= nBo && abs((a + b) / nA - (b + c) / nB) >= abs(dobs) - 1e-12)
      tot <- tot + exp(lchoose(nAo, a) + lchoose(nb, b) + lchoose(nBo, c) - lden)
  }
  min(1, tot)
}
pairs <- bind_rows(lapply(c("down", "up"), function(dir) {
  st <- names(STAGES)
  bind_rows(lapply(combn(st, 2, simplify = FALSE), function(pr) {
    a  <- deg_sets %>% filter(direction == dir, stage == pr[1]); b <- deg_sets %>% filter(direction == dir, stage == pr[2])
    sh <- intersect(a$gene_id, b$gene_id)
    a1 <- a %>% filter(!gene_id %in% sh); b1 <- b %>% filter(!gene_id %in% sh)
    k_sh <- sum(a$bound[a$gene_id %in% sh])
    tb <- matrix(c(sum(a1$bound), sum(!a1$bound), sum(b1$bound), sum(!b1$bound)), 2)
    data.frame(direction = dir, s1 = pr[1], s2 = pr[2], shared = length(sh), n1 = nrow(a), n2 = nrow(b),
               pct1 = 100 * mean(a$bound), pct2 = 100 * mean(b$bound),
               p = exact_p(nrow(a1), length(sh), nrow(b1), sum(a1$bound) + k_sh + sum(b1$bound),
                           mean(a$bound) - mean(b$bound)),
               n1_unique = nrow(a1), n2_unique = nrow(b1),
               pct1_unique = 100 * mean(a1$bound), pct2_unique = 100 * mean(b1$bound),
               p_unique = if (nrow(a1) && nrow(b1)) fisher.test(tb)$p.value else NA)
  }))
})) %>% group_by(direction) %>%
  mutate(q = p.adjust(p, "BH"), q_unique = p.adjust(p_unique, "BH")) %>% ungroup() %>%
  mutate(diff = pct1 - pct2)
# chi-squared statistic of the stage x bound table, genes counted once per stage
chi_stat <- function(M, b) {
  k <- colSums(M * b); n <- colSums(M); e <- n * sum(k) / sum(n)
  sum((k - e)^2 / e + ((n - k) - (n - e))^2 / (n - e))
}
het <- bind_rows(lapply(c("down", "up"), function(dir) {
  s <- deg_sets %>% filter(direction == dir); g <- unique(s$gene_id)
  M <- sapply(names(STAGES), function(st) g %in% s$gene_id[s$stage == st])
  b <- s$bound[match(g, s$gene_id)]
  data.frame(direction = dir, p = perm_p(chi_stat(M, b), replicate(N_PERM, chi_stat(M, sample(b)))))
}))
write.table(pairs, file.path(opt$outdir, "pairwise_bound_share.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(as.data.frame(pairs)); print(as.data.frame(het))

# One reference line for readability. Each row is tested against its own
# stage's other genes, but those shares only span 14-17%, so the mean stands in.
base_pct <- mean(fish$pct_bound_other)
base_rng <- range(fish$pct_bound_other)
at_floor <- function(x) x <= 1.0001 / (N_PERM + 1)

# ---- paper panels: one file each, every label at 10 pt -----------------------
# The draft's titles, subtitles and caption move to the figure legend; the
# q values become stars and the n per side stays in pairwise_bound_share.tsv.
if (opt$paper) {
  PT <- 10; GT <- PT / .pt
  PA_W <- 3.75; PA_H <- 2.6; PB_W <- 2.1; PB_H <- 3.3; PB_ASPECT <- 0.42
  stars <- function(q) case_when(q < 0.001 ~ "***", q < 0.01 ~ "**", q < 0.05 ~ "*", TRUE ~ "")
  th <- theme(text = element_text(size = PT), axis.text = element_text(size = PT, colour = "grey20"),
              axis.title = element_text(size = PT), legend.text = element_text(size = PT),
              legend.title = element_text(size = PT), strip.text = element_text(size = PT, face = "bold"),
              plot.tag = element_text(size = PT, face = "bold"))
  pa <- ggplot(fc, aes(pct_bound_deg, row, colour = direction)) +
    geom_vline(xintercept = base_pct, colour = "grey45", linetype = "22") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.3) +
    geom_point(size = 1.8) +
    # n/N inside the panel, just right of each interval, at 8 pt; a label
    # that would cross the dashed line starts just right of it instead
    geom_text(aes(x = ifelse(ci_hi < base_pct & ci_hi + 13.5 > base_pct, base_pct + 1.5, ci_hi + 1.5), label = sprintf("%d/%d%s", n_bound_deg, n_deg, stars(q))),
              hjust = 0, size = 8 / .pt, colour = "grey20") +
    scale_x_continuous(breaks = seq(0, 60, 20), labels = function(v) paste0(v, "%"),
                       expand = expansion(mult = c(0.02, 0))) +
    coord_cartesian(xlim = c(0, 75)) +
    scale_colour_manual(values = c(up = PAL$up, down = PAL$down), breaks = c("up", "down"), name = "DEGs") +
    labs(x = "DEGs with an HLH-30-bound promoter", y = NULL) +
    theme_bw(base_size = PT) + th +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          # the row labels already say up or down
          legend.position = "none", plot.margin = margin(2, 6, 2, 2))
  # each stage's bound % (as in 02a) under its label, so a cell can be read
  # off as column % minus row %; the % differs per direction, so the axes are free
  pct_of <- setNames(round(fc$pct_bound_deg), paste(fc$stage, fc$direction))
  ax <- function(st, dir) paste0(st, "\n", pct_of[paste(st, dir)], "%")
  ax_levels <- function(sts) unique(unlist(lapply(c("down", "up"), function(d) ax(sts, d))))
  pq <- pairs %>% mutate(
    xs = factor(ax(s1, direction), levels = ax_levels(names(STAGES))),
    ys = factor(ax(s2, direction), levels = ax_levels(rev(names(STAGES)))),
    pfmt = with(het[match(direction, het$direction), ], ifelse(at_floor(p), "< 1e-4", sprintf("= %.1g", p))),
    # colour = the difference, stars = q (as in the SEA grid); no numbers in the cells
    lab = ifelse(is.na(q), "n.t.", stars(q)))
  heat <- function(d, strip_fmt, ...) {
    ggplot(d %>% mutate(dir_lab = sprintf(strip_fmt, direction, pfmt)),
           aes(xs, ys, fill = ifelse(q < 0.05, pmax(pmin(diff, 30), -30), NA))) +
      geom_tile(colour = "white", linewidth = 1) +
      # an asterisk sits at the top of its line box, so drop it to the cell middle
      geom_text(aes(label = lab), size = GT, colour = "grey10", vjust = 0.78) +
      facet_wrap(~ dir_lab, scales = "free", ...) +
      scale_fill_gradient2(low = PAL$div[["low"]], mid = "white", high = PAL$div[["high"]], midpoint = 0,
                           limits = c(-30, 30), breaks = c(-30, 0, 30), na.value = "grey95",
                           name = "difference\n(% points)") +
      scale_x_discrete(position = "top", expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_size = PT) + th +
      theme(panel.grid = element_blank(), legend.key.height = unit(14, "pt"), legend.key.width = unit(8, "pt"),
            legend.title = element_text(size = PT, margin = margin(b = 6)), strip.placement = "outside",
            axis.text.y = element_text(size = PT, colour = "grey20", hjust = 0.5, lineheight = 0.9),
            axis.text.x.top = element_text(size = PT, colour = "grey20", lineheight = 0.9, margin = margin(b = 2)),
            strip.clip = "off", aspect.ratio = PB_ASPECT, panel.spacing = unit(10, "pt"),
            plot.margin = margin(6, 4, 2, 2))
  }
  # one third of the page: facets stacked, short strip titles, row labels on
  # one line so each row is one text line tall, colour bar below
  pq <- pq %>% mutate(ys = factor(sub("\n", " ", ys), levels = sub("\n", " ", levels(ys))))
  pb <- heat(pq, "%s DEGs, p %s", ncol = 1) +
    guides(fill = guide_colourbar(title = "difference (% points)", title.position = "top", title.hjust = 0.5,
                                  theme = theme(legend.key.width = unit(1.1, "in"), legend.key.height = unit(6, "pt")))) +
    theme(aspect.ratio = PB_ASPECT, legend.position = "bottom", legend.direction = "horizontal",
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(8, "pt"),
          plot.margin = margin(4, 8, 2, 2), axis.text.y = element_text(size = PT, colour = "grey20", hjust = 1))
  pdir <- file.path(opt$outdir, "paper"); dir.create(pdir, showWarnings = FALSE)
  unlink(file.path(pdir, paste0(rep(c("02-bound-vs-deg", "02b-bound-share-pairwise-half"), each = 2), c(".pdf", ".png"))))
  save2 <- function(p, name, w, h) {
    ggsave(file.path(pdir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
    # the PNG is rasterised from the PDF: the png() device spaces the axis text differently
    system2("pdftoppm", c("-r", "300", "-png", "-singlefile", shQuote(file.path(pdir, paste0(name, ".pdf"))),
                          shQuote(file.path(pdir, name))))
    message("wrote ", pdir, "/", name, ".{pdf,png} (", w, " x ", h, " in)")
  }
  save2(pa, "02a-bound-vs-deg", PA_W, PA_H)
  save2(pb, "02b-bound-share-pairwise", PB_W, PB_H)
  quit(save = "no")
}

pC <- ggplot(fc, aes(pct_bound_deg, row, colour = direction)) +
  geom_vline(xintercept = base_pct, colour = "grey45", linetype = "22") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25) +
  geom_point(size = 2.2) +
  geom_text(aes(x = 64, label = sprintf("%d/%d   q = %.1g", n_bound_deg, n_deg, q)),
            hjust = 0, size = 2.8, colour = "grey25") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 60, 20), labels = function(v) paste0(v, "%")) +
  coord_cartesian(xlim = c(0, 60), clip = "off") +
  scale_colour_manual(values = c(up = PAL$up, down = PAL$down), name = "DEGs") +
  labs(x = sprintf("%% of the DEGs whose promoter is bound by HLH-30 (95%% CI)\ndashed line = %.0f%%, all other expressed genes (%.0f-%.0f%% by stage)",
                   base_pct, base_rng[1], base_rng[2]),
       y = NULL,
       title = "A  Are HLH-30-bound genes over-represented among DEGs?",
       subtitle = "right: bound DEGs / DEGs; q: Fisher's exact test vs all other expressed genes, BH across the 8 tests") +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
                                   legend.position = "bottom", plot.margin = margin(5, 95, 5, 5))

# the overall permutation p cannot go below 1 / (N_PERM + 1); no shuffle
# reached the observed value there, so show it as an upper bound (at_floor)
pp <-pairs %>% mutate(
  s1 = factor(s1, levels = names(STAGES)), s2 = factor(s2, levels = rev(names(STAGES))),
  dir_lab = sprintf("%s DEGs (overall across stages: p %s)", direction,
                    with(het[match(direction, het$direction), ], ifelse(at_floor(p), "< 1e-4", sprintf("= %.1g", p)))),
  star = case_when(q < 0.001 ~ "***", q < 0.01 ~ "**", q < 0.05 ~ "*", TRUE ~ ""),
  lab = ifelse(is.na(q), sprintf("n = %d vs %d\nno test", n1, n2),
               sprintf("%+.0f pts%s\nq = %.1g\nn = %d vs %d", diff, star, q, n1, n2)))
pB <- ggplot(pp, aes(s1, s2, fill = ifelse(q < 0.05, pmax(pmin(diff, 30), -30), NA))) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = lab), size = 2.7, lineheight = 0.95) +
  facet_wrap(~ dir_lab, nrow = 1) +
  scale_fill_gradient2(low = PAL$div[["low"]], mid = "white", high = PAL$div[["high"]], midpoint = 0,
                       limits = c(-30, 30), na.value = "grey95",
                       name = "column stage\nminus row stage\n(% points, q < 0.05)") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "B  Does the bound share differ between stages?",
       subtitle = paste0("each cell compares the whole DEG lists of two stages (same direction): % bound in the column stage ",
                         "minus % bound in the row stage;\nn = DEGs per side; exact conditional test for lists that share genes, ",
                         "BH over the 6 pairs per direction; grey = not significant")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"), strip.placement = "outside")

cap <- sprintf(paste0("Bound = a summit reproducible in both L4 HLH-30 ChIP replicates (SRX331310, SRX331312; within 100 bp) ",
                      "lies within -%d/+%d bp of the TSS: %d genes from %d summits. DEG: padj < %.2g, |log2FC| >= %g. ",
                      "Panel B: a gene DE at both stages of a pair counts on both sides. The exact test enumerates every way to ",
                      "place the bound genes over the genes of the two lists (A only, both, B only), so it keeps this overlap, ",
                      "which a 2x2 Fisher test cannot. Overall p: %d label shuffles with the chi-squared statistic over all stages."),
               opt$up, opt$down, length(bound), length(summits), opt$padj, opt$lfc, N_PERM)
# references for the statistical choices, printed under the panels
refs <- c(
  "Statistics. A: Fisher's exact test vs all other expressed genes; 95% CI by the Wilson score interval (Wilson 1927).",
  paste("B: overlapping DEG lists are partially overlapping samples, which a test for independent groups must not",
        "treat as unpaired (Derrick et al. 2015); the exact conditional test is the permutation test enumerated in full",
        "(Ernst 2004). A and B: Benjamini-Hochberg FDR (Benjamini & Hochberg 1995) over the tests that support one",
        "conclusion (up and down separately in B; Bender & Lange 2001)."),
  paste("Benjamini Y, Hochberg Y (1995) J R Stat Soc B 57:289-300.",
        "Bender R, Lange S (2001) J Clin Epidemiol 54:343-349.",
        "Derrick B, Dobson-McKittrick A, Toher D, White P (2015) J Appl Quant Methods 10(3):1-14.",
        "Ernst MD (2004) Stat Sci 19:676-685.",
        "Wilson EB (1927) J Am Stat Assoc 22:209-212."))
fig <- pC / pB + plot_layout(heights = c(1.1, 0.9)) +
  plot_annotation(title = "HLH-30 binding vs expression change in hlh-30(tm1978) (draft for discussion)",
                  caption = paste(c(stringr::str_wrap(cap, 160), stringr::str_wrap(refs, 160)), collapse = "\n"),
                  theme = theme(plot.title = element_text(face = "bold", size = 13),
                                plot.caption = element_text(colour = "grey40", hjust = 0, size = 8)))
ggsave(file.path(opt$outdir, "bound_vs_dge.pdf"), fig, width = 10, height = 8.6, device = cairo_pdf)
ggsave(file.path(opt$outdir, "bound_vs_dge.png"), fig, width = 10, height = 8.6, dpi = 150)
message("wrote ", opt$outdir, "/bound_vs_dge.{pdf,png}")
